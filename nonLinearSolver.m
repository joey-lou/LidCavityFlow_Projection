function [STATE, dSTATE] = nonLinearSolver(Parameter,STATE,INPUT)
% inputs:
% INPUT     - vector of velocity profile at top boundary
% STATE     - the stream function that starts this period's simulation from
% Parameter - .N spectral order
%           - .dt simulation time step
%           - .Re base reynolds number
%           - .DT simulation period time (DT/dt = steps taken) 
%              - .T final time horizon for simulation
% outputs:
% STATE:    - vector that represents UV at last time step during the period

persistent Grid Operators LHS

if isempty(Grid)
%     fprintf('Building Grid...\n')
%     tic
    [Grid] = build_Grid(Parameter.N);
%     toc
end
if isempty(Operators)
%     fprintf('Building Operators...\n')
%     tic
    [Operators] = build_Operators(Grid, Parameter);
%     toc
end

%% Normalize inputs
NominalInput= 16*Grid.x.^2.*(1-Grid.x).^2;
VelocityRatio = norm(INPUT,inf)/norm(NominalInput,inf);
U_lid = INPUT./VelocityRatio;
Re = Parameter.Re*VelocityRatio; % adjusted reynolds


%% LHS
if isempty(LHS) % assuming top lid velocity constant
    LHS = Grid.II-0.5*Parameter.dt*Operators.DEL2/Re;
    % clear rows and apply bcs
    LHS([Grid.bd_pts; Grid.bd_pts+Grid.m],:) = 0; 
    LHS([Grid.bd_pts; Grid.bd_pts+Grid.m],[Grid.bd_pts; Grid.bd_pts+Grid.m]) = eye(8*Parameter.N);
    LHS = inv(LHS);
end

%% initialize variables
persistent UV Pressure Adv
if isempty(UV)
    UV.current = STATE;
    UV.previous = STATE;
end

if isempty(Pressure)
    Pressure.current = cal_pressure(UV.current,Operators,Parameter,Grid);
    Pressure.previous = cal_pressure(UV.previous,Operators,Parameter,Grid);
end

if isempty(Adv)
    Adv.current = cal_adv(Operators.DX,Operators.DY,UV.current,Grid.m);
    Adv.previous = cal_adv(Operators.DX,Operators.DY,UV.previous,Grid.m);
end

%% iterate
for i = 1:(Parameter.DT/Parameter.dt)
    % RHS
    RHS = ...
        (Grid.II + ...
        0.5*Parameter.dt/Re*Operators.DEL2)*UV.current - ...
        1.5*Parameter.dt*Adv.current + ...
        0.5*Parameter.dt*Adv.previous;

    % boundary conditions
    temp = 2*Pressure.current-Pressure.previous;
    RHS([Grid.bd_pts; Grid.bd_pts+Grid.m]) = 0;                    
    RHS(Grid.t_pts) = U_lid + ...
        Parameter.dt*Operators.DX(Grid.t_pts,:)*temp; 
    RHS(Grid.b_pts) = ...
        Parameter.dt*Operators.DX(Grid.b_pts,:)*temp; 
    RHS([Grid.l_pts; Grid.r_pts]+Grid.m) = ...
        Parameter.dt*Operators.DY([Grid.l_pts; Grid.r_pts],:)*temp; 
    % normal direction
%     RHS([Grid.l_pts; Grid.r_pts]) = ...
%         Parameter.dt*Operators.DX([Grid.l_pts; Grid.r_pts],:)*temp; 
%     RHS([Grid.t_pts; Grid.b_pts]+Grid.m) = ...
%         Parameter.dt*Operators.DY([Grid.t_pts; Grid.b_pts],:)*temp; 

    UV.previous = UV.current;
    UV.temp = LHS*RHS;
    Pressure.previous = Pressure.current;
    Pressure.current = cal_pressure(UV.temp,Operators,Parameter,Grid);
    UV.current = UV.temp - ...
        Parameter.dt*[Operators.DX*Pressure.current;Operators.DY*Pressure.current];
    Adv.previous = Adv.current;
    Adv.current = cal_adv(Operators.DX,Operators.DY,UV.current,Grid.m);
    
end

STATE = UV.current; % update STATE
%% plot and analysis
W.current = convert_UV2W(UV.current,Operators.DX,Operators.DY,Grid.m);
W.previous = convert_UV2W(UV.previous,Operators.DX,Operators.DY,Grid.m);
dSTATE = norm(W.current-W.previous,2)/norm(W.current,2)/Parameter.dt;
Div = Operators.DX*UV.current(1:Grid.m)+Operators.DY*UV.current(Grid.m+1:end);
DIV = sum(sum(Div(Grid.i_pts).^2))/Parameter.N^2;
fprintf('Change in W: %3.1e, Divergence: %3.1e\n',dSTATE,DIV)
% 
% P_error = Operators.DX*UV.temp(1:Grid.m)+Operators.DY*UV.temp(Grid.m+1:end)-Parameter.dt*Operators.L*Pressure.current;
% nP = Operators.DX([Grid.l_pts; Grid.r_pts],:)*Pressure.current+...
%     Operators.DY([Grid.t_pts; Grid.b_pts],:)*Pressure.current;
% tP = Operators.DY([Grid.l_pts; Grid.r_pts],:)*Pressure.current+...
%     Operators.DX([Grid.t_pts; Grid.b_pts],:)*Pressure.current;
% fprintf('Normal Pressure Boundary Value: %3.1e\n',norm(nP))
% fprintf('Tangent Pressure Boundary Value: %3.1e\n',norm(tP))
% fprintf('Pressure Error: %3.1e\n',norm(P_error))
% subplot(2,2,1)
% Grid.mplot(STATE(1:Grid.m))
% subplot(2,2,2)
% Grid.mplot(STATE(Grid.m+1:end))
% subplot(2,2,3)
% Grid.mplot(Div)
% subplot(2,2,4)
% Grid.cplot(Div)
% pause()

% Grid.mplot(UV(1:Grid.m)-uv_psi(1:Grid.m))
% pause(0.01)



end

%% Helper functions
function Adv = cal_adv(DX,DY,UV,m)
% calculate advection term
u = UV(1:m);
v = UV(m+1:2*m);
adv_1 = u.*(DX*u)+v.*(DY*u);
adv_2 = u.*(DX*v)+v.*(DY*v);
Adv = [adv_1;adv_2];
end

function Pressure = cal_pressure(UV,Operators,Parameter,Grid)
% update pressure term
RHS = (Operators.DX*UV(1:Grid.m)+...
    Operators.DY*UV(Grid.m+1:2*Grid.m))/Parameter.dt;
% temp = RHS;
RHS(Grid.bd_pts) = 0;
Pressure = Operators.pL*RHS;
% Pressure = Operators.oL*RHS;
% fprintf('Pressure error: %3.1e\n',norm(Operators.L*Pressure-temp))
end