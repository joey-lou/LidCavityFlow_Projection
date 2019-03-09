% Generates cavity flow stream function data for DMD purposes
% 
% Set parameters
clear; close all; clc;
clear nonLinearSolver
rng(1)
% load('Data/Re1K_N24.mat') % need improvement

savefilename = 'Data/Re20K_N69.mat';
Parameter.dt = 0.005;     % simulation time step dt
Parameter.Re = 20000;
Parameter.N  = 69;
Parameter.DT = 0.04;      % store state at DT
Parameter.T  = 40;       % final T for each run, # steps = T/DT
Parameter.runs = 1;      % number of runs, each of length T
umax = 1.0; umin = 1.0;  % range for top lid velocity
U = umin + (umax-umin)*rand(Parameter.runs,Parameter.T/Parameter.DT);
[Grid] = build_Grid(Parameter.N);
Default_V = 16*Grid.x.^2.*(1-Grid.x).^2; % specify default top lid velocity
% Default_V = ones(length(Grid.x),1); 
% STATES = Z(:,end);
X = [];
Y = [];
Z = [];

% vector recording actual steps taken in each run
SimLength = ones(1,Parameter.runs)*(Parameter.T/Parameter.DT); 
for i = 1:Parameter.runs
    % start with randomized initial conditions
    b = rand;
    a = [b,1-b]+0.1*randn(1,2);
    fprintf('a = %4.3f, %4.3f\n',a);
%     STATES = a(1)*LimitCycle_Re13k(:,randi(50)) + a(2)*FixedPoint_Re10k;
    STATES = zeros(2*(Parameter.N+1)^2,1);
%     fprintf('Running %i/%i simulation...', i, Parameter.runs)
%     tic
    steps = Parameter.T/Parameter.DT;
    for j = 1:steps
        fprintf('Running %i-%i/%i simulation...', i, j, steps)
        tic
        [STATE, dSTATE] = nonLinearSolver(Parameter,STATES(:,end),Default_V*(U(i,j)));
        if ~isempty(find(isnan(STATE),1))
            fprintf('Nonconvergence in solver occurred at step %i\n',j)
            SimLength(i)=j-1;
            break
        end
        if dSTATE < 1e-8
            fprintf('Solution converged at step %i\n',j)
            SimLength(i)=j-1;
            break
        end
        STATES = [STATES STATE];
        toc
    end
    
    X = [X STATES(:,1:end-1)];
    Y = [Y STATES(:,2:end)];
    Z = [Z STATES];
%     toc
end
% saving for DMD
save(savefilename,'X','Y','Z','Parameter','SimLength')
% save(savefilename,'X','Parameter','SimLength')
