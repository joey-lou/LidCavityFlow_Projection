% test against benchmark
clear; clc; close all;

clear convert_UV2W

% load('Data/Re5K_N32.mat')
load('Data/Re1K_N24.mat')

Grid = build_Grid(Parameter.N);
Operators = build_Operators(Grid,Parameter);
Default_V = 16*Grid.x.^2.*(1-Grid.x).^2;


m = (Parameter.N+1)^2;
UV = Z(:,end);
W = convert_UV2W(UV,Operators.DX,Operators.DY,m);
w = W;
w(Grid.bd_pts) = 0;
PSI = -Operators.oL*w;


% plot average velocity profile at steady state
figure
U = reshape(UV(1:m),Parameter.N+1,Parameter.N+1);
V = reshape(UV(m+1:2*m),Parameter.N+1,Parameter.N+1);
avg_U = mean(U,2);
avg_V = mean(V,1);

plot(avg_U,Grid.y*2-1,'-*')
hold on
plot(Grid.x*2-1,avg_V,'-*')
legend('U','V')
title(['Average Velocity Profile at Re = ',num2str(Parameter.Re)])


figure
subplot(2,2,1)
Grid.cplot(PSI)
subplot(2,2,2)
Grid.cplot(W)
subplot(2,2,3)
Grid.mplot(UV(1:m))
subplot(2,2,4)
Grid.mplot(UV(m+1:end))


fprintf('Top Boundary Error: %5.3e\n',norm(UV(Grid.t_pts)-Default_V))
fprintf('U Error: %5.3e\n',norm(UV([Grid.b_pts,Grid.l_pts,Grid.r_pts])))
fprintf('V Error: %5.3e\n',norm(UV(Grid.bd_pts+Grid.m)))

XX = reshape(Grid.xx,Parameter.N+1,Parameter.N+1);
YY = reshape(Grid.yy,Parameter.N+1,Parameter.N+1);
% 
% fprintf('Top Boundary Error: %5.3e\n',norm(UV(Grid.t_pts)-Default_V))
% fprintf('U Error: %5.3e\n',norm(UV([Grid.b_pts,Grid.l_pts,Grid.r_pts])))
% fprintf('V Error: %5.3e\n',norm(UV(Grid.bd_pts+Grid.m)))
% fprintf('M1 = %6.5f at (%5.4f, %5.4f)\n',M1/2,(XX(row,col)+1)/2,(YY(row,col)+1)/2)
% fprintf('W = %6.5f\n', Omega(row,col))
% fprintf('M2 = %6.5f at (%5.4f, %5.4f)\n',M2*2,(Grid.x(x_2)+1)/2,1)




fprintf('Benchmark Comparison at Re = %d\n', Parameter.Re);
P = min(PSI(Grid.i_pts));
[row.P, col.P] = find(reshape(PSI,Parameter.N+1,Parameter.N+1)==P);
fprintf('Primary Vortex: %7.4f at (%5.3f, %5.3f)\n',P,XX(row.P,col.P),YY(row.P,col.P))

S.br = max(PSI(Grid.br_pts));
[row.br, col.br] = find(reshape(PSI,Parameter.N+1,Parameter.N+1)==S.br);
fprintf('Secondary Vortex (bottom right): %6.4E at (%5.3f, %5.3f)\n',S.br,XX(row.br,col.br),YY(row.br,col.br))
S.bl = max(PSI(Grid.bl_pts));
[row.bl, col.bl] = find(reshape(PSI,Parameter.N+1,Parameter.N+1)==S.bl);
fprintf('Secondary Vortex (bottom left):  %6.4E at (%5.3f, %5.3f)\n',S.bl,XX(row.bl,col.bl),YY(row.bl,col.bl))

if Parameter.Re>=2000
    S.ul = max(PSI(Grid.ul_pts));
    [row.ul, col.ul] = find(reshape(PSI,Parameter.N+1,Parameter.N+1)==S.ul);
    fprintf('Secondary Vortex (bottom left):  %6.4E at (%5.3f, %5.3f)\n',S.ul,XX(row.ul,col.ul),YY(row.ul,col.ul))
end


%% Plotting PSI Evolution
% figure
% for i = 1:floor(size(Z,2)/20):size(Z,2)
%     Grid.cplot(Z(:,i));
%     pause(0.1)
% end