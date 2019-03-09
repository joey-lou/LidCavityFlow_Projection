% POD and DMD analysis
clear; clc; close all;
load('Data/Re20K_N69_R1_DT08_T80.mat')
load('Data/Re20K_N69_R1_DT4_T200.mat')
% clear convert_PSI2UV convert_Hankel
% load('Data/69_Train_20K.mat')
% load('Data/49_Train_10K.mat')
% load('Data/49_Test_10K_1.mat')
% load('Data/Re20K_N69_R1_T40.mat')
% load('Data/69_Test_20K_1.mat')

X = X - mean(X,2);
Y = Y - mean(Y,2);
Z = Z - mean(Z,2);
Grid = build_Grid(Parameter.N);

%% DMD
fprintf('Analyzing DMD modes...\n')
tol = 0.9999;
[A_s, U, Phi, eigs] = DMD(X,Y,tol);
% Plot DMD modes
N = size(X,1)/2;
% for i = 1:9
%     figure(1)
%     subplot(3,3,i)
%     Grid.cplot(real(Phi.exact(1:N,i*2)))
%     title(['U DMD mode ',num2str(i)])
%     figure(2)
%     subplot(3,3,i)
%     Grid.cplot(real(Phi.exact((N+1):end,i*2)))
%     title(['V DMD mode ',num2str(i)])
% end
% Plot DMD spectrum
figure(3)
theta = (0:1:100)*2*pi/100;
plot(cos(theta),sin(theta),'k--') % plot unit circle
hold on, grid on
scatter(real(diag(eigs)),imag(diag(eigs)),'ok')
axis([-1.1 1.1 -1.1 1.1]);
title 'Koopman Eigenvalues'

%% POD
fprintf('Analyzing POD modes...\n')
[PSI,S,V] = svd(Z,'econ');
% Compute DMD (Phi are eigenvectors)
r = length(find(cumsum(diag(S))/sum(diag(S))<=0.9999)); % dominant modes
fprintf('POD trucated at r = %i\n', r);

for k = 1:min(9,r)
    figure(4)
    subplot(3,3,k)
    Grid.cplot(PSI(1:N,k))
    title(['U POD mode ',num2str(k)])
    figure(5)
    subplot(3,3,k)
    Grid.cplot(PSI(N+1:end,k))
    title(['V POD mode ',num2str(k)])
end

%% Autocorrelation
fprintf('Analyzing Autocorrelation...\n')
figure
[R] = calculate_autocorrelation(Z,450,Parameter);
%% PSD
figure
point = randi((Parameter.N+1)^2);
N = 50;
calculate_PSD(X,Parameter,point,N,Grid)

% %% perform predictions
% load('Data/49_Test_10K_2.mat')
% Y_dash = Z(:,1);
% pred_error = [];
% 
% figure(4)
% for i = 1:SimLength
%     Y_dash = [Y_dash U*(A_s*(U'*Y_dash(:,end)))];
%     pred_error = [pred_error norm(Y_dash(:,i+1)-Y(:,i))/norm(Y(:,i))];
% %     Grid.cplot(Y_dash(:,i+1)-Y(:,i))
% %     Grid.cplot(Y_dash(:,i+1))
% %     Grid.cplot(Y(:,i))
% %     pause(0.01)
% end
% 
% % analyze error
% figure(5)
% semilogy((1:SimLength)*Parameter.DT, pred_error,'-*')


% % %% HDMD
% load('Data/Train_10K_1.mat')
% q = 3;
% [HDMD] = convert_Hankel(X,Y,q);
% [A_s, U, Phi, eigs] = DMD(HDMD.X,HDMD.Y,0.9999);
% 
% load('Data/Test_10K_1.mat')
% X_temp = X(:,1:q);
% HX = X_temp(:);
% for i = 1:SimLength
%     HY = U*(A_s*(U'*HX(:,end)));
%     HX = []
%     Y_dash = [Y_dash U*(A_s*(U'*Y_dash(:,end)))];
%     pred_error = [pred_error norm(Y_dash(:,i+1)-Y(:,i))/norm(Y(:,i))];
% 
% end


% UV -test
% figure
% for i = 1:100
%     subplot(2,1,1)
%     Grid.cplot(X(1:4900,i))
%     subplot(2,1,2)
%     Grid.cplot(X(4901:end,i))
%     pause(0.01)
% end