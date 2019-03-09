function calculate_PSD(INPUT,Parameter,point,N,Grid)
% calculates power spectral density
% input: INPUT - [U;V] vector or [PSI] vector
%      : Parameter - .N for dimension
%                  - .T/.DT for sample number
%      : N - data interval for fft

% select only U for PSD calculation
X = INPUT(1:(Parameter.N+1)^2,:);
tau = 1;
PC1n = [];
TSn = [];

% find coordinate of the randomly selected point
dummy = zeros((Parameter.N+1)^2,1);
dummy(point) = 1;
dummy = reshape(dummy,Parameter.N+1,Parameter.N+1);
[x, y] = find(dummy,1);
x = Grid.x(x); y = Grid.x(y);
    
for i = 1:Parameter.runs
    fprintf('Evaluating run %i/%i...\n',i,Parameter.runs)
    Xa = X(:,((i-1)*Parameter.T/Parameter.DT+1):(i*Parameter.T/Parameter.DT));
    [V,~] = eig(Xa*Xa');
    PC1 = (V(:,end))'*Xa; % score?
    PC1 = PC1';           
    TS = squeeze(Xa(point,:))';
    % store PC1 and TS
    for m=1:floor(length(PC1)/N)
        PC1n = [ PC1n fft(PC1(N*(m-1)+1:N*m))];
        TSn  = [ TSn fft(TS(N*(m-1)+1:N*m))];
    end
end

% taking average over all runs
PSD_V = mean(abs(PC1n(1:N/2+1,:)), 2)/mean(abs(PC1n(1,:)), 2);
PSD_T = mean(abs(TSn(1:N/2+1,:)), 2)/mean(abs(TSn(1,:)), 2);
w_V   = ((0:N/2)*2*pi/(N*Parameter.DT))/(2*pi/tau);
w_T   = ((0:N/2)*2*pi/(N*Parameter.DT))/(2*pi/tau);

for i = 1:length(PSD_T)
    if(PSD_T(i) > 1)
        PSD_T(i) = 0.5*(PSD_T(i+1) + PSD_T(i-1));
    end
end

semilogy(w_T, PSD_T, 'linewidth', 2, 'color', 'blue')
hold on
semilogy(w_V, PSD_V, 'linewidth', 2, 'color', 'red')
hold on
% xlim([0 5]);
% ylim([5e-3 1])
xlabel('$\omega/\tilde{\omega}$', 'interpreter', 'latex');
ylabel('$\mathrm{PSD}$', 'interpreter', 'latex');
legend(['STATE at (',num2str(x),',',num2str(y),')'] , 'PC1');
set(gca, 'Fontsize', 14);

end