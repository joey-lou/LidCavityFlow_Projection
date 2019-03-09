function [R] = calculate_autocorrelation(X,N,Parameter)
% calculate autocorrelation for given data matrix
% each column is a snapshot
% input: N - length of each section for SVD procedure
%        must be smaller than total length of X
X0 = X(:,1:N);
% X0 = X0-mean(X0,2);
[~,~,V0] = svd(X0,'econ');
v0 = V0(:,1);
if max(v0)~=norm(v0,inf)
        v0 = -v0;
end
R = 1;
for i=1:30
    fprintf('Lag = %i, ',i)
    Xi = X(:,1+i:min(N+i,size(X,2)));
%     Xi = Xi-mean(Xi,2);
    [~,~,Vi] = svd(Xi,'econ');
    vi = Vi(:,1);
    if max(vi)~=norm(vi,inf)
        vi = -vi;
    end
    temp = (v0-mean(v0))'*(vi-mean(vi))/((v0-mean(v0))'*(v0-mean(v0)));
    fprintf('correlation = %5.3f\n',temp);
    if temp<=0
        R(i+1) = temp;
        break
    end
    R(i+1) = abs(temp);
    % stop if R goes below 0
end
f2 = fit(Parameter.DT*(0:length(R)-1)',R','exp1');
plot(f2,Parameter.DT*(0:length(R)-1),R,'b*')
% plot(Parameter.DT*(0:length(R)-1),R,'b-*')
hold on
line0 = linspace(0,Parameter.DT*(length(R)),100);
plot(line0,zeros(100,1),'r.')
title(['Autocorrelation N=',num2str(N)])
xlabel '\tau (s)'
ylabel 'autocorrelation'
end