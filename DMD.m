function [A_s, U, Phi, eigs] = DMD(X,Y,tol)
% tol specifies the percentage of singular values
[U,S,V] = svd(X,'econ');
% Compute DMD (Phi are eigenvectors)
r = length(find(cumsum(diag(S))/sum(diag(S))<=tol)); % dominant modes
fprintf('DMD trucated at r = %i\n', r);
U = U(:,1:r);
S = S(1:r,1:r);
V = V(:,1:r);
A_s = U'*Y*V*inv(S);
% A_s = U'*Y*V*diag((1./diag(S)));
[W,eigs] = eig(A_s);
Phi.exact = Y*V*inv(S)*W; % Exact Modes
% Phi.projected = U*W; % Projected Modes
end