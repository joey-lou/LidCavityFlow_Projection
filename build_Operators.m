function [Operators] = build_Operators(Grid,Parameter)
% differentiation matrices for UV
Operators.DX  = kron(Grid.D,Grid.I); Operators.DY = kron(Grid.I,Grid.D); 
D2 = Grid.D^2;
L = kron(D2,Grid.I)+kron(Grid.I,D2);
Operators.DEL2 = ...
    [L, zeros(Grid.m)
    zeros(Grid.m), L];
% differentiation matrix for Pressure (L)
L([Grid.t_pts; Grid.b_pts],:) = Operators.DY([Grid.t_pts; Grid.b_pts],:);
L([Grid.l_pts; Grid.r_pts],:) = Operators.DX([Grid.l_pts; Grid.r_pts],:);
% L_new = [L; [1  zeros(1,size(L,1)-1)]];
Operators.pL = pinv(L);
% operator for converting omega to psi
L(Grid.bd_pts,:) = 0;
L(Grid.bd_pts,Grid.bd_pts) = eye(4*Parameter.N);
Operators.oL = inv(L);
Operators.L = L;
end