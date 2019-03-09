function W = convert_UV2W(UV,DX,DY,m)
% calculate omega from u,v
u = UV(1:m,:);
v = UV(m+1:2*m,:);
W = DX*v-DY*u;
end