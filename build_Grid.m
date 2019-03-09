function [Grid] = build_Grid(N)
% Matrix constructions
Grid.m = (1+N)^2;                % grid number
[D,x] = cheb(N); 
Grid.D = 2*D;
Grid.x = (1+x)/2; Grid.y=Grid.x;           
Grid.I = eye(N+1);
[xx,yy] = meshgrid(Grid.x,Grid.y); xx=xx(:); yy=yy(:);
% find points of interest
Grid.bd_pts = find(xx==0 | xx==1 | yy==0 | yy==1); 
Grid.t_pts  = find(yy==1);   Grid.l_pts = find(xx==0);
Grid.r_pts  = find(xx==1);   Grid.b_pts = find(yy==0);
Grid.i_pts  = find(xx~=0 & xx~=1 & yy~=0 & yy~=1 );
Grid.br_pts = find(xx>=0.9 & xx~=1 & yy<=0.1 & yy~=0 ); % bottom right corner
Grid.bl_pts = find(xx<=0.1 & xx~=0 & yy<=0.1 & yy~=0 ); % bottom left corner
Grid.ul_pts = find(xx<=0.1 & xx~=0 & yy>=0.9 & yy~=1 ); % upper  right corner

Grid.cplot = @(X) contourf(reshape(xx,N+1,N+1),reshape(yy,N+1,N+1),...
    reshape(X,N+1,N+1),50,'LineStyle','None');
Grid.mplot = @(X) mesh(reshape(xx,N+1,N+1),reshape(yy,N+1,N+1),...
    reshape(X,N+1,N+1));
Grid.xx = xx; Grid.yy = yy;
Grid.II = eye(Grid.m*2);

% mid_val=Grid.x(floor(N/2)+1);
% Grid.xmid_pts = find(xx==mid_val);
% Grid.ymid_pts = find(yy==mid_val);
% Grid.D4 = (diag(1-Grid.x.^2)*Grid.D^4 - 8*diag(Grid.x)*Grid.D^3 - 12*Grid.D^2); 


end