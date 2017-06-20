function [u,v,p] = Updateboundaries(Nx,Ny,u,v,p)
% u = [u_(1,1) u_(1,2) u_(1,3) .... u_(1,Ny) (dummy row)
%      u_(2,1) u_(2,2) u_(2,3) .... u_(2,Ny) (bc put on this row)
%       ...
%      u_(Nx,1) u_(Nx,2) u_(Nx,3) .... u_(Nx,Ny)]
% v = [v_(1,1) v_(1,2) v_(1,3) .... v_(1,Ny)
%      v_(2,1) v_(2,2) v_(2,3) .... v_(2,Ny)
%       ...
%      v_(Nx,1) v_(Nx,2) v_(Nx,3) .... v_(Nx,Ny)] (first column is dummy column, bc is put on second column)

% Three zero gradient boundary conditions
u(:,1)      =  u(:,2);          % u_{i,1}  = u_{i+1,2}   for i = 1,..Nx
u(:,Ny)     =  u(:,Ny-1);       % u_{i,Ny} = u_{i,Ny-1}  for i = 1,..Nx
u(Nx,:)     =  u(Nx-1,:);       % u_{Nx,J} = u_{Nx-1,J}  for J = 1,..Ny (hence u_Inf comes via first row in field)

v(:,2)      =  v(:,3);
v(:,1)      =  v(:,2);          % v_{I,1}  = v_{I,2}    for I = 1,..Nx
v(:,Ny)     =  v(:,Ny-1);       % v_{I,Ny} = v_{I,Ny-1} for I = 1,..Nx
v(Nx,:)     =  v(Nx-1,:);       % v_{Nx,j} = v_{Nx,j}   for j = 1,..Ny (hence v_Inf comes via row in field)

p(:,2)      =  p(:,3);          % Trick to make pressure field nice
p(:,Ny-1)   =  p(:,Ny-2);       % Trick to make pressure field nice
p(Nx-1,:)   =  p(Nx-2,:);       % Trick to make pressure field nice
p(:,1)      =  p(:,2);          % p_{i,1}  = p_{i+1,2}   for i = 1,..Nx
p(:,Ny)     =  p(:,Ny-1);       % p_{i,Ny} = p_{i,Ny-1}  for i = 1,..Nx
p(Nx,:)     =  p(Nx-1,:);       % p_{Nx,J} = p_{Nx-1,J}  for J = 1,..Ny

% p = p + (rand(size(p))-.5)/.5*1e-2;
% v = v + (rand(size(u))-.5)/.5*5e-1;

end
