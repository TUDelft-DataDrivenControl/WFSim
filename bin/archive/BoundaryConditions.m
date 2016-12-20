function [ax,ay,bx,by] = BoundaryConditions(Wp,ax,ay,u,v)

Nx = Wp.Nx;
Ny = Wp.Ny;

% Zero gradient outflow
% for u-direction
%ax.aP((Nx-1),(1:Ny))  = ax.aP((Nx-1),(1:Ny)) - ax.aE((Nx),(1:Ny));
%ax.aP((1:Nx),(Ny-1))  = ax.aP((1:Nx),(Ny-1)) - ax.aN((1:Nx),(Ny));
%ax.aP((1:Nx),(2))     = ax.aP((1:Nx),2)      - ax.aS((1:Nx),1);

% for v-direction
%ay.aP((Nx-1),(1:Ny))  = ay.aP((Nx-1),(1:Ny)) - ay.aE((Nx),(1:Ny));
%ay.aP((1:Nx),(Ny-1))  = ay.aP((1:Nx),(Ny-1)) - ay.aN((1:Nx),(Ny));
%ay.aP((1:Nx),(2))     = ay.aP((1:Nx),2)      - ay.aS((1:Nx),1);

% Inflow boundary for non linear model
bx      = kron([1;zeros(Nx-4,1)],(ax.aW(2,2:end-1).*u(2,2:end-1))');
by      = [v(1,3:Ny-1)'.*ay.aW(1,3:Ny-1)';zeros((Nx-3)*(Ny-3),1)];

end