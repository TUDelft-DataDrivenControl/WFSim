xline  = Wp.mesh.xline;
yline  = Wp.mesh.yline;
Drotor = Wp.turbine.Drotor;
N      = Wp.turbine.N;

% Determine mixing length distribution in the field

mixing_length = ConstructLmu(Wp.mesh.ldxx2,Wp.mesh.ldyy,tan(Wp.site.v_Inf/Wp.site.u_Inf),...
    Wp.turbine.Crx,Wp.turbine.Cry,Wp.turbine.Drotor,...
    Wp.site.d_lower,Wp.site.d_upper,Wp.site.lm_slope);

% figure; surf(Wp.mesh.ldyy2,Wp.mesh.ldxx,mixing_length); axis equal; xlabel('y'); ylabel('x');

% include turbulence model in equations
% For u-momentum equation
ax.Tnx              = zeros(Nx,Ny);
ax.Tsx              = zeros(Nx,Ny);

ax.Tnx(2:Nx,1:Ny-1) = Rho*(mixing_length(2:Nx,1:Ny-1).^2).*(dxx(2:Nx,1:Ny-1)./(dyy(2:Nx,2:Ny).^2)).*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1));
ax.Tsx(1:Nx-1,2:Ny) = Rho*(mixing_length(1:Nx-1,2:Ny).^2).*(dxx(2:Nx,2:Ny)./(dyy(2:Nx,2:Ny).^2)).*abs(u(2:Nx,1:Ny-1)-u(2:Nx,2:Ny));

ax.aN             = ax.aN + ax.Tnx;
ax.aS             = ax.aS + ax.Tsx;
ax.aP             = ax.aP + ax.Tnx + ax.Tsx;

% For v-momentum equation
ay.Tey            = zeros(Nx,Ny);
ay.Twy            = zeros(Nx,Ny);

ay.Tey(1:Nx-1,1:Ny) = Rho*(mixing_length(1:Nx-1,1:Ny).^2).*(dyy(1:Nx-1,1:Ny)./(dxx(1:Nx-1,1:Ny).^2)).*abs(v(2:Nx,1:Ny)-v(1:Nx-1,1:Ny));
ay.Twy(2:Nx,1:Ny)   = Rho*(mixing_length(2:Nx,1:Ny).^2).*(dyy(2:Nx,1:Ny)./(dxx(2:Nx,1:Ny).^2)).*abs(v(1:Nx-1,1:Ny)-v(2:Nx,1:Ny));

ay.aE             = ay.aE + ay.Tey;
ay.aW             = ay.aW + ay.Twy;
ay.aP             = ay.aP + ay.Tey + ay.Twy;


if Linearversion
    
    % For u-momentum equation
    %dTsxd1           = zeros(Nx,Ny);
    %dTsxd2           = zeros(Nx,Ny);
    %dTnxd1           = zeros(Nx,Ny);
    %dTnxd2           = zeros(Nx,Ny);
    
    % dTsx/du_(i,J)
    %dTsxd1(1:Nx,2:Ny)   = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx2(1:Nx,2:Ny)./(dyy(1:Nx,2:Ny).^2)).*sign((u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny)));
    % dTsx/du_(i,J-1)
    %dTsxd2(1:Nx,2:Ny)   = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx2(1:Nx,2:Ny)./(dyy(1:Nx,2:Ny).^2)).*-sign((u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny)));
    % dTnx/du_(i,J)
    %dTnxd1(1:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx2(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2)).*-sign((u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1)));
    % dTnx/du_(i,J+1)
    %dTnxd2(1:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx2(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2)).*sign((u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1)));
    
    dax.S(1:Nx,2:Ny)   = dax.S(1:Nx,2:Ny)   + ax.Tsx(1:Nx,2:Ny) ;
    dax.N(1:Nx,1:Ny-1) = dax.N(1:Nx,1:Ny-1) + ax.Tnx(1:Nx,1:Ny-1);
    dax.P(1:Nx,1:Ny-1) = dax.P(1:Nx,1:Ny-1) + ax.Tnx(1:Nx,1:Ny-1) + ax.Tsx(1:Nx,1:Ny-1) ;
    
    %dax.S(1:Nx,2:Ny)   = dax.S(1:Nx,2:Ny)+dTsxd2(1:Nx,2:Ny).*u(1:Nx,2:Ny)-dTsxd2(1:Nx,2:Ny).*u(1:Nx,1:Ny-1);
    %dax.N(1:Nx,1:Ny-1) = dax.N(1:Nx,1:Ny-1)-dTnxd2(1:Nx,1:Ny-1).*u(1:Nx,1:Ny-1)+dTnxd2(1:Nx,1:Ny-1).*u(1:Nx,2:Ny);
    %dax.P(1:Nx,1:Ny-1) = dax.P(1:Nx,1:Ny-1)+dTnxd1(1:Nx,1:Ny-1).*u(1:Nx,1:Ny-1)+dTsxd1(1:Nx,1:Ny-1).*u(1:Nx,1:Ny-1) ...
    %-dTnxd1(1:Nx,1:Ny-1).*u(1:Nx,2:Ny)-dTsxd1(1:Nx,1:Ny-1).*u(1:Nx,2:Ny);
    
    % For v-momentum equation
    day.E(1:Nx,2:Ny)   = day.E(1:Nx,2:Ny)   + ay.Tey(1:Nx,2:Ny) ;
    day.W(1:Nx,1:Ny-1) = day.W(1:Nx,1:Ny-1) + ay.Twy(1:Nx,1:Ny-1);
    day.P(1:Nx,1:Ny-1) = day.P(1:Nx,1:Ny-1) + ay.Tey(1:Nx,1:Ny-1) + ay.Twy(1:Nx,1:Ny-1) ;
end;
