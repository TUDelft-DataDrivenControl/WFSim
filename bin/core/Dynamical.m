function [StrucDiscretization,StrucDynamical] = Dynamical(Wp,StrucDiscretization,sol,dt,Linearversion)

Rho  = Wp.site.Rho;
dxx  = Wp.mesh.dxx;
dyy  = Wp.mesh.dyy;
dxx2 = Wp.mesh.dxx2;
dyy2 = Wp.mesh.dyy2;
Nx   = Wp.mesh.Nx;
Ny   = Wp.mesh.Ny;
u    = sol.uk;
v    = sol.vk;

% Fully implicit (page 248 Versteeg) See also page 257
StrucDiscretization.ax.aP   = StrucDiscretization.ax.aP + Rho.*dxx.*dyy2./dt;       % Rho.*dxx.*dyy2./dt = a_P^0
StrucDiscretization.ay.aP   = StrucDiscretization.ay.aP + Rho.*dxx2.*dyy./dt;

StrucDynamical.ccx    = vec(Rho.*dxx(3:end-1,2:end-1)'.*dyy2(3:end-1,2:end-1)'./dt);
StrucDynamical.cx     = StrucDynamical.ccx.*vec(u(3:end-1,2:end-1)');
StrucDynamical.ccy    = vec(Rho.*dxx2(2:end-1,3:end-1)'.*dyy(2:end-1,3:end-1)'./dt);
StrucDynamical.cy     = StrucDynamical.ccy.*vec(v(2:end-1,3:end-1)');

if Linearversion
    StrucDiscretization.dax.P   = StrucDiscretization.dax.P + Rho.*dxx.*dyy2./dt;       
    StrucDiscretization.day.P   = StrucDiscretization.day.P + Rho.*dxx2.*dyy./dt;
    StrucDynamical.dcdx         = blkdiag(spdiags(StrucDynamical.ccx,0,length(StrucDynamical.ccx),length(StrucDynamical.ccx)),...
                                             spdiags(StrucDynamical.ccy,0,length(StrucDynamical.ccy),length(StrucDynamical.ccy)),...
                                             sparse((Ny-2)*(Nx-2),(Ny-2)*(Nx-2)) );
end;

