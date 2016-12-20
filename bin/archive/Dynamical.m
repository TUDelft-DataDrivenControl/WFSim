function [ax,ay,cx,cy,ccx,ccy] = Dynamical(Wp,ax,ay,uk,vk,dt,Rho)
% Fully implicit (page 248 Versteeg) See also page 257
ax.aP   = ax.aP + Rho.*Wp.dxx.*Wp.dyy2./dt;       % Rho.*dxx.*dyy2./dt = a_P^0
ay.aP   = ay.aP + Rho.*Wp.dxx2.*Wp.dyy./dt;

% Use uk, vk since they are from the previous time step
cx      = vec(Rho.*Wp.dxx(3:end-1,2:end-1)'.*Wp.dyy2(3:end-1,2:end-1)'./dt).*vec(uk(3:end-1,2:end-1)');
cy      = vec(Rho.*Wp.dxx2(2:end-1,3:end-1)'.*Wp.dyy(2:end-1,3:end-1)'./dt).*vec(vk(2:end-1,3:end-1)');
ccx     = vec(Rho.*Wp.dxx(3:end-1,2:end-1)'.*Wp.dyy2(3:end-1,2:end-1)'./dt);
ccy     = vec(Rho.*Wp.dxx2(2:end-1,3:end-1)'.*Wp.dyy(2:end-1,3:end-1)'./dt);


