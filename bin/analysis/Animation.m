Dr     = Wp.turbine.Drotor;
Nx     = Wp.mesh.Nx;
Ny     = Wp.mesh.Ny;
dxx    = Wp.mesh.dxx;
dyy    = Wp.mesh.dyy;
dxx2   = Wp.mesh.dxx2;
dyy2   = Wp.mesh.dyy2;
ldxx   = Wp.mesh.ldxx;
ldyy   = Wp.mesh.ldyy;
ldxx2  = Wp.mesh.ldxx2;
ldyy2  = Wp.mesh.ldyy2;
xline  = Wp.mesh.xline;
yline  = Wp.mesh.yline;

Rho    = Wp.site.Rho;
mu     = Wp.site.mu;
lmu    = Wp.site.lmu;
Turb   = Wp.site.turbul;
u_Inf  = Wp.site.u_Inf;
v_Inf  = Wp.site.v_Inf;

Drotor = Wp.turbine.Drotor;
N      = Wp.turbine.N;
Cry    = Wp.turbine.Cry;
Crx    = Wp.turbine.Crx;

time   = Wp.sim.time;

yaw_angles = .5*Dr*exp(1i*input{k}.phi*pi/180);  % Yaw angles

%% Plot u velocity flow component
figure(hfig);
subplot(2,2,[1 3]);
contourf(ldyy(1,:),ldxx2(:,1)',sol.u,(0:0.1:u_Inf*1.2),'Linecolor','none');  colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.04]);  hold all; colorbar;
axis equal; axis tight;

if N==1
    % Compute center line wake
    kd = .15;
    ad = -4.5;
    bd = -0.01;
    
    X        = ldxx(xline(1),1);
    x        = ldxx(xline(1):end,1);
    temp     = 2*kd*(x-X)/Dr +1;
    xinit    = .5*cos(input{k}.phi(1)*pi/180)^2*sin(input{k}.phi(1)*pi/180)*CT(k);
    xw(:,k)  = (xinit*15*temp.^4 + xinit^2) ./(30*kd/Dr*temp.^5)-...
        xinit*Dr*(15+xinit^2)/(30*kd);%+ad+bd*(x-X);
    
    % Compute theoretical velocity deficit for one turbine
    ke = .065;
    MU = 5.5;   % For mixing region 5.5, near wake region .5 and far wale region 1
    aU = 5;
    bU = 1.66;
    kf = .075;
    
    X     = ldxx(xline(1),1);
    x     = ldxx(xline(1):end,1);
    mU    = MU/(cos(aU+bU*input{k}.phi(1)*pi/180));
    
    Vwake   = Ueffect(1,k)*(1-((1-(1-2*a(1,k)))./(1+2*kf*(x-X)./Dr).^2)); % Jensen
    
    % Compute theoretical wake diameter
    Dwm(:,k) = max(Dr+2*ke*1*(x-X),0);     % mixing region
    Dwf(:,k) = max(Dr+2*ke*.22*(x-X),0);   % far wake
    Dwn(:,k) = max(Dr+2*ke*-.5*(x-X),0);   % near wake
    
    % Plot center line wake
    plot(Cry(1)-xw(:,1),ldxx(xline(1):end,1),'k--');
    plot(Cry(1)-xw(:,1)+Dwf(:,1)/2,ldxx(xline(1):end,1),'k--');
    plot(Cry(1)-xw(:,1)-Dwf(:,1)/2,ldxx(xline(1):end,1),'k--');
    vline(ldyy(1,yline{1}(1)),'b:')
    vline(ldyy(1,yline{1}(end)),'b:')
end

% Plot the turbines in the field
for kk=1:N
    Qy     = (Cry(kk)-real(yaw_angles(kk))):1:(Cry(kk)+real(yaw_angles(kk)));
    Qx     = linspace(Crx(kk)-imag(yaw_angles(kk)),Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    plot(Qy,Qx,'k','linewidth',1)
end
text(0,ldxx2(end,end)+80,['Time ', num2str(time(k),'%.1f'), 's']);
xlabel('y [m]')
ylabel('x [m]');
title('u [m/s]');
hold off;

%% Plot the v velocity flow component
subplot(2,2,2);
contourf(ldyy(1,:),ldxx2(:,1)',min(sol.v,u_Inf*1.2),min(sol.v(:)):0.1:max(sol.v(:)),'Linecolor','none');  colormap(hot);   hold all
% contourf(ldyy(1,:),ldxx2(:,1)',min(v,u_Inf*1.2),'Linecolor','none');  colormap(hot);   hold all
colorbar;
for kk=1:N
    Qy     = (Cry(kk)-real(yaw_angles(kk))):1:(Cry(kk)+real(yaw_angles(kk)));
    Qx     = linspace(Crx(kk)-imag(yaw_angles(kk)),Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    plot(Qy,Qx,'k','linewidth',1)
end
axis equal; axis tight
xlabel('y [m]')
ylabel('x [m]');
title('v [m/s]')
hold off;
drawnow;

%% Wake mean centreline first turbine
D_ind    = yline{1};
cent_ind = ceil(0.5*(D_ind(1)+D_ind(end)));
up(:,k)  = mean(sol.u(:,D_ind),2);

figure(hfig);
subplot(2,2,4)
plot(ldxx2(:,1)',up(:,k),'k','Linewidth',2);
xlabel('x [m]');ylabel('U_c [m/s]');grid;
ylim([0 u_Inf+1]);xlim([ldxx2(1,1) ldxx2(end,1)]);
for kk=1:N
    vline(Crx(kk),'r--')
    vline(Crx(kk)+6*Dr,'b:','6D')
end
if N==1
    hline(u_Inf,'k--','u_b');
    hline(u_Inf*(1-1/3),'k--','u_b(1-a)');
    hline(u_Inf*(1-2*1/3),'k--','u_b(1-2a)');
end
hold off;
drawnow;
