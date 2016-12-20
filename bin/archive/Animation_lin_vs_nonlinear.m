Dr         = Wp.turbine.Drotor;
yaw_angles = .5*Dr*exp(1i*Phi(:,k)*pi/180);  % Yaw angles
%% Plot u velocity flow component
figure(1);
subplot(2,4,1); 
contourf(Wp.ldyy(1,:),Wp.ldxx2(:,1)',u,(0:0.1:u_Inf*1.2),'Linecolor','none');  colormap(jet); caxis([min(min(u))-2 u_Inf*1.2]);  hold all; colorbar;
axis equal; axis tight;

% Plot center line wake
if (Dynamic==0 && Wp.N==1 && Verification==1)
    for kk=1:Wp.N
        plot(Wp.Cry(kk)-xw(:,kk),Wp.ldxx(Wp.xline(kk):end,1),'k--');
        plot(Wp.Cry(kk)-xw(:,kk)+Dwf(:,kk)/2,Wp.ldxx(Wp.xline(kk):end,1),'k--');
        plot(Wp.Cry(kk)-xw(:,kk)-Dwf(:,kk)/2,Wp.ldxx(Wp.xline(kk):end,1),'k--');
    end
end

% Plot the turbines in the field
for kk=1:Wp.N
    Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
    Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    plot(Qy,Qx,'k','linewidth',1)
end
text(0,Wp.ldxx2(end,end)+230,['Time ', num2str(time(k)), 's']);
ylabel('x [m]');
title('u [m/s]')
hold off;

Dr         = Wp.turbine.Drotor;
yaw_angles = .5*Dr*exp(1i*Phi(:,k)*pi/180);  % Yaw angles

%% Plot u velocity flow component
figure(1);
subplot(2,4,2);
contourf(Wp.ldyy(1,:),Wp.ldxx2(:,1)',ul,(0:0.1:u_Inf*1.2),'Linecolor','none');  colormap(jet); caxis([min(min(u))-2 u_Inf*1.2]);  hold all; colorbar;
axis equal; axis tight;

% Plot center line wake
if (Dynamic==0 && Wp.N==1 && Verification==1)
    for kk=1:Wp.N
        plot(Wp.Cry(kk)-xw(:,kk),Wp.ldxx(Wp.xline(kk):end,1),'k--');
        plot(Wp.Cry(kk)-xw(:,kk)+Dwf(:,kk)/2,Wp.ldxx(Wp.xline(kk):end,1),'k--');
        plot(Wp.Cry(kk)-xw(:,kk)-Dwf(:,kk)/2,Wp.ldxx(Wp.xline(kk):end,1),'k--');
    end
end

% Plot the turbines in the field
for kk=1:Wp.N
    Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
    Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    plot(Qy,Qx,'k','linewidth',1)
end
title('u_l [m/s]')
hold off;drawnow

%% Plot v velocity flow component
figure(1);
subplot(2,4,3);
contourf(Wp.ldyy(1,:),Wp.ldxx2(:,1)',min(v,u_Inf*1.2),(min(min(v)):0.1:max(max(v))),'Linecolor','none');  colormap(jet);   hold all
colorbar;
for kk=1:Wp.N
    Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
    Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    plot(Qy,Qx,'k','linewidth',1)
end
axis equal; axis tight
title('v [m/s]')
hold off;
drawnow;

%% Plot v velocity flow component
figure(1);
subplot(2,4,4);
contourf(Wp.ldyy(1,:),Wp.ldxx2(:,1)',min(vl,u_Inf*1.2),(min(min(vl)):0.1:max(max(vl))),'Linecolor','none');  colormap(jet);   hold all
colorbar;
for kk=1:Wp.N
    Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
    Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    plot(Qy,Qx,'k','linewidth',1)
end
axis equal; axis tight
title('v_l [m/s]')
hold off;
drawnow;

figure(1);
subplot(2,4,5);
contourf(Wp.ldyy(1,:),Wp.ldxx2(:,1)',min(du,u_Inf*1.2),(min(min(du)):0.1:max(max(v))),'Linecolor','none');  colormap(jet);   hold all
colorbar;
for kk=1:Wp.N
    Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
    Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    plot(Qy,Qx,'k','linewidth',1)
end
axis equal; axis tight
title('du [m/s]')
xlabel('y [m]')
ylabel('x [m]')
hold off;
drawnow;

figure(1);
subplot(2,4,6);
contourf(Wp.ldyy(1,:),Wp.ldxx2(:,1)',min(dv,u_Inf*1.2),(min(min(dv)):0.1:max(max(v))),'Linecolor','none');  colormap(jet);   hold all
colorbar;
for kk=1:Wp.N
    Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
    Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    plot(Qy,Qx,'k','linewidth',1)
end
axis equal; axis tight
title('dv [m/s]')
xlabel('y [m]')
hold off;
drawnow;


figure(1);
subplot(2,4,7);
contourf(Wp.ldyy(1,:),Wp.ldxx2(:,1)',min(u-ul,u_Inf*1.2),(min(min(v)):0.1:max(max(v))),'Linecolor','none');  colormap(jet);   hold all
colorbar;
for kk=1:Wp.N
    Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
    Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    plot(Qy,Qx,'k','linewidth',1)
end
axis equal; axis tight
xlabel('y [m]')
title('e_u [m/s]')
hold off;
drawnow;

figure(1);
subplot(2,4,8);
contourf(Wp.ldyy(1,:),Wp.ldxx2(:,1)',min(v-vl,u_Inf*1.2),(min(min(v)):0.1:max(max(v))),'Linecolor','none');  colormap(jet);   hold all
colorbar;
for kk=1:Wp.N
    Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
    Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    plot(Qy,Qx,'k','linewidth',1)
end
axis equal; axis tight
xlabel('y [m]')
title('e_v [m/s]')
hold off;
drawnow;


