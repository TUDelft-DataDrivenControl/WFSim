%% Plot u velocity flow component
figure(1);
subplot(1,2,1);
contourf(Wp.ldyy2(1,:),Wp.ldxx2(:,1)',u,(0:0.1:u_Inf*1.2),'Linecolor','none');  colormap(jet); caxis([min(min(u))-2 u_Inf*1.2]);  hold all; colorbar;
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
    %Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
    %Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    Qy  = [Wp.Cry(kk)-Wp.turbine.Drotor/2; Wp.Cry(kk)+Wp.turbine.Drotor/2];
    Qx  = [Wp.Crx(kk); Wp.Crx(kk)];
    
    plot(Qy,Qx,'k','linewidth',1)
end
text(0,Wp.ldxx2(end,end)+80,['Time ', num2str(dt*k), 's']);
xlabel('y-direction [m]')
ylabel('x-direction [m]');
title('u-velocity')

% Plot the streamlines
if streamlines==1;
    for kk=1:Wp.N
        XY = stream2(Wp.ldyy2(1,:),Wp.ldxx2(:,1)',v,u,[(Wp.Cry(kk)-0.5*Ar);(Wp.Cry(kk)+0.5*Ar)],[Wp.Crx(kk);Wp.Crx(kk)]);
        streamline(XY)
        XY1 = stream2(Wp.ldyy2(1,:),Wp.ldxx2(:,1)',-v,-u,[(Wp.Cry(kk)-0.5*Ar);(Wp.Cry(kk)+0.5*Ar)],[Wp.Crx(kk);Wp.Crx(kk)]);
        streamline(XY1)
    end
end
hold off;

%% Plot the v velocity flow component
subplot(1,2,2);
contourf(Wp.ldyy2(1,:),Wp.ldxx2(:,1)',min(v,u_Inf*1.2),(min(min(v)):0.1:max(max(v))),'Linecolor','none');  colormap(jet);   hold all
colorbar;
for kk=1:Wp.N
%     Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
%     Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
    Qy  = [Wp.Cry(kk)-Wp.turbine.Drotor/2; Wp.Cry(kk)+Wp.turbine.Drotor/2];
    Qx  = [Wp.Crx(kk); Wp.Crx(kk)];
    plot(Qy,Qx,'k','linewidth',1)
end
axis equal; axis tight
text(0,Wp.ldxx2(end,end)+80,['Time ', num2str(dt*k), 's']);
xlabel('y-direction [m]')
ylabel('x-direction [m]');
title('v-velocity')
hold off;
drawnow;

% figure(2);
% for kk=1:Wp.N
%     Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
%     Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
%     plot(Qy,Qx,'k','linewidth',1);hold on;
% end
% quiver(Wp.ldyy2,Wp.ldxx2,v,u);grid;
% xlabel('y-direction [m]');ylabel('x-direction [m]');
% axis tight
% drawnow
