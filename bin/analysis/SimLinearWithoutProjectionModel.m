%% Simulate the linear model obtained without projection

% Compute steady state solution with WFSim (Projection=0&&Linearversion=1)
% and save the workspace. Then run this script to simulate linear model.

clear; clc; close all;

load('..\..\Data_WFSim\LinearWithoutProjection')

Animate   = 5;
Movie     = 0;

n1      = (Wp.mesh.Nx-3)+(Wp.mesh.Nx-2);    % # cells in x-direction for one y position
n2      = size(Wp.mesh.yline{1},2);         % # cells in y-direction where the turbine is
L       = 250;
h       = Wp.sim.h;
time    = (0:h:L);
NN      = length(time);

x       = zeros(size(sys.Al,2),NN);         % [u;v] 

% Control inputs
%deltabeta   = 0.2*ones(1,NN);
deltayaw    = 0*ones(1,NN);
for kk=1:NN;if sin(2*pi*5/L*time(kk))>=0;deltabeta(kk) = .1;else deltabeta(kk) = -.1;end;end;

w           = [deltabeta;-deltabeta; ...
    deltayaw;-deltayaw]; % input signal

% Flow fields
[du,dv,ul,vl]     = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny,NN));

% ss(sys.Al,sys.Bl,I,0,sys.A,h); % Model from w to [u;v;p]


% Figure
if Animate>1
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
    if Movie==1;
        ha=gcf;
        aviobj = VideoWriter(['..\..\Data_WFSim\Compare_linearsys_with_truncated','.avi']);
        aviobj.FrameRate=5;
        open(aviobj)
    end
end

% Time loop
for kk=1:NN-1
        
    % Update plant states
    x(:,kk+1)               = sys.A\(sys.Al*x(:,kk)+sys.Bl*w(:,kk));
        
    du(3:end-1,2:end-1,kk)  = reshape(x(1:Wp.Nu,kk),Wp.mesh.Ny-2,Wp.mesh.Nx-3)';
    dv(2:end-1,3:end-1,kk)  = reshape(x(Wp.Nu+1:Wp.Nu+Wp.Nv,kk),Wp.mesh.Ny-3,Wp.mesh.Nx-2)';
        
    ul(:,:,kk)  = sol.u + du(:,:,kk);
    vl(:,:,kk)  = sol.v + dv(:,:,kk);
   
    if Animate>0
        if ~rem(kk,Animate)
            
            yaw_angles = .5*Wp.turbine.Drotor*exp(1i*w(Wp.turbine.N+1:end,kk)*pi/180);  % Yaw angles

            subplot(2,2,[1 3]);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',ul(:,:,kk),'Linecolor','none');  colormap(hot);
            caxis([min(min(ul(:,:,kk))) max(max(ul(:,:,kk)))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            text(0,Wp.mesh.ldxx2(end,end)+80,['Time ', num2str(time(kk),'%.1f'), 's']);
            xlabel('y [m]')
            ylabel('x [m]');
            title('u_l [m/s]');
            hold off;
                       
            % v-velocity component
            subplot(2,2,[2 4]);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',vl(:,:,kk),'Linecolor','none');  colormap(hot);
            caxis([min(min(vl(:,:,kk))) max(max(vl(:,:,kk)))]);  hold all; colorbar;
            axis equal; axis tight
            for ll=1:Wp.turbine.N
               Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
               Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
               plot(Qy,Qx,'k','linewidth',1)
            end
            xlabel('y [m]')
            ylabel('x [m]');
            title('v_l [m/s]')
            hold off;
            
            drawnow;
            
            if Movie==1; F = getframe(hfig); writeVideo(aviobj,F); end
        end
    end
        
end

if Movie==1;close(aviobj);end

%save('..\..\Data_WFSim\linearsol','ul','vl','z')