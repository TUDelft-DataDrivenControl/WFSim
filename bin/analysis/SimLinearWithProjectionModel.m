%% Simulate the linear model

% Compute steady state solution with WFSim (Projection=1&&Linearversion=1)
% and save the workspace using: clear hfig;save('Data_WFSim\nonlinearsol').
% Then run this script to simulate linear model.

clear; clc; close all;

load('..\..\Data_WFSim\nonlinearsol')

Animate   = 5;
Movie     = 0;

nu      = Wp.Nu;
nv      = Wp.Nv;
np      = Wp.Np;
n1      = (Wp.mesh.Nx-3)+(Wp.mesh.Nx-2);    % # cells in x-direction for one y position
n2      = size(Wp.mesh.yline{1},2);         % # cells in y-direction where the turbine is
L       = 250;
h       = Wp.sim.h;
time    = (0:h:L);
NN      = length(time);

x       = zeros(size(sys.Qsp,2),NN);            % alpha state
x1      = zeros(size(sys.Qsp,2),NN);            % alpha state truncated sys
x2      = zeros(np+nv+nu,NN);                   % state for power equation
y       = zeros(nu+nv,NN);                      % [u;v]
y1      = zeros(nu+nv,NN);                      % [u;v] truncated sys
z       = zeros(1,NN);                          % Performance channel with power signals
nz      = size(z,1);

% Control inputs
deltabeta   = 0.2*ones(1,NN);
deltayaw    = 0*ones(1,NN);
%for kk=1:NN;if sin(2*pi*5/L*time(kk))>=0;deltabeta(kk) = .1;else deltabeta(kk) = -.1;end;end;

w           = [deltabeta;-deltabeta; ...
    deltayaw;-deltayaw]; % input signal

% Flow fields
[du,dv,ul,vl]     = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny,NN));
[du1,dv1,ul1,vl1] = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny,NN));

% System matrices
A          = sys.Etl\sys.Atl;
B          = sys.Etl\sys.Btl;
A(A<10e-2) = 0;

Ctl        = sys.Qsp;
Cz         = zeros(1,nu+nv);
Cz((Wp.mesh.xline(2)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{2}(1)-1:...
    (Wp.mesh.xline(2)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{2}(end)-1) = 1; % Extracts sum of rotor velocities second turbine
C          = Cz*Ctl;

% ss(Atl,Btl,Ctl,0,Etl,h); % Model from w to y
% ss(A,B,Ctl,0,h); % Model from w to y1
% ss(A,B,C,0,h); % Model from w to z


% Power model (have to be checked)
%Em    = (derivatives.A);
%Am    = blkdiag(spdiags(ccx,0,length(ccx),length(ccx)),spdiags(ccy,0,length(ccy),length(ccy)),0.*speye((Wp.Nx-2)*(Wp.Ny-2)-2));
%Am    = Am+derivatives.dSmdx-derivatives.dAdx;
%Bm    = derivatives.dSmdbeta;
%Cp    = -derivatives.dJdx';
%Dp    = -derivatives.dJdbeta;

%Ap = Em\Am;
%Bp = Em\Bm;

%Ap(Ap<10e-2)=0;
%Model = ss(Am, Bm, Cm,Dm,Em,h); % Model from w to z

%[Qsp,Bsp]  = Solution_space(B1,B2,bc);
%qq         = length(ccx)+length(ccy);
%Model2     = ss(Qsp'*Am(1:qq,1:qq)*Qsp, Qsp'*Bm(1:qq,:), Cm(:,1:qq)*Qsp,Dm,Qsp'*Em(1:qq,1:qq)*Qsp,h);

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
       
    y(:,kk)    = Ctl*x(:,kk);   % Output full system
    y1(:,kk)   = Ctl*x1(:,kk);  % Output truncated system
    z(:,kk)    = Cz*y1(:,kk);   % Sum rotor velocities turbine 2
    
    % Update plant states
    x(:,kk+1)  = sys.Etl\(sys.Atl*x(:,kk) + sys.Btl*w(:,kk));
    x1(:,kk+1) = A*x1(:,kk) + B*w(:,kk);
    
    
    du(3:end-1,2:end-1,kk)  = reshape(y(1:nu,kk),Wp.mesh.Ny-2,Wp.mesh.Nx-3)';
    dv(2:end-1,3:end-1,kk)  = reshape(y(nu+1:nu+nv,kk),Wp.mesh.Ny-3,Wp.mesh.Nx-2)';
    du1(3:end-1,2:end-1,kk) = reshape(y1(1:nu,kk),Wp.mesh.Ny-2,Wp.mesh.Nx-3)';
    dv1(2:end-1,3:end-1,kk) = reshape(y1(nu+1:nu+nv,kk),Wp.mesh.Ny-3,Wp.mesh.Nx-2)';
    
    ul(:,:,kk)  = sol.u + du(:,:,kk);
    vl(:,:,kk)  = sol.v + dv(:,:,kk);
    ul1(:,:,kk) = sol.u + du1(:,:,kk);
    vl1(:,:,kk) = sol.v + dv1(:,:,kk);
    
    %z(:,kk) - sum(du1(Wp.mesh.xline(2),Wp.mesh.yline{2},kk)) % Has to be zero
    
    if Animate>0
        if ~rem(kk,Animate)
            
            yaw_angles = .5*Wp.turbine.Drotor*exp(1i*w(Wp.turbine.N+1:end,kk)*pi/180);  % Yaw angles

            subplot(2,3,[1 4]);
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
            
            subplot(2,3,[2 5]);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',ul1(:,:,kk),'Linecolor','none');  colormap(hot);
            caxis([min(min(ul1(:,:,kk))) max(max(ul1(:,:,kk)))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            xlabel('y [m]')
            ylabel('x [m]');
            title('u_l [m/s] truncated sys');
            hold off;
            
            subplot(2,3,[3 6]);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',ul(:,:,kk)-ul1(:,:,kk),'Linecolor','none');  colormap(hot);
            caxis([min(min(ul(:,:,kk)-ul1(:,:,kk))) max(max(ul(:,:,kk)-ul1(:,:,kk)))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            xlabel('y [m]')
            ylabel('x [m]');
            title('error [m/s]');
            hold off;
            
            % The following is for the v-velocity component
            %contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',vl(:,:,kk),'Linecolor','none');  colormap(jet);
            %caxis([min(min(vl(:,:,kk))) max(max(vl(:,:,kk)))]);  hold all; colorbar;
            %axis equal; axis tight
            %for ll=1:Wp.turbine.N
            %    Qy     = (Wp.mesh.Cry(ll)-real(yaw_angles(ll))):1:(Wp.mesh.Cry(ll)+real(yaw_angles(ll)));
            %    Qx     = linspace(Wp.mesh.Crx(ll)-imag(yaw_angles(ll)),Wp.mesh.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
            %    plot(Qy,Qx,'k','linewidth',1)
            %end
            %xlabel('y [m]')
            %ylabel('x [m]');
            %title('v_l [m/s]')
            %hold off;
            
            drawnow;
            
            if Movie==1; F = getframe(hfig); writeVideo(aviobj,F); end
        end
    end
    
    RMSE(kk)               = rms(x(:,kk)-x1(:,kk));
    [maxe(kk),maxeloc(kk)] = max(abs(x(:,kk)-x1(:,kk)));
    
end

if Movie==1;close(aviobj);end

save('..\..\Data_WFSim\linearsol','ul','vl','z')

figure(2);clf
subplot(2,1,1)
plot(time(1:end-1),RMSE);hold on;
plot(time(1:end-1),maxe,'r');grid;
ylabel('RMSE and max');
title('\color{blue}RMSE, \color{red}max')
subplot(2,1,2)
plot(time,deltabeta);hold on;
plot(time,-deltabeta,'r');grid;
ylabel('\delta \beta');xlabel('Time [s]');
title('\color{blue}\delta \beta_1, \color{red}\delta \beta_2')
