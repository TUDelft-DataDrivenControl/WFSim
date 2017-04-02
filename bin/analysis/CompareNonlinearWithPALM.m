clear; clc; close all;

[WFSimFolder, ~, ~] = fileparts(which([mfilename '.m']));   % Get WFSim directory
addpath(genpath('..\..\Data_PALM\'));                       % Add PALM data

% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 1;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name             = 'wfcontrol_2turb';
Wp.Turbulencemodel  = 'WFSim1';

Animate       = 25;                      % Show 2D flow fields every x iterations (0: no plots)
plotMesh      = 0;                      % Show meshing and turbine locations
conv_eps      = 1e-6;                   % Convergence threshold
max_it_dyn    = 1;                      % Maximum number of iterations for k > 1

if options.startUniform==1
    max_it = 1;
else
    max_it = 50;
end

% WFSim general initialization script
[Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] ...
    = InitWFSim(Wp,options,plotMesh);

% Initialize variables and figure specific to this script
uk = Wp.site.u_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
vk = Wp.site.v_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
pk = Wp.site.p_init*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);

if Animate > 0
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
end

%% Get PALM data
filename = 'wfcontrol_2turb.nc'; nc_dump(filename);
u        = double(nc_varget(filename,'u'));
v        = double(nc_varget(filename,'v'));
w        = double(nc_varget(filename,'w'));

Turbine_settings_wfcontrol_2turb;

%% Loop
for k=1:Wp.sim.NN
    tic
    it        = 0;
    eps       = 1e19;
    epss      = 1e20;
    
    % Write flow field solutions to a 3D matrix
    uk(:,:,k) = sol.u;
    vk(:,:,k) = sol.v;
    pk(:,:,k) = sol.p;
    
    while ( eps>conv_eps && it<max_it && eps<epss );
        it   = it+1;
        epss = eps;
        
        if k>1
            max_it = max_it_dyn;
        end
        
        [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k)] = ...
            Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options); % Create system matrices
        [sol,sys] = Computesol(sys,input{k},sol,k,it,options);                   % Compute solution
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);         % Map solution to field
        Phi(:,k)  = input{k}.phi;
    end
    toc
    
    uPALM = reshape(u(k,:,:,:),77,221)';
    vPALM = reshape(v(k,:,:,:),77,221)';
    
    eu                    = vec(sol.u-uPALM); eu(isnan(eu)) = [];
    ev                    = vec(sol.v-vPALM); ev(isnan(ev)) = [];
    RMSE(k)               = rms([eu;ev]);
    [maxe(k),maxeloc(k)]  = max(abs(eu));
    
    PowerPALM(:,k)       = [Turbine1(k,8);Turbine2(k,8)];
    
    
    if Animate > 0
        if ~rem(k,Animate)
            
            yaw_angles = .5*Wp.turbine.Drotor*exp(1i*input{k}.phi*pi/180);  % Yaw angles
            
            subplot(2,3,[1 4]);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.u,'Linecolor','none');  colormap(hot);
            caxis([min(min(sol.u)) max(max(sol.u))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            text(0,Wp.mesh.ldxx2(end,end)+250,['Time ', num2str(Wp.sim.time(k),'%.1f'), 's']);
            ylabel('x [m]');
            title('WFSim u [m/s]');
            hold off;
            
            subplot(2,3,[2 5]);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',uPALM,'Linecolor','none');  colormap(hot);
            caxis([min(min(uPALM)) max(max(uPALM))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            title('PALM u [m/s]');
            hold off;
            
            subplot(2,3,[3 6]);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.u-uPALM,'Linecolor','none');  colormap(hot);
            caxis([min(min(sol.u-uPALM)) max(max(sol.u-uPALM))]);  hold all; colorbar;
            axis equal; axis tight;
            ldyyv = Wp.mesh.ldyy(:); ldxx2v = Wp.mesh.ldxx2(:);
            plot(ldyyv(maxeloc(k)),ldxx2v(maxeloc(k)),'whiteo','LineWidth',1,'MarkerSize',8,'DisplayName','Maximum error location');
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            title('error [m/s]');
            hold off;
            
%             subplot(2,3,4);
%             contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.v,'Linecolor','none');  colormap(hot);
%             caxis([min(min(sol.v)) max(max(sol.v))]);  hold all; colorbar;
%             axis equal; axis tight;
%             for ll=1:Wp.turbine.N
%                 Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
%                 Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
%                 plot(Qy,Qx,'k','linewidth',1)
%             end
%             xlabel('y [m]')
%             ylabel('x [m]');
%             title('WFSim v [m/s]');
%             hold off;
%             
%             subplot(2,3,5);
%             contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',vPALM,'Linecolor','none');  colormap(hot);
%             caxis([min(min(vPALM)) max(max(vPALM))]);  hold all; colorbar;
%             axis equal; axis tight;
%             for ll=1:Wp.turbine.N
%                 Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
%                 Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
%                 plot(Qy,Qx,'k','linewidth',1)
%             end
%             xlabel('y [m]')
%             title('PALM v [m/s]');
%             hold off;
%             
%             subplot(2,3,6);
%             contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.v-vPALM,'Linecolor','none');  colormap(hot);
%             caxis([min(min(sol.v-vPALM)) max(max(sol.v-vPALM))]);  hold all; colorbar;
%             axis equal; axis tight;
%             for ll=1:Wp.turbine.N
%                 Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
%                 Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
%                 plot(Qy,Qx,'k','linewidth',1)
%             end
%             xlabel('y [m]')
%             title('error [m/s]');
%             hold off;
            drawnow;
            
        end;
    end;
end;

%%
figure(2);clf
plot(Wp.sim.time(1:end-1),RMSE);hold on;
plot(Wp.sim.time(1:end-1),maxe,'r');grid;
ylabel('RMSE and max');
title(['{\color{blue}{RMSE}}, {\color{red}{max}} and meanRMSE = ',num2str(mean(RMSE),3)])

%% Wake centreline
D_ind    = Wp.mesh.yline{1};
indices  = [50 125 200 300];

for k=indices
    up(:,k)      = mean(uk(:,D_ind,k),2);
    uPALM        = reshape(u(k,:,:,:),77,221)';
    upPALM(:,k)  = mean(uPALM(:,D_ind),2);
    VAF(:,k)     = vaf(upPALM(:,k),up(:,k));
end


figure(3);clf;
subplot(2,2,1)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(1)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(1)),'b','Linewidth',1);grid;
ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(7)); end
title( ['VAF = ',num2str(VAF(indices(1)),3), '\% at $k$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
subplot(2,2,2)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(2)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(2)),'b','Linewidth',1);grid;
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(7)); end
title( ['VAF = ',num2str(VAF(indices(2)),3), '\% at $k$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
subplot(2,2,3)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(3)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(3)),'b','Linewidth',1);grid;
xlabel('$x$ [m]','interpreter','latex');ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(7)); end
title( ['VAF = ',num2str(VAF(indices(3)),3), '\% at $k$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
subplot(2,2,4)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(4)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(4)),'b','Linewidth',1);grid;
xlabel('$x$ [m]','interpreter','latex');
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(7)); end
title( ['VAF = ',num2str(VAF(indices(4)),3), '\% at $k$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
if Wp.turbine.N==9
text( -2560, 26, 'First row: WFSim (black) and PALM (blue)','interpreter','latex') ;
%suptitle('First row: WFSim (black) and PALM (blue)')
end

%% Plot cross-sections
% if Wp.turbine.N==2
% indices  = [250 500 750 999];
% 
% for k=indices
%     SOWFAdata     = load([sourcepath num2str(list{k})]);
%     usowfa(:,:,k) = SOWFAdata.uq;
%     VAF(k)        = vaf(usowfa(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,k),...
%     uk(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,k));
% end
% 
% figure(6);clf;
% subplot(2,2,1)
% plot(Wp.mesh.ldyy2(1,:)',uk(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(1)),'k','Linewidth',1);hold on;
% plot(Wp.mesh.ldyy2(1,:)',usowfa(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(1)),'b','Linewidth',1);grid;
% vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(1)));vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(end)))
% ylabel('$u$ [m/s]','interpreter','latex');
% %ylim([3 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
% title( ['VAF = ',num2str(VAF(indices(1)),3), '\% at $k$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
% subplot(2,2,2)
% plot(Wp.mesh.ldyy2(1,:)',uk(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(2)),'k','Linewidth',1);hold on;
% plot(Wp.mesh.ldyy2(1,:)',usowfa(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(2)),'b','Linewidth',1);grid;
% vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(1)));vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(end)))
% %ylim([3 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
% title( ['VAF = ',num2str(VAF(indices(2)),3), '\% at $k$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
% subplot(2,2,3)
% plot(Wp.mesh.ldyy2(1,:)',uk(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(3)),'k','Linewidth',1);hold on;
% plot(Wp.mesh.ldyy2(1,:)',usowfa(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(3)),'b','Linewidth',1);grid;
% vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(1)));vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(end)))
% xlabel('$y$ [m]','interpreter','latex');ylabel('$u$ [m/s]','interpreter','latex');
% %ylim([3 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
% title( ['VAF = ',num2str(VAF(indices(3)),3), '\% at $k$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
% subplot(2,2,4)
% plot(Wp.mesh.ldyy2(1,:)',uk(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(4)),'k','Linewidth',1);hold on;
% plot(Wp.mesh.ldyy2(1,:)',usowfa(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(4)),'b','Linewidth',1);grid;
% vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(1)));vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(end)))
% xlabel('$y$ [m]','interpreter','latex');
% %ylim([3 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
% title( ['VAF = ',num2str(VAF(indices(4)),3), '\% at $k$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
% end

%% Plot control signals (Numbering according topology presented in paper)
if Wp.turbine.N==2
figure(9);clf;
subplot(2,1,1)
plot(Wp.sim.time(1:end-1),a(1,:),'b');hold on;
plot(Wp.sim.time(1:end-1),a(2,:),'k');grid;
ylabel('$a$','interpreter','latex');
title('$a^1$ (blue), $a^2$ (black)','interpreter','latex')
%legend({'$a_1$','$a_2$'},'interpreter','latex')
subplot(2,1,2)
plot(Wp.sim.time(1:end-1),Phi(1,:),'b');hold on;
plot(Wp.sim.time(1:end-1),Phi(2,:),'k');grid;
ylabel('$\gamma$ [deg]','interpreter','latex');xlabel('$k$ [s]','interpreter','latex');
title('$\gamma^1$ (blue), $\gamma^2$ (black)','interpreter','latex')
%legend({'$\gamma_1$','$\gamma_2$'},'interpreter','latex')
end

%% Power
figure(10);clf;
subplot(2,1,1)
plot(Power(1,1:end));hold on;
plot(PowerPALM(1,1:end),'r--');
grid;ylabel('Power T_1');legend('WFSim','PALM')
subplot(2,1,2)
plot(Power(2,1:end));hold on;
plot(PowerPALM(2,1:end),'r--');
grid;xlabel('Time [s]');ylabel('Power T_2')




