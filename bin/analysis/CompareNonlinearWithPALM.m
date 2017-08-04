clear; clc; close all;

run('../../WFSim_addpaths');

[WFSimFolder, ~, ~] = fileparts(which([mfilename '.m']));   % Get WFSim directory

% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 1;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name             = '2turb_adm';
Wp.Turbulencemodel  = 'WFSim3';

Animate       = 100;                      % Show 2D flow fields every x iterations (0: no plots)
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

% Load PALM data
M1  = load('../../Data_PALM/2turb_adm/example_2turb_adm_matlab_turbine_parameters01.txt');
M2  = load('../../Data_PALM/2turb_adm/example_2turb_adm_matlab_turbine_parameters02.txt');
%M = [Time   UR  Uinf  Ct_adm  a Yaw Thrust Power  WFPower]

filename = '../../Data_PALM/2turb_adm/example_2turb_adm_matlab_m01.nc';
u        = double(nc_varget(filename,'u'));
v        = double(nc_varget(filename,'v'));
x        = double(nc_varget(filename,'x'));
y        = double(nc_varget(filename,'y'));
xu       = double(nc_varget(filename,'xu'));
yv       = double(nc_varget(filename,'yv'));
zw_3d    = double(nc_varget(filename,'zw_3d'));
nz       = 4;

%% Loop
for k=1:size(u,1)
    
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
    
    uPALM                 = reshape(u(k,nz,:,:),size(u,3),size(u,4))';  % u(k,z,y,x)
    vPALM                 = reshape(v(k,nz,:,:),size(v,3),size(v,4))';
    % Interpolate PALM data on WFSim grid
    targetSize            = [Wp.mesh.Nx Wp.mesh.Ny];
    sourceSize            = size(uPALM);
    [X_samples,Y_samples] = meshgrid(linspace(1,sourceSize(2),targetSize(2)), linspace(1,sourceSize(1),targetSize(1)));
    uPALM                 = interp2(uPALM, X_samples, Y_samples);
    vPALM                 = interp2(vPALM, X_samples, Y_samples);
    
    eu                    = vec(sol.u-uPALM); eu(isnan(eu)) = [];
    ev                    = vec(sol.v-vPALM); ev(isnan(ev)) = [];
    RMSE(k)               = rms([eu;ev]);
    [maxe(k),maxeloc(k)]  = max(abs(eu));
    
    PowerPALM(:,k)       = [M1(k,8);M2(k,8)];
    
    if Animate > 0
        if ~rem(k,Animate)
            
            yaw_angles = .5*Wp.turbine.Drotor*exp(1i*input{k}.phi*pi/180);  % Yaw angles
            
            subplot(2,3,1);
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
            
            subplot(2,3,2);
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
            
            subplot(2,3,3);
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
            
            subplot(2,3,4);
            plot(Wp.sim.time(1:k),PowerPALM(1,1:k));hold on
            plot(Wp.sim.time(1:k),Power(1,1:k),'r');
            title('$P$ [W]','interpreter','latex');
            axis([0,Wp.sim.time(size(u,1)) 0 max(max(PowerPALM(:,1:end)))+10^5])
            title('upwind: blue PALM, red WFSim')
            grid;hold off;
            drawnow
 
            subplot(2,3,5);
            plot(Wp.sim.time(1:k),PowerPALM(1,1:k));hold on
            plot(Wp.sim.time(1:k),Power(2,1:k),'r');
            title('$P$ [W]','interpreter','latex');
            axis([0,Wp.sim.time(size(u,1)) 0 max(max(Power(:,1:end)))+10^5]);
            title('downwind: blue PALM, red WFSim')
            grid;hold off;
            drawnow
            
            subplot(2,3,6)
            plot(Wp.sim.time(1:k),PowerPALM(1,1:k)+PowerPALM(2,1:k));hold on
            plot(Wp.sim.time(1:k),Power(1,1:k)+Power(2,1:k),'r');
            title('WF power [W]','interpreter','latex');
            axis([0,Wp.sim.time(size(u,1)) min(min(Power(:,1:k)+Power(:,1:k)))-10^3 max(max(Power(:,1:k)+Power(:,1:k)))-10^3]);
            grid;hold off;
        end
        drawnow
               
    end;
end;

%%
figure(2);clf
plot(Wp.sim.time(1:size(u,1)),RMSE);hold on;
plot(Wp.sim.time(1:size(u,1)),maxe,'r');grid;
ylabel('RMSE and max');
title(['{\color{blue}{RMSE}}, {\color{red}{max}} and meanRMSE = ',num2str(mean(RMSE),3)])

figure(3);clf            
plot(Wp.sim.time(1:size(u,1)),a(1,1:size(u,1)));hold on
plot(Wp.sim.time(1:size(u,1)),a(2,1:size(u,1)),'r');
ylabel('$a$','interpreter','latex')
title('axial induction','interpreter','latex');
axis([0,Wp.sim.time(size(u,1)) 0 max(max(a(:,1:size(u,1))))+.2]);
grid;hold off;
            
            
%% Wake centreline
D_ind    = Wp.mesh.yline{1};
indices  = [300 400 700 800];

for k=indices
    up(:,k)      = mean(uk(:,D_ind,k),2);
    uPALM        = reshape(u(k,nz,:,:),size(u,3),size(u,4))';
    
    % Interpolate PALM data on WFSim grid
    targetSize            = [Wp.mesh.Nx Wp.mesh.Ny];
    sourceSize            = size(uPALM);
    [X_samples,Y_samples] = meshgrid(linspace(1,sourceSize(2),targetSize(2)), linspace(1,sourceSize(1),targetSize(1)));
    uPALM                 = interp2(uPALM, X_samples, Y_samples);
    
    upPALM(:,k)  = mean(uPALM(:,D_ind),2);
    VAF(:,k)     = vaf(upPALM(:,k),up(:,k));
end


figure(4);clf;
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

