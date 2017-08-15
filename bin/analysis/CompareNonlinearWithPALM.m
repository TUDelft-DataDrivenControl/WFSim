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

Wp.name             = '1turb_adm';
Wp.Turbulencemodel  = 'WFSim3';

if strcmp('9turb_adm',Wp.name) || strcmp('9turb_adm_turbulence',Wp.name)
    options.startUniform  = 0;
end

Animate       = 50;                      % Show 2D flow fields every x iterations (0: no plots)
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
if Wp.turbine.N==1
    M1  = load('../../Data_PALM/1turb_adm/example_1turb_adm_matlab_turbine_parameters01.txt');
    
    filename = '../../Data_PALM/1turb_adm/example_1turb_adm_matlab_m01.nc';
elseif Wp.turbine.N==2
    M1  = load('../../Data_PALM/2turb_adm/example_2turb_adm_matlab_turbine_parameters01.txt');
    M2  = load('../../Data_PALM/2turb_adm/example_2turb_adm_matlab_turbine_parameters02.txt');
    
    filename = '../../Data_PALM/2turb_adm/example_2turb_adm_matlab_m01.nc';
elseif Wp.turbine.N==9 && strcmp('9turb_adm',Wp.name)
    M1  = load('../../Data_PALM/9turb_adm/example_9turb_adm_matlab_turbine_parameters01.txt');
    M2  = load('../../Data_PALM/9turb_adm/example_9turb_adm_matlab_turbine_parameters02.txt');
    M3  = load('../../Data_PALM/9turb_adm/example_9turb_adm_matlab_turbine_parameters03.txt');
    M4  = load('../../Data_PALM/9turb_adm/example_9turb_adm_matlab_turbine_parameters04.txt');
    M5  = load('../../Data_PALM/9turb_adm/example_9turb_adm_matlab_turbine_parameters05.txt');
    M6  = load('../../Data_PALM/9turb_adm/example_9turb_adm_matlab_turbine_parameters06.txt');
    M7  = load('../../Data_PALM/9turb_adm/example_9turb_adm_matlab_turbine_parameters07.txt');
    M8  = load('../../Data_PALM/9turb_adm/example_9turb_adm_matlab_turbine_parameters08.txt');
    M9  = load('../../Data_PALM/9turb_adm/example_9turb_adm_matlab_turbine_parameters09.txt');
    
    filename = '../../Data_PALM/9turb_adm/example_9turb_adm_matlab_m01.nc';
elseif Wp.turbine.N==9 && strcmp('9turb_adm_turbulence',Wp.name)
    M1  = load('../../Data_PALM/9turb_adm_turbulence/example_9turb_adm_matlab_turbine_parameters01.txt');
    M2  = load('../../Data_PALM/9turb_adm_turbulence/example_9turb_adm_matlab_turbine_parameters02.txt');
    M3  = load('../../Data_PALM/9turb_adm_turbulence/example_9turb_adm_matlab_turbine_parameters03.txt');
    M4  = load('../../Data_PALM/9turb_adm_turbulence/example_9turb_adm_matlab_turbine_parameters04.txt');
    M5  = load('../../Data_PALM/9turb_adm_turbulence/example_9turb_adm_matlab_turbine_parameters05.txt');
    M6  = load('../../Data_PALM/9turb_adm_turbulence/example_9turb_adm_matlab_turbine_parameters06.txt');
    M7  = load('../../Data_PALM/9turb_adm_turbulence/example_9turb_adm_matlab_turbine_parameters07.txt');
    M8  = load('../../Data_PALM/9turb_adm_turbulence/example_9turb_adm_matlab_turbine_parameters08.txt');
    M9  = load('../../Data_PALM/9turb_adm_turbulence/example_9turb_adm_matlab_turbine_parameters09.txt');
    
    filename = '../../Data_PALM/9turb_adm_turbulence/example_9turb_adm_matlab_m01.nc';
end

%M = [Time   UR  Uinf  Ct_adm  a Yaw Thrust Power  WFPower]

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
    
    uPALM                 = reshape(u(k,nz,:,:),size(u,3),size(u,4))';  % u(k,z,y,x)
    vPALM                 = reshape(v(k,nz,:,:),size(v,3),size(v,4))';
    % Interpolate PALM data on WFSim grid
    targetSize            = [Wp.mesh.Nx Wp.mesh.Ny];
    sourceSize            = size(uPALM);
    [X_samples,Y_samples] = meshgrid(linspace(1,sourceSize(2),targetSize(2)), linspace(1,sourceSize(1),targetSize(1)));
    uPALM                 = interp2(uPALM, X_samples, Y_samples);
    vPALM                 = interp2(vPALM, X_samples, Y_samples);
    
    
    %if k==2
    %    sol.u = uPALM;
    %    sol.v = vPALM;
    %end
    if k==300
        
    end
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
    
    eu                    = vec(sol.u-uPALM); eu(isnan(eu)) = [];
    ev                    = vec(sol.v-vPALM); ev(isnan(ev)) = [];
    RMSE(k)               = rms([eu;ev]);
    [maxe(k),maxeloc(k)]  = max(abs(eu));
    
    if Wp.turbine.N==1
        PowerPALM(:,k)       = M1(k,8);
    elseif Wp.turbine.N==2
        PowerPALM(:,k)       = [M1(k,8);M2(k,8)];
    elseif Wp.turbine.N==9
        PowerPALM(:,k)       = [M1(k,8);M2(k,8);M3(k,8);M4(k,8);M5(k,8);M6(k,8);M7(k,8);M8(k,8);M9(k,8)];
    end
    
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
            title('$T_1$: blue PALM, red WFSim', 'interpreter','latex')
            grid;hold off;
            drawnow
            
            if Wp.turbine.N>1
                subplot(2,3,5);
                plot(Wp.sim.time(1:k),PowerPALM(2,1:k));hold on
                plot(Wp.sim.time(1:k),Power(2,1:k),'r');
                title('$P$ [W]','interpreter','latex');
                axis([0,Wp.sim.time(size(u,1)) 0 max(max(Power(:,1:end)))+10^5]);
                title('$T_2$: blue PALM, red WFSim', 'interpreter','latex')
                grid;hold off;
                drawnow
                
                subplot(2,3,6)
                plot(Wp.sim.time(1:k),sum(PowerPALM(:,1:k)),'b--');hold on
                plot(Wp.sim.time(1:k),sum(Power(:,1:k)),'k');
                title('WF power [W]','interpreter','latex');
                axis([0,Wp.sim.time(size(u,1)) min(min(sum(PowerPALM(:,1:k))))-10^3 max(max(sum(PowerPALM(:,1:k))))+10^3]);
                grid;hold off;
            end
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
plot(Wp.sim.time(1:size(u,1)),sum(Power(:,1:size(u,1))),'k','Linewidth',1);hold on;
plot(Wp.sim.time(1:size(u,1)),sum(PowerPALM(:,1:size(u,1))),'b--');
grid;xlabel('$t [s]$','interpreter','latex');
ylabel('$\Sigma P$ [W]','interpreter','latex');
title('WF Power: PALM (blue dashed), WFSim (black) ','interpreter','latex');

if Wp.turbine.N==2
    figure(4);clf
    plot(Wp.sim.time(1:size(u,1)),CT(1,1:size(u,1)));hold on
    plot(Wp.sim.time(1:size(u,1)),CT(2,1:size(u,1)),'r');
    ylabel('$CT^{\prime}$','interpreter','latex');
    xlabel('$t [s]$','interpreter','latex');
    title('$T_1$ (blue), $T_2$ (red) ','interpreter','latex');
    axis([0,Wp.sim.time(size(u,1)) 0 max(max(CT(:,1:size(u,1))))+.2]);
    grid;hold off;
elseif Wp.turbine.N==9
    figure(4);clf;
    subplot(1,3,1)
    plot(Wp.sim.time(1:size(u,1)),CT(1,1:size(u,1)),'b');hold on;
    plot(Wp.sim.time(1:size(u,1)),CT(2,1:size(u,1)),'k');
    plot(Wp.sim.time(1:size(u,1)),CT(3,1:size(u,1)),'r');grid;
    xlim([0 Wp.sim.time(size(u,1))])
    ylabel('$CT^{\prime}$','interpreter','latex');
    xlabel('$t$ [s]','interpreter','latex');
    title('$CT^{\prime}_1$ (blue), $CT^{\prime}_2$ (black), $CT^{\prime}_3$ (red)','interpreter','latex')
    subplot(1,3,2)
    plot(Wp.sim.time(1:size(u,1)),CT(4,1:size(u,1)),'b');hold on;
    plot(Wp.sim.time(1:size(u,1)),CT(5,1:size(u,1)),'k');
    plot(Wp.sim.time(1:size(u,1)),CT(6,1:size(u,1)),'r');grid;
    xlim([0 Wp.sim.time(size(u,1))])
    xlabel('$t$ [s]','interpreter','latex');
    title('$CT^{\prime}_4$ (blue), $CT^{\prime}_5$ (black), $CT^{\prime}_6$ (red)','interpreter','latex')
    subplot(1,3,3)
    plot(Wp.sim.time(1:size(u,1)),CT(7,1:size(u,1)),'b');hold on;
    plot(Wp.sim.time(1:size(u,1)),CT(8,1:size(u,1)),'k');
    plot(Wp.sim.time(1:size(u,1)),CT(9,1:size(u,1)),'r');grid;
    xlim([0 Wp.sim.time(size(u,1))])
    title('$CT^{\prime}_7$ (blue), $CT^{\prime}_8$ (black), $CT^{\prime}_9$ (red)','interpreter','latex')
    xlabel('$t$ [s]','interpreter','latex');
end

% Wake centreline
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


figure(5);clf;
subplot(2,2,1)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(1)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(1)),'b--','Linewidth',1);grid;
ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(3)); end
title( ['VAF = ',num2str(VAF(indices(1)),3), '\% at $t$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
subplot(2,2,2)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(2)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(2)),'b--','Linewidth',1);grid;
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(3)); end
title( ['VAF = ',num2str(VAF(indices(2)),3), '\% at $t$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
subplot(2,2,3)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(3)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(3)),'b--','Linewidth',1);grid;
xlabel('$x$ [m]','interpreter','latex');ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(3)); end
title( ['VAF = ',num2str(VAF(indices(3)),3), '\% at $t$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
subplot(2,2,4)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(4)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(4)),'b--','Linewidth',1);grid;
xlabel('$x$ [m]','interpreter','latex');
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(3));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(7)); end
title( ['VAF = ',num2str(VAF(indices(4)),3), '\% at $t$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
if Wp.turbine.N==2
    text( -1300, 20.4, 'WFSim (black) and PALM (blue)','interpreter','latex') ;
    %suptitle('First row: WFSim (black) and PALM (blue)')
end
if Wp.turbine.N==9
    text( -1800, 20.4, 'First row: WFSim (black) and PALM (blue)','interpreter','latex') ;
    %suptitle('First row: WFSim (black) and PALM (blue)')
end

if Wp.turbine.N==9
    
    yline    = Wp.mesh.yline{4};
    m        = size(yline,2);
    if rem(m,2)
        ind  = ceil(m/2);
    else
        ind  = [m/2 m/2+1];
    end
    yline    = yline(ind);
    D_ind    = yline;
    
    clear up upPALM
    for k=indices
        up(:,k)       = mean(uk(:,D_ind,k),2);
        uPALM         = reshape(u(k,nz,:,:),size(u,3),size(u,4))';
        % Interpolate PALM data on WFSim grid
        targetSize            = [Wp.mesh.Nx Wp.mesh.Ny];
        sourceSize            = size(uPALM);
        [X_samples,Y_samples] = meshgrid(linspace(1,sourceSize(2),targetSize(2)), linspace(1,sourceSize(1),targetSize(1)));
        uPALM                 = interp2(uPALM, X_samples, Y_samples);
        upPALM(:,k)  = mean(uPALM(:,D_ind),2);
        VAF_2(:,k)     = vaf(upPALM(:,k),up(:,k));
    end
    
    figure(6);clf;
    subplot(2,2,1)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(1)),'k','Linewidth',1.5);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(1)),'b--','Linewidth',1);grid;
    ylabel('$U^c$ [m/s]','interpreter','latex');
    ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(3));
    title( ['VAF = ',num2str(VAF_2(indices(1)),3), '\% at $t$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,2)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(2)),'k','Linewidth',1.5);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(2)),'b--','Linewidth',1);grid;
    ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(3));
    title( ['VAF = ',num2str(VAF_2(indices(2)),3), '\% at $t$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,3)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(3)),'k','Linewidth',1.5);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(3)),'b--','Linewidth',1);grid;
    xlabel('$x$ [m]','interpreter','latex');ylabel('$U^c$ [m/s]','interpreter','latex');
    ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(3));
    title( ['VAF = ',num2str(VAF_2(indices(3)),3), '\% at $t$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,4)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(4)),'k','Linewidth',1.5);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(4)),'b--','Linewidth',1);grid;
    xlabel('$x$ [m]','interpreter','latex');
    ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(3));
    title( ['VAF = ',num2str(VAF_2(indices(4)),3), '\% at $t$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
    text( -1800, 20.4, 'Second row: WFSim (black) and PALM (blue dashed)','interpreter','latex') ;
    %suptitle('Second row: WFSim (black) and PALM (blue)')
    
    yline    = Wp.mesh.yline{7};
    m        = size(yline,2);
    if rem(m,2)
        ind  = ceil(m/2);
    else
        ind  = [m/2 m/2+1];
    end
    yline    = yline(ind);
    D_ind    = yline;
    
    clear up upPALM
    for k=indices
        up(:,k)               = mean(uk(:,D_ind,k),2);
        uPALM                 = reshape(u(k,nz,:,:),size(u,3),size(u,4))';
        % Interpolate PALM data on WFSim grid
        targetSize            = [Wp.mesh.Nx Wp.mesh.Ny];
        sourceSize            = size(uPALM);
        [X_samples,Y_samples] = meshgrid(linspace(1,sourceSize(2),targetSize(2)), linspace(1,sourceSize(1),targetSize(1)));
        uPALM                 = interp2(uPALM, X_samples, Y_samples);
        upPALM(:,k)           = mean(uPALM(:,D_ind),2);
        VAF_3(:,k)            = vaf(upPALM(:,k),up(:,k));
    end
    
    figure(7);clf;
    subplot(2,2,1)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(1)),'k','Linewidth',1.5);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(1)),'b--','Linewidth',1);grid;
    ylabel('$U^c$ [m/s]','interpreter','latex');
    ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));vline(Wp.turbine.Crx(3));
    title( ['VAF = ',num2str(VAF_3(indices(1)),3), '\% at $t$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,2)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(2)),'k','Linewidth',1.5);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(2)),'b--','Linewidth',1);grid;
    ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));vline(Wp.turbine.Crx(3));
    title( ['VAF = ',num2str(VAF_3(indices(2)),3), '\% at $t$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,3)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(3)),'k','Linewidth',1.5);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(3)),'b--','Linewidth',1);grid;
    xlabel('$x$ [m]','interpreter','latex');ylabel('$U^c$ [m/s]','interpreter','latex');
    ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));vline(Wp.turbine.Crx(3));
    title( ['VAF = ',num2str(VAF_3(indices(3)),3), '\% at $t$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,4)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(4)),'k','Linewidth',1.5);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(4)),'b--','Linewidth',1);grid;
    xlabel('$x$ [m]','interpreter','latex');
    ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));vline(Wp.turbine.Crx(3));
    title( ['VAF = ',num2str(VAF_3(indices(4)),3), '\% at $t$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
    %suptitle('Third row: WFSim (black) and PALM (blue)')
    text( -1800, 20.4, 'Third row: WFSim (black) and PALM (blue dashed)','interpreter','latex') ;
    
    n = 2;
    figure(10);clf;
    subplot(3,3,1)
    plot(Power(1,1:end),'k','Linewidth',1);hold on;
    plot(PowerPALM(1,1:end),'b--');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    grid;ylabel('$P_1$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(size(u,1))])
    subplot(3,3,2)
    plot(Power(2,1:end),'k','Linewidth',1);hold on;
    plot(PowerPALM(2,1:end),'b--');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    grid;ylabel('$P_2$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(size(u,1))])
    subplot(3,3,3)
    plot(Power(3,1:end),'k','Linewidth',1);hold on;
    plot(PowerPALM(3,1:end),'b--');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    grid;ylabel('$P_3$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(size(u,1))])
    subplot(3,3,4)
    plot(Power(4,1:end),'k','Linewidth',1);hold on;
    plot(PowerPALM(4,1:end),'b--');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    grid;ylabel('$P_4$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(size(u,1))])
    subplot(3,3,5)
    plot(Power(5,1:end),'k','Linewidth',1);hold on;
    plot(PowerPALM(5,1:end),'b--');
    grid;ylabel('$P_5$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(size(u,1))])
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    subplot(3,3,6)
    plot(Power(6,1:end),'k','Linewidth',1);hold on;
    plot(PowerPALM(6,1:end),'b--');
    grid;ylabel('$P_6$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(size(u,1))])
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    subplot(3,3,7)
    plot(Power(7,1:end),'k','Linewidth',1);hold on;
    plot(PowerPALM(7,1:end),'b--');
    grid;xlabel('$t [s]$','interpreter','latex');ylabel('$P_7$','interpreter','latex');ylim([0 n*10^6])
    xlim([0 Wp.sim.time(size(u,1))])
    subplot(3,3,8)
    plot(Power(8,1:end),'k','Linewidth',1);hold on;
    plot(PowerPALM(8,1:end),'b--');
    grid;xlabel('$t [s]$','interpreter','latex');ylabel('$P_8$','interpreter','latex');ylim([0 n*10^6]);
    xlim([0 Wp.sim.time(size(u,1))])
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    subplot(3,3,9)
    plot(Power(9,1:end),'k','Linewidth',1);hold on;
    plot(PowerPALM(9,1:end),'b--');
    grid;xlabel('$t [s]$','interpreter','latex');ylabel('$P_9$','interpreter','latex');ylim([0 n*10^6]);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    xlim([0 Wp.sim.time(size(u,1))])
end




