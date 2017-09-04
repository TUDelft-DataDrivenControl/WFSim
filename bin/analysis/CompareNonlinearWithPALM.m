clear; clc; close all;

%% Define script settings
% Model settings
scriptOptions.Projection        = 0;        % Use projection (true/false)
scriptOptions.Linearversion     = 0;        % Provide linear variant of WFSim (true/false)
scriptOptions.exportLinearSol   = 0;        % Calculate linear solution of WFSim
scriptOptions.Derivatives       = 0;        % Compute derivatives
scriptOptions.startUniform      = 1;        % Start from a uniform flowfield (true) or a steady-state solution (false)
scriptOptions.exportPressures   = ~scriptOptions.Projection;   % Calculate pressure fields

% Convergence settings
scriptOptions.conv_eps          = 1e-6;     % Convergence threshold
scriptOptions.max_it_dyn        = 1;        % Maximum number of iterations for k > 1
if scriptOptions.startUniform==1
    scriptOptions.max_it = 1;
else
    scriptOptions.max_it = 50;
end

% Display and visualization settings
scriptOptions.printProgress     = 1;  % Print progress every timestep
scriptOptions.printConvergence  = 1;  % Print convergence parameters every timestep
scriptOptions.Animate           = 10;  % Show 2D flow fields every x iterations (0: no plots)
scriptOptions.plotMesh          = 0;  % Show meshing and turbine locations


%%%------------------------------------------------------------------------%%%%

%% Script core
% WFSim: call initialization script
Wp.name      = '2turb_adm_noturb';

run('../../WFSim_addpaths'); % Add paths
[Wp,sol,sys] = InitWFSim(Wp,scriptOptions);

% Initialize variables and figure specific to this script
sol_array = {};
CPUTime   = zeros(Wp.sim.NN,1);
if scriptOptions.Animate > 0
    %scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
end

load(Wp.sim.measurementFile)
sol.u = squeeze(u(1,:,:));
sol.v = squeeze(v(1,:,:));

% Performing timestepping until end
disp(['Performing ' num2str(Wp.sim.NN) ' forward simulations..']);
%% Loop
while sol.k < Wp.sim.NN
    tic;         % Intialize timer
    
    [sol,sys]      = WFSim_timestepping(sol,sys,Wp,scriptOptions); % forward timestep with WFSim
    CPUTime(sol.k) = toc; % Take time
    
    eu                            = vec(sol.u-squeeze(u(sol.k,:,:))); eu(isnan(eu)) = [];
    ev                            = vec(sol.v-squeeze(v(sol.k,:,:))); ev(isnan(ev)) = [];
    RMSE(sol.k)                   = rms([eu;ev]);
    [maxe(sol.k),maxeloc(sol.k)]  = max(abs(eu));
       
    % Save sol to a cell array
    sol_array{sol.k} = sol;
    Power(:,sol.k)   = sol_array{sol.k}.turbine.power;
    
    % Display progress and animations
    if scriptOptions.printProgress
        disp(['Simulated t(' num2str(sol.k) ') = ' num2str(sol.time) ' s. CPU: ' num2str(CPUTime(sol.k)*1e3,3) ' ms.']);
    end;
    
    if scriptOptions.Animate > 0
        if ~rem(sol.k,scriptOptions.Animate)
            
            turb_coord = .5*Wp.turbine.Drotor*exp(1i*turbData.phi(sol.k,:)*pi/180);  % Yaw angles
            
            subplot(2,3,1);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.u,'Linecolor','none');  colormap(hot);
            caxis([min(min(sol.u)) max(max(sol.u))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(turb_coord(ll))):1:(Wp.turbine.Cry(ll)+real(turb_coord(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(turb_coord(ll)),Wp.turbine.Crx(ll)+imag(turb_coord(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            text(0,Wp.mesh.ldxx2(end,end)+250,['Time ', num2str(Wp.sim.time(sol.k),'%.1f'), 's']);
            ylabel('x [m]');
            title('WFSim u [m/s]');
            hold off;
            subplot(2,3,2);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',squeeze(u(sol.k,:,:)),'Linecolor','none');  colormap(hot);
            caxis([min(min(squeeze(u(sol.k,:,:)))) max(max(squeeze(u(sol.k,:,:))))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(turb_coord(ll))):1:(Wp.turbine.Cry(ll)+real(turb_coord(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(turb_coord(ll)),Wp.turbine.Crx(ll)+imag(turb_coord(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            title('PALM u [m/s]');
            hold off;
            
            subplot(2,3,3);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.u-squeeze(u(sol.k,:,:)),'Linecolor','none');  colormap(hot);
            caxis([min(min(sol.u-squeeze(u(sol.k,:,:)))) max(max(sol.u-squeeze(u(sol.k,:,:))))]);  hold all; colorbar;
            axis equal; axis tight;
            ldyyv = Wp.mesh.ldyy(:); ldxx2v = Wp.mesh.ldxx2(:);
            plot(ldyyv(maxeloc(sol.k)),ldxx2v(maxeloc(sol.k)),'whiteo','LineWidth',1,'MarkerSize',8,'DisplayName','Maximum error location');
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(turb_coord(ll))):1:(Wp.turbine.Cry(ll)+real(turb_coord(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(turb_coord(ll)),Wp.turbine.Crx(ll)+imag(turb_coord(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            title('error [m/s]');
            hold off;           
            subplot(2,3,4);
            plot(Wp.sim.time(1:sol.k),turbData.power(1:sol.k,1));hold on
            plot(Wp.sim.time(1:sol.k),turbData.power(1:sol.k,2),'r');
            title('$P$ [W]','interpreter','latex');
            axis([0,Wp.sim.NN 0 max(max(turbData.power(1:end,:)))+10^5])
            title('Power PALM')
            grid;hold off;
            drawnow
            
            subplot(2,3,5);
            plot(Wp.sim.time(1:sol.k),Power(1,1:sol.k));hold on
            plot(Wp.sim.time(1:sol.k),Power(2,1:sol.k),'r');
            title('$P$ [W]','interpreter','latex');
            axis([0,Wp.sim.NN 0 max(max(turbData.power(1:end,1)))+10^5]);
            title('Power WFSim')
            grid;hold off;
            drawnow
            
            subplot(2,3,6)
            plot(Wp.sim.time(1:sol.k),Power(1,1:sol.k)'-turbData.power(1:sol.k,1));hold on
            plot(Wp.sim.time(1:sol.k),Power(2,1:sol.k)'-turbData.power(1:sol.k,2),'r');
            title('error [W]','interpreter','latex');
            axis([0,Wp.sim.NN min(min(Power(:,1:sol.k)'-turbData.power(1:sol.k,:)))-.2 max(max(Power(:,1:sol.k)'-turbData.power(1:sol.k,:)))+.2]);
            grid;hold off;
        end
        drawnow
        
    end;

end;
turbData.power = turbData.power';
%%
figure(2);clf
plot(Wp.sim.time(1:Wp.sim.NN),RMSE);hold on;
plot(Wp.sim.time(1:Wp.sim.NN),maxe,'r');grid;
ylabel('RMSE and max');
title(['{\color{blue}{RMSE}}, {\color{red}{max}} and meanRMSE = ',num2str(mean(RMSE),3)])

Nt = Wp.sim.NN; % until what time you want to plot

figure(3);clf
plot(Wp.sim.time(1:Nt),sum(Power(:,1:Nt)),'k','Linewidth',1);hold on;
plot(Wp.sim.time(1:Nt),sum(turbData.power(:,1:Nt)),'b--');
grid;xlabel('$t [s]$','interpreter','latex');
ylabel('$P$ [W]','interpreter','latex');
title('Wind farm power: WFSim (black) PALM (blue dashed)','interpreter','latex');
xlim([0 Nt])
if Wp.turbine.N==2
    figure(4);clf
    subplot(2,1,1)
    plot(Wp.sim.time(1:Nt),turbData.CT(1:Nt,1),'linewidth',1.5);hold on
    plot(Wp.sim.time(1:Nt),turbData.CT(1:Nt,2),'r--','linewidth',1.5);
    ylabel('$CT^{\prime}$','interpreter','latex');
    xlabel('$t [s]$','interpreter','latex');
    title('$T_1$ (blue), $T_2$ (red dashed) ','interpreter','latex');
    axis([0,Wp.sim.time(Nt) 0 max(max(turbData.CT(1:Nt,:)))+.2]);
    grid;hold off;xlim([0 Nt]);
    subplot(2,1,2)
    plot(Wp.sim.time(1:Nt),turbData.phi(1:Nt,1));hold on
    plot(Wp.sim.time(1:Nt),turbData.phi(1:Nt,1),'r');
    ylabel('$\gamma$','interpreter','latex');
    xlabel('$t [s]$','interpreter','latex');
    title('$T_1$ (blue), $T_2$ (red) ','interpreter','latex');
    axis([0,Wp.sim.time(Nt) min(min(turbData.phi(1:Nt,:)))-5 max(max(turbData.phi(1:Nt,:)))+5]);
    grid;hold off;  
elseif Wp.turbine.N==9
    figure(4);clf;
    subplot(1,3,1)
    plot(Wp.sim.time(1:Nt),turbData.CT(1:Nt,:),'b');hold on;
    plot(Wp.sim.time(1:Nt),turbData.CT(2:Nt,:),'k');
    plot(Wp.sim.time(1:Nt),turbData.CT(3:Nt,:),'r');grid;
    xlim([0 Wp.sim.time(Nt)])
    ylabel('$CT^{\prime}$','interpreter','latex');
    xlabel('$t$ [s]','interpreter','latex');
    title('$CT^{\prime}_1$ (blue), $CT^{\prime}_2$ (black), $CT^{\prime}_3$ (red)','interpreter','latex')
    subplot(1,3,2)
    plot(Wp.sim.time(1:Nt),turbData.CT(4:Nt,:),'b');hold on;
    plot(Wp.sim.time(1:Nt),turbData.CT(5:Nt,:),'k');
    plot(Wp.sim.time(1:Nt),turbData.CT(6:Nt,:),'r');grid;
    xlim([0 Wp.sim.time(Nt)])
    xlabel('$t$ [s]','interpreter','latex');
    title('$CT^{\prime}_4$ (blue), $CT^{\prime}_5$ (black), $CT^{\prime}_6$ (red)','interpreter','latex')
    subplot(1,3,3)
    plot(Wp.sim.time(1:Nt),turbData.CT(7:Nt,:),'b');hold on;
    plot(Wp.sim.time(1:Nt),turbData.CT(8:Nt,:),'k');
    plot(Wp.sim.time(1:Nt),turbData.CT(9:Nt,:),'r');grid;
    xlim([0 Wp.sim.time(Nt)])
    title('$CT^{\prime}_7$ (blue), $CT^{\prime}_8$ (black), $CT^{\prime}_9$ (red)','interpreter','latex')
    xlabel('$t$ [s]','interpreter','latex');
end

% Wake centreline
D_ind    = Wp.mesh.yline{1};
indices  = [300 400 700 800];
%indices  = [200 300 400 500];

for k=indices
    up(:,k)      = mean(sol_array{k}.u(:,D_ind),2); 
    uPALM        = squeeze(u(k,:,:));   
    upPALM(:,k)  = mean(uPALM(:,D_ind),2);
    VAF(:,k)     = vaf(upPALM(:,k),up(:,k));
end

figure(5);clf;
subplot(2,2,1)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(1)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(1)),'b--','Linewidth',1);grid;
ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));
if Wp.turbine.N>1
vline(Wp.turbine.Crx(2));
end
if Wp.turbine.N==9; vline(Wp.turbine.Crx(3)); end
title( ['VAF = ',num2str(VAF(indices(1)),3), '\% at $t$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
subplot(2,2,2)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(2)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(2)),'b--','Linewidth',1);grid;
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));
if Wp.turbine.N>1
vline(Wp.turbine.Crx(2));
end
if Wp.turbine.N==9; vline(Wp.turbine.Crx(3)); end
title( ['VAF = ',num2str(VAF(indices(2)),3), '\% at $t$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
subplot(2,2,3)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(3)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(3)),'b--','Linewidth',1);grid;
xlabel('$x$ [m]','interpreter','latex');ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));
if Wp.turbine.N>1
vline(Wp.turbine.Crx(2));
end
if Wp.turbine.N==9; vline(Wp.turbine.Crx(3)); end
title( ['VAF = ',num2str(VAF(indices(3)),3), '\% at $t$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
subplot(2,2,4)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(4)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upPALM(:,indices(4)),'b--','Linewidth',1);grid;
xlabel('$x$ [m]','interpreter','latex');
ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));
if Wp.turbine.N>1
    vline(Wp.turbine.Crx(2));
end
if Wp.turbine.N>2
    vline(Wp.turbine.Crx(3));
end
title( ['VAF = ',num2str(VAF(indices(4)),3), '\% at $t$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
if Wp.turbine.N==2
    text( -1550, 20.4, 'WFSim (black) and PALM (blue dashed)','interpreter','latex') ;
    %suptitle('First row: WFSim (black) and PALM (blue)')
end
if Wp.turbine.N==9
    text( -1800, 20.4, 'First row: WFSim (black) and PALM (blue dashed)','interpreter','latex') ;
    %suptitle('First row: WFSim (black) and PALM (blue)')
end
%
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
        up(:,k)      = mean(sol_array{k}.u(:,D_ind),2);
        uPALM        = squeeze(u(k,:,:));
        upPALM(:,k)  = mean(uPALM(:,D_ind),2);
        VAF_2(:,k)   = vaf(upPALM(:,k),up(:,k));
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
    up(:,k)      = mean(sol_array{k}.u(:,D_ind),2); 
    uPALM        = squeeze(u(k,:,:));   
    upPALM(:,k)  = mean(uPALM(:,D_ind),2);
    VAF_3(:,k)   = vaf(upPALM(:,k),up(:,k));
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
    plot(turbData.power(1,1:Nt),'b--');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    grid;ylabel('$P_1$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(Nt)])
    subplot(3,3,2)
    plot(Power(2,1:end),'k','Linewidth',1);hold on;
    plot(turbData.power(2,1:Nt),'b--');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    grid;ylabel('$P_2$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(Nt)])
    title('Power: WFSim (black) PALM (dashed blue)','interpreter','latex')
    subplot(3,3,3)
    plot(Power(3,1:end),'k','Linewidth',1);hold on;
    plot(turbData.power(3,1:Nt),'b--');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    grid;ylabel('$P_3$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(Nt)])
    subplot(3,3,4)
    plot(Power(4,1:end),'k','Linewidth',1);hold on;
    plot(turbData.power(4,1:Nt),'b--');
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    grid;ylabel('$P_4$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(Nt)])
    subplot(3,3,5)
    plot(Power(5,1:end),'k','Linewidth',1);hold on;
    plot(turbData.power(5,1:Nt),'b--');
    grid;ylabel('$P_5$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(Nt)])
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    subplot(3,3,6)
    plot(Power(6,1:end),'k','Linewidth',1);hold on;
    plot(turbData.power(6,1:Nt),'b--');
    grid;ylabel('$P_6$','interpreter','latex');ylim([0 n*10^6]);xlim([0 Wp.sim.time(Nt)])
    set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
    subplot(3,3,7)
    plot(Power(7,1:end),'k','Linewidth',1);hold on;
    plot(turbData.power(7,1:Nt),'b--');
    grid;xlabel('$t [s]$','interpreter','latex');ylabel('$P_7$','interpreter','latex');ylim([0 n*10^6])
    xlim([0 Wp.sim.time(Nt)])
    subplot(3,3,8)
    plot(Power(8,1:end),'k','Linewidth',1);hold on;
    plot(turbData.power(8,1:Nt),'b--');
    grid;xlabel('$t [s]$','interpreter','latex');ylabel('$P_8$','interpreter','latex');ylim([0 n*10^6]);
    xlim([0 Wp.sim.time(Nt)])
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    subplot(3,3,9)
    plot(Power(9,1:end),'k','Linewidth',1);hold on;
    plot(turbData.power(9,1:Nt),'b--');
    grid;xlabel('$t [s]$','interpreter','latex');ylabel('$P_9$','interpreter','latex');ylim([0 n*10^6]);
    set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
    xlim([0 Wp.sim.time(Nt)])
end

% Individual powers
figure(101);
subplot(2,1,1)
plot(Wp.sim.time(1:Nt),turbData.power(1,1:Nt),'b--');hold on
plot(Wp.sim.time(1:Nt),Power(1,1:Nt),'k');
axis([0,Wp.sim.time(Nt) 0 max(max(turbData.power(:,1:end)))+10^5])
ylabel('$P_1$ [W]', 'interpreter','latex')
xlabel('$t [s]$','interpreter','latex');
title('$T_1$: WFSim (black) PALM (blue dashed)', 'interpreter','latex')
grid;hold off;

subplot(2,1,2)
plot(Wp.sim.time(1:Nt),turbData.power(2,1:Nt),'b--');hold on
plot(Wp.sim.time(1:Nt),Power(2,1:Nt),'k');
axis([0,Wp.sim.time(Nt) 0 max(max(turbData.power(:,1:end)))+10^5])
ylabel('$P_2$ [W]', 'interpreter','latex')
xlabel('$t [s]$','interpreter','latex');
title('$T_2$: WFSim (black) PALM (blue dashed)', 'interpreter','latex')
grid;hold off;

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
