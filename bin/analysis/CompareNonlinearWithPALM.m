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
Wp.name      = '6turb_adm_turb';

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
sol.v = 0*squeeze(v(1,:,:));

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
                text(Wp.turbine.Cry(ll)+50,Wp.turbine.Crx(ll),num2str(ll))
            end
            text(0,Wp.mesh.ldxx2(end,end)+250,['Time ', num2str(Wp.sim.time(sol.k),'%.1f'), 's']);
            ylabel('$x$ [m]','interpreter','latex');
            title('WFSim $u$ [m/s]','interpreter','latex');
            hold off;
            subplot(2,3,2);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',squeeze(u(sol.k,:,:)),'Linecolor','none');  colormap(hot);
            caxis([min(min(squeeze(u(sol.k,:,:)))) max(max(squeeze(u(sol.k,:,:))))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(turb_coord(ll))):1:(Wp.turbine.Cry(ll)+real(turb_coord(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(turb_coord(ll)),Wp.turbine.Crx(ll)+imag(turb_coord(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
                text(Wp.turbine.Cry(ll)+50,Wp.turbine.Crx(ll),num2str(ll))
            end
            title('PALM $u$ [m/s]','interpreter','latex');
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
                text(Wp.turbine.Cry(ll)+50,Wp.turbine.Crx(ll),num2str(ll))
            end
            title('error [m/s]','interpreter','latex');
            hold off;   
            
            subplot(2,3,[4 5]);
            plot(Wp.sim.time(1:sol.k),sum(turbData.power(1:sol.k,:),2)','b');hold on
            plot(Wp.sim.time(1:sol.k),sum(Power(:,1:sol.k)),'r--');
            ylabel('$\Sigma P_i$ [W]','interpreter','latex');
            axis([0,Wp.sim.time(sol.k) 0 max(max(sum(Power(:,1:end))))+10^5])
            title('WF power PALM (blue) WFSim (red dashed)','interpreter','latex')
            grid;hold off;         
            
            subplot(2,3,6)
            plot(Wp.sim.time(1:sol.k),sum(turbData.power(1:sol.k,:),2)'-sum(Power(:,1:sol.k)));hold on
            title('error [W]','interpreter','latex');
            axis([0,Wp.sim.time(sol.k) min(min(sum(turbData.power(1:sol.k,:),2)'-sum(Power(:,1:sol.k))))-.2 max(max(sum(turbData.power(1:sol.k,:),2)'-sum(Power(:,1:sol.k))))+.2]);
            grid;hold off;
        end
        drawnow
        
    end;

end;


%% Post analysis

Nt             = Wp.sim.NN; % until what time you want to plot
turbData.power = turbData.power';
for kk=1:Wp.turbine.N
   yline(kk,:) = Wp.mesh.yline{kk}; 
end

% Plotting    
figure(2);clf
plot(Wp.sim.time(1:Wp.sim.NN),RMSE);hold on;
plot(Wp.sim.time(1:Wp.sim.NN),maxe,'r');grid;
ylabel('RMSE and max');
title(['{\color{blue}{RMSE}}, {\color{red}{max}} and meanRMSE = ',num2str(mean(RMSE),3)])

% wind farm power
figure(3);clf
plot(Wp.sim.time(1:Nt),sum(Power(:,1:Nt)),'k','Linewidth',1);hold on;
plot(Wp.sim.time(1:Nt),sum(turbData.power(:,1:Nt)),'b--');
grid;xlabel('$t [s]$','interpreter','latex');
ylabel('$P$ [W]','interpreter','latex');
title('Wind farm power: WFSim (black) PALM (blue dashed)','interpreter','latex');
xlim([0 Wp.sim.time(Nt)])

% wind turbine powers
n = size(unique(Wp.mesh.xline),1);
m = size(unique(yline,'rows'),1);

figure(4);clf
for kk=1:Wp.turbine.N
    subplot(n,m,kk)
    plot(Wp.sim.time(1:Nt),turbData.power(kk,1:Nt),'b--');hold on
    plot(Wp.sim.time(1:Nt),Power(kk,1:Nt),'k');
    axis([0,Wp.sim.time(Nt) 0 max(max(turbData.power(:,1:end)))+10^5])
    ylabel(['$P$' num2str(kk) ' [W]'], 'interpreter','latex')
    xlabel('$t [s]$','interpreter','latex');
    title('WFSim (black) PALM (blue dashed)', 'interpreter','latex')
    grid;hold off;
end

% control signals CT'
figure(5);clf
for kk=1:Wp.turbine.N
    subplot(n,m,kk)
    plot(Wp.sim.time(1:Nt),turbData.CT(1:Nt,kk));
    ylabel(['$CT^{\prime}$' num2str(kk)],'interpreter','latex');
    xlabel('$t [s]$','interpreter','latex');
    title('Thrust coefficients','interpreter','latex');
    axis([0,Wp.sim.time(Nt) 0 max(max(turbData.CT(1:Nt,:)))+.2]);
    grid;hold off;
end


% Wake centrelines
yline    = unique(yline,'rows');

indices  = [300 400 700 800]/Wp.sim.h; % sample times at which the centerline will be plotted 
                                       % select always four     

for kk=1:length(indices)
    for ll=1:m
        up(:,kk,ll)     = mean(sol_array{indices(kk)}.u(:,yline(ll,:)),2);
        temp            = squeeze(u(indices(kk),:,:));
        upPALM(:,kk,ll) = mean(temp(:,yline(ll,:)),2);
    end
end


for ll=1:m 
    figure(5+ll);clf
    for kk=1:length(indices)
        subplot(2,2,kk)
        plot(Wp.mesh.ldxx2(:,1)',up(:,kk,ll),'k','Linewidth',1);hold on;
        plot(Wp.mesh.ldxx2(:,1)',upPALM(:,kk,ll),'b--','Linewidth',1);grid;
        xlabel('$x$ [m]','interpreter','latex');ylabel(['$U^c$ [m/s] at $y$ = ' num2str(mean(Wp.mesh.ldyy2(1,yline(ll,:))),3) ' [m]'],'interpreter','latex');
        title(['t = ' num2str(indices(kk)*Wp.sim.h) ' [s]'],'interpreter','latex')
        ylim([2 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    end
end

