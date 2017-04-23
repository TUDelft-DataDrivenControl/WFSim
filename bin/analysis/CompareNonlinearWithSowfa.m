clear; clc; close all;

% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 1;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name       = 'APC_3x3turb_noyaw_9turb_100x50_lin';    % Meshing name (see "\bin\core\meshing.m")
%Wp.name       = 'NoPrecursor_2turb_60x30_lin';         % Meshing name (see "\bin\core\meshing.m")
%Wp.name       = 'YawCase3_50x50_lin';                  % Meshing name (see "\bin\core\meshing.m")
%Wp.name       = 'Yawcase1_2turb_100x50_lin';
Wp.Turbulencemodel  = 'WFSim3';


Animate       = 0;                      % Show 2D flow fields every x iterations (0: no plots)
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

sourcepath = ['..\..\Data_SOWFA\' char(strtok(char(Wp.name), '_')) '\' ...
    num2str(Wp.turbine.N) 'turb_' num2str(Wp.mesh.Nx) 'x' num2str(Wp.mesh.Ny) '_'...
    char(Wp.mesh.type) '\'];
list       = dir(fullfile(sourcepath, '*.mat'));
list       = {list.name};
list       = natsort(list);

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
    
    SOWFAdata             = load([sourcepath num2str(list{k})]);
    eu                    = vec(sol.u-SOWFAdata.uq); eu(isnan(eu)) = [];
    ev                    = vec(sol.v-SOWFAdata.vq); ev(isnan(ev)) = [];
    RMSE(k)               = rms([eu;ev]);
    [maxe(k),maxeloc(k)]  = max(abs(eu));
    
    if isfield(SOWFAdata,'power')
        Powersowfa(:,k)       = SOWFAdata.power;
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
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',SOWFAdata.uq,'Linecolor','none');  colormap(hot);
            caxis([min(min(SOWFAdata.uq)) max(max(SOWFAdata.uq))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            title('SOWFA u [m/s]');
            hold off;
            
            subplot(2,3,3);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.u-SOWFAdata.uq,'Linecolor','none');  colormap(hot);
            caxis([min(min(sol.u-SOWFAdata.uq)) max(max(sol.u-SOWFAdata.uq))]);  hold all; colorbar;
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
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.v,'Linecolor','none');  colormap(hot);
            caxis([min(min(sol.v)) max(max(sol.v))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            xlabel('y [m]')
            ylabel('x [m]');
            title('WFSim v [m/s]');
            hold off;
            
            subplot(2,3,5);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',SOWFAdata.vq,'Linecolor','none');  colormap(hot);
            caxis([min(min(SOWFAdata.vq)) max(max(SOWFAdata.vq))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            xlabel('y [m]')
            title('SOWFA v [m/s]');
            hold off;
            
            subplot(2,3,6);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.v-SOWFAdata.vq,'Linecolor','none');  colormap(hot);
            caxis([min(min(sol.v-SOWFAdata.vq)) max(max(sol.v-SOWFAdata.vq))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            xlabel('y [m]')
            title('error [m/s]');
            hold off;
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
indices  = [250 500 750 999];

for k=indices
    up(:,k)       = mean(uk(:,D_ind,k),2);
    SOWFAdata     = load([sourcepath num2str(list{k})]);
    upsowfa(:,k)  = mean(SOWFAdata.uq(:,D_ind),2);
    VAF(:,k)      = vaf(upsowfa(:,k),up(:,k));
end


figure(3);clf;
subplot(2,2,1)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(1)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(1)),'b','Linewidth',1);grid;
ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(7)); end
title( ['VAF = ',num2str(VAF(indices(1)),3), '\% at $k$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
subplot(2,2,2)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(2)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(2)),'b','Linewidth',1);grid;
ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(7)); end
title( ['VAF = ',num2str(VAF(indices(2)),3), '\% at $k$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
subplot(2,2,3)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(3)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(3)),'b','Linewidth',1);grid;
xlabel('$x$ [m]','interpreter','latex');ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(7)); end
title( ['VAF = ',num2str(VAF(indices(3)),3), '\% at $k$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
subplot(2,2,4)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(4)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(4)),'b','Linewidth',1);grid;
xlabel('$x$ [m]','interpreter','latex');
ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));
if Wp.turbine.N==9; vline(Wp.turbine.Crx(7)); end
title( ['VAF = ',num2str(VAF(indices(4)),3), '\% at $k$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
if Wp.turbine.N==9
    text( -2560, 26, 'First row: WFSim (black) and SOWFA (blue)','interpreter','latex') ;
    %suptitle('First row: WFSim (black) and SOWFA (blue)')
end

%% Plot cross-sections
if Wp.turbine.N==2
    indices  = [250 500 750 999];
    
    for k=indices
        SOWFAdata     = load([sourcepath num2str(list{k})]);
        usowfa(:,:,k) = SOWFAdata.uq;
        VAF(k)        = vaf(usowfa(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,k),...
            uk(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,k));
    end
    
    figure(6);clf;
    subplot(2,2,1)
    plot(Wp.mesh.ldyy2(1,:)',uk(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(1)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldyy2(1,:)',usowfa(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(1)),'b','Linewidth',1);grid;
    vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(1)));vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(end)))
    ylabel('$u$ [m/s]','interpreter','latex');
    %ylim([3 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    title( ['VAF = ',num2str(VAF(indices(1)),3), '\% at $k$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,2)
    plot(Wp.mesh.ldyy2(1,:)',uk(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(2)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldyy2(1,:)',usowfa(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(2)),'b','Linewidth',1);grid;
    vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(1)));vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(end)))
    %ylim([3 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    title( ['VAF = ',num2str(VAF(indices(2)),3), '\% at $k$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,3)
    plot(Wp.mesh.ldyy2(1,:)',uk(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(3)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldyy2(1,:)',usowfa(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(3)),'b','Linewidth',1);grid;
    vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(1)));vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(end)))
    xlabel('$y$ [m]','interpreter','latex');ylabel('$u$ [m/s]','interpreter','latex');
    %ylim([3 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    title( ['VAF = ',num2str(VAF(indices(3)),3), '\% at $k$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,4)
    plot(Wp.mesh.ldyy2(1,:)',uk(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(4)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldyy2(1,:)',usowfa(Wp.mesh.xline(1)+round(5*length(Wp.mesh.yline{1})),:,indices(4)),'b','Linewidth',1);grid;
    vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(1)));vline(Wp.mesh.ldyy2(1,Wp.mesh.yline{1}(end)))
    xlabel('$y$ [m]','interpreter','latex');
    %ylim([3 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    title( ['VAF = ',num2str(VAF(indices(4)),3), '\% at $k$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
end

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

if Wp.turbine.N==9
    figure(9);clf;
    subplot(1,3,1)
    plot(Wp.sim.time(1:end-1),a(3,:),'b');hold on;
    plot(Wp.sim.time(1:end-1),a(6,:),'k');
    plot(Wp.sim.time(1:end-1),a(1,:),'r');grid;
    ylabel('$a$','interpreter','latex');
    xlabel('$k$ [s]','interpreter','latex');
    title('$a^1$ (blue), $a^2$ (black), $a^3$ (red)','interpreter','latex')
    subplot(1,3,2)
    plot(Wp.sim.time(1:end-1),a(2,:),'b');hold on;
    plot(Wp.sim.time(1:end-1),a(4,:),'k');
    plot(Wp.sim.time(1:end-1),a(5,:),'r');grid;
    xlabel('$k$ [s]','interpreter','latex');
    title('$a^4$ (blue), $a^5$ (black), $a^6$ (red)','interpreter','latex')
    subplot(1,3,3)
    plot(Wp.sim.time(1:end-1),a(7,:),'b');hold on;
    plot(Wp.sim.time(1:end-1),a(8,:),'k');
    plot(Wp.sim.time(1:end-1),a(9,:),'r');grid;
    title('$a^7$ (blue), $a^8$ (black), $a^9$ (red)','interpreter','latex')
    xlabel('$k$ [s]','interpreter','latex');
end

%% Mean wake centrelines
if Wp.turbine.N==9
    
    D_ind    = Wp.mesh.yline{4};
    
    clear up upsowfa VAF
    for k=indices
        up(:,k)       = mean(uk(:,D_ind,k),2);
        SOWFAdata     = load([sourcepath num2str(list{k})]);
        upsowfa(:,k)  = mean(SOWFAdata.uq(:,D_ind),2);
        VAF(:,k)      = vaf(upsowfa(:,k),up(:,k));
    end
    
    figure(4);clf;
    subplot(2,2,1)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(1)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(1)),'b','Linewidth',1);grid;
    ylabel('$U^c$ [m/s]','interpreter','latex');
    ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));
    title( ['VAF = ',num2str(VAF(indices(1)),3), '\% at $k$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,2)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(2)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(2)),'b','Linewidth',1);grid;
    ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));
    title( ['VAF = ',num2str(VAF(indices(2)),3), '\% at $k$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,3)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(3)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(3)),'b','Linewidth',1);grid;
    xlabel('$x$ [m]','interpreter','latex');ylabel('$U^c$ [m/s]','interpreter','latex');
    ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));
    title( ['VAF = ',num2str(VAF(indices(3)),3), '\% at $k$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,4)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(4)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(4)),'b','Linewidth',1);grid;
    xlabel('$x$ [m]','interpreter','latex');
    ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));
    title( ['VAF = ',num2str(VAF(indices(4)),3), '\% at $k$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
    text( -2500, 26, 'Second row: WFSim (black) and SOWFA (blue)','interpreter','latex') ;
    %suptitle('Second row: WFSim (black) and SOWFA (blue)')
    
    D_ind    = Wp.mesh.yline{3};
    
    clear up upsowfa VAF
    for k=indices
        up(:,k)       = mean(uk(:,D_ind,k),2);
        SOWFAdata     = load([sourcepath num2str(list{k})]);
        upsowfa(:,k)  = mean(SOWFAdata.uq(:,D_ind),2);
        temp          = up(:,k);
        tempsowfa     = upsowfa(:,k);
        [row, col]    = find(isnan(tempsowfa));
        temp(row)     = [];
        tempsowfa(row)= [];
        VAF(:,k)      = vaf(tempsowfa,temp);
    end
    figure(5);clf;
    subplot(2,2,1)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(1)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(1)),'b','Linewidth',1);grid;
    ylabel('$U^c$ [m/s]','interpreter','latex');
    ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));
    title( ['VAF = ',num2str(VAF(indices(1)),3), '\% at $k$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,2)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(2)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(2)),'b','Linewidth',1);grid;
    ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));
    title( ['VAF = ',num2str(VAF(indices(2)),3), '\% at $k$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,3)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(3)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(3)),'b','Linewidth',1);grid;
    xlabel('$x$ [m]','interpreter','latex');ylabel('$U^c$ [m/s]','interpreter','latex');
    ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));
    title( ['VAF = ',num2str(VAF(indices(3)),3), '\% at $k$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
    subplot(2,2,4)
    plot(Wp.mesh.ldxx2(:,1)',up(:,indices(4)),'k','Linewidth',1);hold on;
    plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(4)),'b','Linewidth',1);grid;
    xlabel('$x$ [m]','interpreter','latex');
    ylim([5 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
    vline(Wp.turbine.Crx(1));vline(Wp.turbine.Crx(2));vline(Wp.turbine.Crx(7));
    title( ['VAF = ',num2str(VAF(indices(4)),3), '\% at $k$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
    %suptitle('Third row: WFSim (black) and SOWFA (blue)')
    text( -2500, 26, 'Third row: WFSim (black) and SOWFA (blue)','interpreter','latex') ;
end

%% Power
if isfield(SOWFAdata,'power')
    if Wp.turbine.N<9
        figure(10);clf;
        subplot(2,1,1)
        plot(Power(1,1:end));hold on;
        plot(Powersowfa(1,1:end),'r--');
        grid;ylabel('Power');legend('WFSim','SOWFA')
        subplot(2,1,2)
        plot(Power(2,1:end));hold on;
        plot(Powersowfa(2,1:end),'r--');
        grid;xlabel('Time [s]');ylabel('Power')
    else
        %Power       = Power/max(max(Powersowfa));
        %Powersowfa  = Powersowfa/max(max(Powersowfa));
        clear VAF
        for k=1:9
        e(k,:)   = Power(k,:)-Powersowfa(k,:);
        VAF(:,k) = vaf(Power(k,:),Powersowfa(k,:));
        RMSE_power(:,k) = rms(e(k,:));
        end       
        
        figure(10);clf;
        subplot(3,3,1)
        plot(Power(3,1:end));hold on;
        plot(Powersowfa(3,1:end),'r');
        set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        grid;ylabel('$P_1$','interpreter','latex');ylim([0 4*10^6])
        subplot(3,3,2)
        plot(Power(6,1:end));hold on;
        plot(Powersowfa(6,1:end),'r');
        set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        grid;ylim([0 4*10^6]);ylabel('$P_2$','interpreter','latex')
        subplot(3,3,3)
        plot(Power(1,1:end));hold on;
        plot(Powersowfa(1,1:end),'r');
        set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        grid;ylim([0 4*10^6]);ylabel('$P_3$','interpreter','latex')
        subplot(3,3,4)
        plot(Power(2,1:end));hold on;
        plot(Powersowfa(2,1:end),'r');
        set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        grid;ylabel('Power');ylim([0 4*10^6]);ylabel('$P_4$','interpreter','latex')
        subplot(3,3,5)
        plot(Power(4,1:end));hold on;
        plot(Powersowfa(4,1:end),'r');
        grid;ylim([0 4*10^6]);ylabel('$P_5$','interpreter','latex')
        set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        subplot(3,3,6)
        plot(Power(5,1:end));hold on;
        plot(Powersowfa(5,1:end),'r');
        grid;ylim([0 4*10^6]);ylabel('$P_6$','interpreter','latex')
        set(gca, 'XTickLabelMode', 'manual', 'XTickLabel', []);
        set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        subplot(3,3,7)
        plot(Power(7,1:end));hold on;
        plot(Powersowfa(7,1:end),'r');
        grid;xlabel('$k [s]$','interpreter','latex');ylabel('Power');ylim([0 4*10^6]);ylabel('$P_7$','interpreter','latex')
        subplot(3,3,8)
        plot(Power(8,1:end));hold on;
        plot(Powersowfa(8,1:end),'r');
        grid;xlabel('$k [s]$','interpreter','latex');ylim([0 4*10^6]);ylabel('$P_8$','interpreter','latex')
        set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        subplot(3,3,9)
        plot(Power(9,1:end));hold on;
        plot(Powersowfa(9,1:end),'r');
        grid;xlabel('$k [s]$','interpreter','latex');ylim([0 4*10^6]);ylabel('$P_9$','interpreter','latex')
        set(gca, 'YTickLabelMode', 'manual', 'YTickLabel', []);
        %suptitle('SOWFA (red) and WFSim (blue)')
    end
    
end



