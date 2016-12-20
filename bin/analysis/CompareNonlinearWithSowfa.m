clear; clc; close all;

%% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 1;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name       = 'APC_3x3turb_noyaw_9turb_100x50_lin';    % Meshing name (see "\bin\core\meshing.m")
%Wp.name       = 'NoPrecursor_2turb_60x30_lin';         % Meshing name (see "\bin\core\meshing.m")
%Wp.name       = 'YawCase3_50x50_lin';                  % Meshing name (see "\bin\core\meshing.m")

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

B2  = 2*B2;

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
        
    end
    toc
    
    SOWFAdata             = load([sourcepath num2str(list{k})]);
    eu                    = vec(sol.u-SOWFAdata.uq); eu(isnan(eu)) = [];
    ev                    = vec(sol.v-SOWFAdata.vq); ev(isnan(ev)) = [];
    RMSE(k)               = rms([eu;ev]);
    [maxe(k),maxeloc(k)]  = max(abs(eu));
    
    if Wp.turbine.N~=9
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
disp('Completed simulations.');

figure(2);clf
plot(Wp.sim.time(1:end-1),RMSE);hold on;
plot(Wp.sim.time(1:end-1),maxe,'r');grid;
ylabel('RMSE and max');
title('\color{blue}RMSE, \color{red}max')

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
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(1)),'r--','Linewidth',2);grid;
ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
title( ['VAF = ',num2str(VAF(indices(1)),3), '\% at $t$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
subplot(2,2,2)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(2)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(2)),'r--','Linewidth',2);grid;
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
title( ['VAF = ',num2str(VAF(indices(2)),3), '\% at $t$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
subplot(2,2,3)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(3)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(3)),'r--','Linewidth',2);grid;
xlabel('x [m]','interpreter','latex');ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
title( ['VAF = ',num2str(VAF(indices(3)),3), '\% at $t$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
subplot(2,2,4)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(4)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(4)),'r--','Linewidth',2);grid;
xlabel('x [m]','interpreter','latex');
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
title( ['VAF = ',num2str(VAF(indices(4)),3), '\% at $t$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
if Wp.turbine.N==9
suptitle('First column')
end
mean(RMSE)

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
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(1)),'r--','Linewidth',2);grid;
ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
title( ['VAF = ',num2str(VAF(indices(1)),3), '\% at $t$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
subplot(2,2,2)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(2)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(2)),'r--','Linewidth',2);grid;
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
title( ['VAF = ',num2str(VAF(indices(2)),3), '\% at $t$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
subplot(2,2,3)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(3)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(3)),'r--','Linewidth',2);grid;
xlabel('x [m]','interpreter','latex');ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
title( ['VAF = ',num2str(VAF(indices(3)),3), '\% at $t$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
subplot(2,2,4)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(4)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(4)),'r--','Linewidth',2);grid;
xlabel('x [m]','interpreter','latex');
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]); 
title( ['VAF = ',num2str(VAF(indices(4)),3), '\% at $t$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
suptitle('Second column')

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
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(1)),'r--','Linewidth',2);grid;
ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
title( ['VAF = ',num2str(VAF(indices(1)),3), '\% at $t$ = ', num2str(indices(1)), ' [s]'] , 'interpreter','latex')
subplot(2,2,2)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(2)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(2)),'r--','Linewidth',2);grid;
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
title( ['VAF = ',num2str(VAF(indices(2)),3), '\% at $t$ = ', num2str(indices(2)), ' [s]'] , 'interpreter','latex')
subplot(2,2,3)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(3)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(3)),'r--','Linewidth',2);grid;
xlabel('x [m]','interpreter','latex');ylabel('$U^c$ [m/s]','interpreter','latex');
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]);
title( ['VAF = ',num2str(VAF(indices(3)),3), '\% at $t$ = ', num2str(indices(3)), ' [s]'] , 'interpreter','latex')
subplot(2,2,4)
plot(Wp.mesh.ldxx2(:,1)',up(:,indices(4)),'k','Linewidth',1);hold on;
plot(Wp.mesh.ldxx2(:,1)',upsowfa(:,indices(4)),'r--','Linewidth',2);grid;
xlabel('x [m]','interpreter','latex');
ylim([0 Wp.site.u_Inf+1]);xlim([Wp.mesh.ldxx2(1,1) Wp.mesh.ldxx2(end,1)]); 
title( ['VAF = ',num2str(VAF(indices(4)),3), '\% at $t$ = ', num2str(indices(4)), ' [s]'] , 'interpreter','latex')
suptitle('Third column')

end

%% Power
if Wp.turbine.N~=9    
    figure(6);clf;
    subplot(2,1,1)
    plot(Power(1,1:end));hold on;
    plot(Powersowfa(1,1:end),'r--');
    grid;ylabel('Power');legend('WFSim','SOWFA')
    subplot(2,1,2)
    plot(Power(2,1:end));hold on;
    plot(Powersowfa(2,1:end),'r--');
    grid;xlabel('Time [s]');ylabel('Power')
end
%%
if 1
    yaw_angles = .5*Wp.turbine.Drotor*exp(1i*input{1}.phi*pi/180);  % Yaw angles
    
    figure(7);clf
    contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',SOWFAdata.uq,'Linecolor','none');  colormap(jet);
    caxis([min(min(SOWFAdata.uq)) max(max(SOWFAdata.uq))]);  hold all; colorbar;
    axis equal; axis tight;
    for ll=1:Wp.turbine.N
        Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
        Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
        plot(Qy,Qx,'k','linewidth',1)
    end
    title('Simulation set-up');
end