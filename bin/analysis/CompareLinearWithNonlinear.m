%% Compare the linear model with the nonlinear model
% This script to simulates linear model and the nonlinear model
% with equivalent delta perturbations.
clear; clc; close all;

%% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 1;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 1;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name       = 'NoPrecursor_2turb_60x30_lin';   % Meshing name (see "\bin\core\meshing.m")

Wp.Turbulencemodel  = 'WFSim3';


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

%% Loop to get state-space matrices linear model 
for k=1:2
    tic
    it        = 0;
    eps       = 1e19;
    epss      = 1e20;
    
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
    
end;

%% Simulate linear and nonlinear model with same delta perturbation
Linearversion = 0;
Animate       = 10;

L       = 250;
h       = Wp.sim.h;
time    = (0:h:L);
NN      = length(time);

% System matrices
if options.Projection>0
    % ss(Atl,Btl,sys.Qsp,0,Etl,h);
    A   = sys.Etl\sys.Atl;
    nx  = size(A,1);
    B   = sys.Etl\sys.Btl;
    nw  = size(B,2);
    C   = sys.Qsp;
    Cz  = zeros(2,Wp.Nu+Wp.Nv);
    Cz(1,(Wp.mesh.xline(1)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{1}(1)-1:...
    (Wp.mesh.xline(1)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{1}(end)-1) = 1; % Extracts sum of rotor velocities first turbine 
    Cz(2,(Wp.mesh.xline(2)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{2}(1)-1:...
    (Wp.mesh.xline(2)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{2}(end)-1) = 1; % Extracts sum of rotor velocities second turbine 
    ny  = size(C,1);
else
    % ss(sys.Al,sys.Bl,I,0,sys.A,h);
    A   = sys.A\sys.Al;
    nx  = size(A,1);
    B   = sys.A\sys.Bl;
    nw  = size(B,2);
    C   = spdiags(ones(nx,1),0,nx,nx);
    Cz  = zeros(2,nx);
    Cz(1,(Wp.mesh.xline(1)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{1}(1)-1:...
    (Wp.mesh.xline(1)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{1}(end)-1) = 1/length(Wp.mesh.yline{1}); % Extracts sum of rotor velocities first turbine 
    Cz(2,(Wp.mesh.xline(2)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{2}(1)-1:...
    (Wp.mesh.xline(2)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{2}(end)-1) = 1/length(Wp.mesh.yline{2}); % Extracts sum of rotor velocities second turbine  
    ny  = size(C,1);
end

x       = zeros(nx,NN);             % State
y       = zeros(ny,NN);             % Output
zl      = zeros(Wp.turbine.N,NN);   % Normalised rotor velocities linear model
z       = zeros(Wp.turbine.N,NN);   % Normalised rotor velocities nonlinear model
w       = zeros(nw,NN);             % Input

% Define perturbations
deltabeta   = 0.1*ones(1,NN);
deltayaw    = 0*ones(1,NN);
%for kk=1:NN;if sin(2*pi*5/L*time(kk))>=0;deltabeta(kk) = .1;else deltabeta(kk) = -.1;end;end;

w(1,:)      = deltabeta;
w(2,:)      = 0*-deltabeta;
w(3,:)      = deltayaw;
w(4,:)      = -deltayaw;

for k=2:NN
    input{k}.beta = input{1}.beta + w(1:Wp.turbine.N,k);     
    input{k}.phi  = input{1}.phi  + w(Wp.turbine.N+1:end,k);
end

% Flow fields
[du,dv,ul,vl]  = deal(zeros(Wp.mesh.Nx,Wp.mesh.Ny));
uss            = sol.u;
vss            = sol.v;
Ur1            = mean(uss(Wp.mesh.xline(1),Wp.mesh.yline{1}));
Ur2            = mean(uss(Wp.mesh.xline(2),Wp.mesh.yline{2}));
z(:,1)         = [Ur1;Ur2];
% Figure
if Animate>0
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
end

% Time loop
for k=1:NN-1
       
    % Solve for linear model
    y(:,k)              = C*x(:,k);
    zl(:,k)             = Cz*y(:,k)+[Ur1;Ur2];
    x(:,k+1)            = A*x(:,k) + B*w(:,k);
    du(3:end-1,2:end-1) = reshape(y(1:Wp.Nu,k),Wp.mesh.Ny-2,Wp.mesh.Nx-3)';
    dv(2:end-1,3:end-1) = reshape(y(Wp.Nu+1:Wp.Nu+Wp.Nv,k),Wp.mesh.Ny-3,Wp.mesh.Nx-2)';
    ul                  = uss + du;
    vl                  = vss + dv;
    
    eu                   = vec(sol.u-ul); eu(isnan(eu)) = [];
    ev                   = vec(sol.v-vl); ev(isnan(ev)) = [];
    RMSE(k)              = rms([eu;ev]);
    [maxe(k),maxeloc(k)] = max(abs(eu));
    
    % Solve for nonlinear model
    [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k)] = ...
        Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options);        % Create system matrices
    [sol,sys] = Computesol(sys,input{k},sol,k,it,options);                     % Compute solution
    [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);        % Map solution to field
     
    z(:,k+1)             = Cz*sol.x;

    
    %% All plot related stuff from here
    if Animate>0
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
            text(0,Wp.mesh.ldxx2(end,end)+320,['Time ', num2str(time(k),'%.1f'), 's']);
            xlabel('y [m]')
            ylabel('x [m]');
            title('u [m/s]');
            hold off;
            
            subplot(2,3,2);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',ul,'Linecolor','none');  colormap(hot);
            caxis([min(min(ul)) max(max(ul))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            xlabel('y [m]')
            ylabel('x [m]');
            title('u_l [m/s]');
            hold off;
            
            subplot(2,3,3);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.u-ul,'Linecolor','none');  colormap(hot);
            caxis([min(min(sol.u-ul)) max(max(sol.u-ul))]);  hold all; colorbar;
            axis equal; axis tight;
            ldyyv = Wp.mesh.ldyy(:); ldxx2v = Wp.mesh.ldxx2(:);
            plot(ldyyv(maxeloc(k)),ldxx2v(maxeloc(k)),'whiteo','LineWidth',1,'MarkerSize',8,'DisplayName','Maximum error location');
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            xlabel('y [m]')
            ylabel('x [m]');
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
            title('v [m/s]');
            hold off;
            
            subplot(2,3,5);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',vl,'Linecolor','none');  colormap(hot);
            caxis([min(min(vl)) max(max(vl))]);  hold all; colorbar;
            axis equal; axis tight;
            for ll=1:Wp.turbine.N
                Qy     = (Wp.turbine.Cry(ll)-real(yaw_angles(ll))):1:(Wp.turbine.Cry(ll)+real(yaw_angles(ll)));
                Qx     = linspace(Wp.turbine.Crx(ll)-imag(yaw_angles(ll)),Wp.turbine.Crx(ll)+imag(yaw_angles(ll)),length(Qy));
                plot(Qy,Qx,'k','linewidth',1)
            end
            xlabel('y [m]')
            ylabel('x [m]');
            title('v_l [m/s]');
            hold off;
            
            subplot(2,3,6);
            contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.v-vl,'Linecolor','none');  colormap(hot);
            caxis([min(min(sol.v-vl)) max(max(sol.v-vl))]);  hold all; colorbar;
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
            drawnow;
        end
    end
        
end
%%
figure(2);clf
subplot(3,1,1)
plot(time(1:end-1),RMSE);hold on;
plot(time(1:end-1),maxe,'r');grid;
ylabel('RMSE and max');
title('\color{blue}RMSE, \color{red}max')
subplot(3,1,2)
plot(time,w(1,:));hold on;
plot(time,w(2,:),'r');
plot(time,w(3,:),'k');
plot(time,w(4,:),'m');grid;
ylabel('\delta \beta, \delta \gamma');xlabel('Time [s]');
ylim([1.1*min(-max(deltayaw),-max(deltabeta)) 1.1*max(max(deltayaw),max(deltabeta))])
title('\color{blue}\delta \beta_1, \color{red}\delta \beta_2, \color{black}\delta \gamma_1, \color{magenta}\delta \gamma_2')
subplot(3,1,3)
plot(time(1:end-1),z(1,1:end-1));hold on;
plot(time(1:end-1),z(2,1:end-1),'r');
plot(time(1:end-1),zl(1,1:end-1),'b--');
plot(time(1:end-1),zl(2,1:end-1),'r--');grid;
ylabel('U_{r_1}, U_{r_2}');xlabel('Time [s]');
title('\color{blue} U_{r_1}, \color{red} U_{r_2}');
ylim([0 8])