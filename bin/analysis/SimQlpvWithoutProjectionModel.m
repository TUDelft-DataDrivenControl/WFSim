%% Simulate the qlpv model obtained without projection

clear; clc; close all;

%% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

%Wp.name       = 'SingleTurbine_50x50_lin';   % Meshing name (see "\bin\core\meshing.m")
Wp.name       = 'NoPrecursor_2turb_60x30_lin';   % Meshing name (see "\bin\core\meshing.m")

Wp.Turbulencemodel  = 'WFSim3';

Animate       = 5;                      % Show 2D flow fields every x iterations (0: no plots)
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
uk    = Wp.site.u_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
vk    = Wp.site.v_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
pk    = Wp.site.p_init*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);

if Animate > 0
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
end


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
        
        
        if (k==1 && it==1)
            sol.x = [Wp.site.u_Inf*ones(Wp.Nu,1);Wp.site.v_Inf*ones(Wp.Nv,1);Wp.site.p_init*ones(Wp.Np,1)];
        end
        F                 = sys.M(sys.pRCM,sys.pRCM)*sol.x(sys.pRCM,1)+sys.m(sys.pRCM);
        sol.x(sys.pRCM,1) = sys.A(sys.pRCM,sys.pRCM)\F;
        
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);         % Map solution to field
        
    end
    toc
    
    if Animate > 0
        if ~rem(k,Animate)
            Animation;
        end;
    end;
end;
disp('Completed simulations.');