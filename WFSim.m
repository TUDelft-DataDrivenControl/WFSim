%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Types of systems when Projection=0:
% Linearversion = 0;
% sys.A * x_{k+1}  = sys.b(x_k,u_k)
% sys.A * x_{k+1}  = sys.M * x_k + sys.m(x_k,u_k)
% Linearversion = 1;
% sys.A * dx_{k+1} = sys.Al * dx_k + sys.bl(du_k)
% sys.A * dx_{k+1} = sys.Al * dx_k + sys.Bl * du_k
% sys.Ac * \dot{x} = sys.Alc * dx + sys.Bl * u

clear; clc; %close all;

%% Define script settings
scriptOptions.Projection        = 0;        % Use projection (true/false)
scriptOptions.Linearversion     = 0;        % Provide linear variant of WFSim (true/false)
scriptOptions.exportLinearSol   = 0;        % Calculate linear solution of WFSim
scriptOptions.Derivatives       = 0;        % Compute derivatives
scriptOptions.startUniform      = 1;        % Start from a uniform flowfield (true) or a steady-state solution (false)
scriptOptions.Animate           = 10;       % Show 2D flow fields every x iterations (0: no plots)
scriptOptions.plotMesh          = 0;        % Show meshing and turbine locations
scriptOptions.conv_eps          = 1e-6;     % Convergence threshold
scriptOptions.max_it_dyn        = 1;        % Maximum number of iterations for k > 1
scriptOptions.exportPressures   = ~scriptOptions.Projection;   % Calculate pressure fields

if scriptOptions.startUniform==1
    scriptOptions.max_it = 1; 
else
    scriptOptions.max_it = 50;
end

% WFSim general initialization script
Wp.name      = 'YawCase3_50x50_lin_OBS';
[Wp,sol,sys] = InitWFSim(Wp,scriptOptions);

% Initialize variables and figure specific to this script
uk = Wp.site.u_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
vk = Wp.site.v_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
pk = Wp.site.p_init*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);

if scriptOptions.Animate > 0
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
           [0 0 1 1],'ToolBar','none','visible', 'on');
end


%% Loop
CPUTime = zeros(Wp.sim.NN,1);
for k=1:Wp.sim.NN
    tic;
    [sol,sys] = WFSim_timestepping(sol,sys,Wp,scriptOptions);
    CPUTime(k) = toc;

    if scriptOptions.Animate > 0
        if ~rem(k,scriptOptions.Animate)
            hfig = WFSim_animation(Wp,sol,hfig); 
        end; 
    end; 
end;
disp('Completed simulations.');