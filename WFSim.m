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

WFSim_addpaths

%% Initialize script
options.Projection     = 0;                      % Use projection (true/false)
options.Linearversion  = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                      % Calculate linear solution of WFSim
options.Derivatives    = 0;                      % Compute derivatives
options.startUniform   = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;    % Calculate pressure fields

Wp.name             = 'ThreeTurbine_Ampc';      % Meshing name (see "\bin\core\meshing.m")

Wp.Turbulencemodel  = 'WFSim3';


Animate       = 20;                      % Show 2D flow fields every x iterations (0: no plots)
plotMesh      = 0;                      % Show meshing and turbine locations
conv_eps      = 1e-6;                   % Convergence threshold
max_it_dyn    = 1;                      % Maximum number of iterations for k > 1

if options.startUniform==1
    max_it = 1; 
else
    max_it = 50;
end

% WFSim general initialization script
[Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] = InitWFSim(Wp,options,plotMesh);

% Initialize variables and figure specific to this script
uk = Wp.site.u_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
vk = Wp.site.v_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
pk = Wp.site.p_init*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);

if Animate > 0
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
           [0 0 1 1],'ToolBar','none','visible', 'on');
end

%% Loop
CPUTime = zeros(Wp.sim.NN,1);
for k=1:Wp.sim.NN
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
                    Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options);              % Create system matrices
        [sol,sys] = Computesol(sys,input{k},sol,k,it,options);                      % Compute solution
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);            % Map solution to field
        
        %display(['k ',num2str(k,'%-1000.1f'),', It ',num2str(it,'%-1000.0f'),', Nv=', num2str(Normv{k}(it),'%10.2e'), ', Nu=', num2str(Normu{k}(it),'%10.2e'), ', TN=',num2str(eps,'%10.2e'),', Np=','Mean effective=',num2str(mean(Ueffect(1,k)),'%-1000.2f')]) ;
    end
    CPUTime(k) = toc;

    if Animate > 0
        if ~rem(k,Animate)
            Animation; 
        end; 
    end; 
end;
disp('Completed simulations.');