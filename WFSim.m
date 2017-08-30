clear; clc; close all; %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%      WIND FARM SIMULATOR (WFSIM) by S. Boersma and B. Doekemeijer
%                 Delft University of Technology, 2017
%              Repo: https://github.com/Bartdoekemeijer/WFSim
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%
%%   Quick use:
%     1. Specify the wind farm you would like to simulate on line 68.
%        A list of all the wind farm scenarios can be found in
%        'bin/core/meshing.m'. You can also create your own wind farm
%        scenario here.
%     2. Set up the scriptOptions settings in lines 71-91. Leave the
%        default if you are unfamiliar with the code.
%     3. Press start.
%
%%   Relevant input/output variables
%     - scriptOptions: this struct contains all simulation settings, not
%     related to the wind farm itself (solution methodology, outputs, etc.)
%
%     - Wp: this struct contains all the simulation settings related to the
%           wind farm, the turbine inputs, the atmospheric properties, etc.
%         Wp.Nu:      Number of model states concerning longitudinal flow.
%         Wp.Nv:      Number of model states concerning lateral flow.
%         Wp.Np:      Number of model states concerning pressure terms.
%         Wp.sim:     Substruct containing timestep and simulation length.
%         Wp.turbine: Substruct containing turbine properties and settings.
%         Wp.site:    Substruct containing freestream atmospheric properties.
%         Wp.mesh:    Substruct containing topology and meshing settings.
%
%     - sol: this struct contains the system states at a certain timestep.
%         sol.k:     Discrete timestep  to which these system states belong
%         sol.time:  Actual time (in s) to which these system states belong
%         sol.x:     True system state (basically flow field excluding bcs)
%         sol.u:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.v:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.p:     Instantaneous pressure field over the mesh (in Pa)
%         sol.uu:    Same as sol.u, used for convergence
%         sol.vv:    Same as sol.v, used for convergence
%         sol.pp:    Same as sol.p, used for convergence
%         sol.turbine: a struct containing relevant turbine outputs such as
%         the ax. ind. factor, the generated power, and the ct coefficient
%         sol.measuredData: a struct containing the true flow field and the
%         measurements used for estimation
%         sol.score: a struct containing estimator performance measures
%         such as the magnitude and location of the maximum flow estimation
%         error, the RMSE, and the computational cost (CPU time)
%
%     - sys: this struct contains the system matrices at a certain timestep.
%         sys.A:     System matrix A in the grand picture: A*sol.x = b
%         sys.b:     System vector b in the grand picture: A*sol.x = b
%         sys.pRCM:  Reverse Cuthill-McKee algorithm for solving A*x=b faster.
%         sys.B1:    Important matrix in the boundary conditions.
%         sys.B2:    Important matrix in the boundary conditions.
%         sys.bc:    Important vector in the boundary conditions.
%
%%   Debugging and contributing:
%     - First, try to locate any errors by turning all possible outputs
%       on (printProgress, printConvergence, Animate, plotMesh).
%     - If you cannot solve your problems, reach out on the Github.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
<<<<<<< HEAD
% Types of systems when Projection=0:
% Linearversion = 0;
% sys.A * x_{k+1}  = sys.b(x_k,u_k)
% sys.A * x_{k+1}  = sys.M * x_k + sys.m(x_k,u_k)
% Linearversion = 1;
% sys.A * dx_{k+1} = sys.Al * dx_k + sys.bl(du_k)
% sys.A * dx_{k+1} = sys.Al * dx_k + sys.Bl * du_k

clear; clc; close all;

WFSim_addpaths

%% Initialize script
options.Projection     = 0;                      % Use projection (true/false)
options.Linearversion  = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                      % Calculate linear solution of WFSim
options.Derivatives    = 0;                      % Compute derivatives
options.startUniform   = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;    % Calculate pressure fields
  
%Wp.name             = 'WP_CPUTime';      % Meshing name (see "\bin\core\meshing.m")
%Wp.name             = 'WP_CPUTime';
Wp.name             = '2turb_adm';

Wp.Turbulencemodel  = 'WFSim4';

Animate       = 0;                      % Show 2D flow fields every x iterations (0: no plots)
plotMesh      = 0;                      % Show meshing and turbine locations
conv_eps      = 1e-6;                   % Convergence threshold
max_it_dyn    = 1;                      % Maximum number of iterations for k > 1
=======


%% Define script settings
Wp.name      = 'apc_9turb_alm_turb';    % Choose which scenario to simulate. See 'bin/core/meshing.m' for the full list.

% Model settings (recommended: leave default)
scriptOptions.Projection        = 0;        % Solve WFSim by projecting away the continuity equation (bool). Default: false.
scriptOptions.Linearversion     = 0;        % Calculate linear system matrices of WFSim (bool).              Default: false.
scriptOptions.exportLinearSol   = 0;        % Calculate linear solution of WFSim (bool).                     Default: false.
scriptOptions.Derivatives       = 0;        % Compute derivatives, useful for predictive control (bool).     Default: false.
scriptOptions.startUniform      = 1;        % Start with a uniform flowfield (true) or with fully developed wakes (false).
scriptOptions.exportPressures   = ~scriptOptions.Projection;   % Calculate pressure fields. Default: '~scriptOptions.Projection'
>>>>>>> master

% Convergence settings (recommended: leave default)
scriptOptions.conv_eps          = 1e-6;     % Convergence threshold. Default: 1e-6.
scriptOptions.max_it_dyn        = 1;        % Maximum number of iterations for k > 1. Default: 1.
if scriptOptions.startUniform==1
    scriptOptions.max_it = 1;               % Maximum n.o. of iterations for k == 1, when startUniform = 1.
else
    scriptOptions.max_it = 50;              % Maximum n.o. of iterations for k == 1, when startUniform = 0.
end

<<<<<<< HEAD
%profile on

%vid = VideoWriter('flow.avi');
%open(vid);

inc         = 1;
meanCPUTime = zeros(inc,2);
for ll=1:inc
   
    Wp.ll = ll;

% WFSim general initialization script
[Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] = InitWFSim(Wp,options,plotMesh);

if Animate > 0
    scrsz = get(0,'ScreenSize');
=======
% Display and visualization settings
scriptOptions.printProgress     = 1;    % Print progress in cmd window every timestep. Default: true.
scriptOptions.printConvergence  = 0;    % Print convergence values every timestep.     Default: false.
scriptOptions.Animate           = 10;   % Plot flow fields every [X] iterations (0: no plots). Default: 10.
scriptOptions.plotMesh          = 1;    % Plot mesh, turbine locations, and print grid offset values. Default: false.



%% Script core functions
run('WFSim_addpaths.m');                    % Add essential paths to MATLABs environment
[Wp,sol,sys] = InitWFSim(Wp,scriptOptions); % Initialize WFSim model

% Initialize variables and figure specific to this script
sol_array = {}; % Create empty array to save 'sol' to at each time instant
CPUTime   = zeros(Wp.sim.NN,1); % Create empty matrix to save CPU timings
if scriptOptions.Animate > 0 % Create empty figure if Animation is on
>>>>>>> master
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
           [0 0 1 1],'ToolBar','none','visible', 'on');
end

<<<<<<< HEAD
CPUTime = zeros(Wp.sim.NN,1);
%% Loop 
for k=1:Wp.sim.NN
    tic
    it        = 0;
    eps       = 1e19;
    epss      = 1e20;
    sol.uk    = sol.u; 
    sol.vk    = sol.v;
    u(:,:,k)  = sol.u;
    v(:,:,k)  = sol.v;
    
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
        
        display(['k ',num2str(k,'%-1000.1f'),', It ',num2str(it,'%-1000.0f'),]);
    end
    CPUTime(k) = toc;
    
    %if k==50
        %Wp.site.u_Inf           = 12;
        %[B1,B2,bc]              = Compute_B1_B2_bc(Wp);
        %sol.u(1:2,1:Wp.mesh.Ny) = Wp.site.u_Inf;
        %Wp.site.v_Inf           = 2;    
        %sol.v(1:2,1:Wp.mesh.Ny) = Wp.site.v_Inf;
    %end
    
    if Animate > 0
        if ~rem(k,Animate)
            Animation; 
        end; 
    end; 
    
    %fname = sprintf('D:/sboersma3/My Documents/MATLAB/WFSim/data_WFSim/Change axial induction/%d',k);
    %fname = sprintf('D:/sboersma3/My Documents/MATLAB/WFSim/data_WFSim/Change yaw/%d',k);
    %save(fname,'sol','Wp','k','input')
    
    %frame = getframe(gcf);
    %writeVideo(vid,frame);
    q(k,:) = sys.q;
end
%close(vid);

meanCPUTime(ll,:) = [mean(CPUTime(2:end)) size(sol.x,1)];
%meanCPUTimeWFSim3 = meanCPUTime; save('WFSim3extended.mat','meanCPUTimeWFSim3extended')
end
%profile viewer

disp('Completed simulations.');
=======
% Performing forward time propagations
disp(['Performing ' num2str(Wp.sim.NN) ' forward simulations..']);
while sol.k < Wp.sim.NN
    tic;                    % Start stopwatch
    [sol,sys]      = WFSim_timestepping(sol,sys,Wp,scriptOptions); % forward timestep: x_k+1 = f(x_k)
    CPUTime(sol.k) = toc;   % Stop stopwatch
    
    % Save sol to cell array
    sol_array{sol.k} = sol; 
    
    % Print progress, if necessary
    if scriptOptions.printProgress
        disp(['Simulated t(' num2str(sol.k) ') = ' num2str(sol.time) ...
              ' s. CPU: ' num2str(CPUTime(sol.k)*1e3,3) ' ms.']);
    end;
    
    % Plot animations, if necessary
    if scriptOptions.Animate > 0
        if ~rem(sol.k,scriptOptions.Animate)
            hfig = WFSim_animation(Wp,sol,hfig); 
        end; 
    end; 
end;
disp(['Completed ' num2str(Wp.sim.NN) ' forward simulations. Average CPU time: ' num2str(mean(CPUTime)*10^3,3) ' ms.']);
>>>>>>> master
