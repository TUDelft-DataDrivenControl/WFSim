clear; clc; close all; %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%             WIND FARM SIMULATOR (WFSIM) by S. Boersma
%                 Delft University of Technology, 2017
%              Repo: https://github.com/Bartdoekemeijer/WFSim
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%
%%   Quick use:
%     1. Specify the wind farm you would like to simulate on line 64.
%        A list of all the wind farm scenarios can be found in
%        'bin/core/meshing.m'. You can also create your own wind farm
%        scenario here.
%     2. Set up the scriptOptions settings in lines 67-87. Leave the
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
%         sol.u:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.v:     Instantaneous longitudinal flow field over the mesh (in m/s)
%         sol.p:     Instantaneous pressure field over the mesh (in Pa)
%         sol.uu:    Same as sol.u, used for convergence.
%         sol.vv:    Same as sol.v, used for convergence.
%         sol.pp:    Same as sol.p, used for convergence.
%         sol.a:     Axial induction factor of each turbine at time sol.k.
%         sol.power: Generated power (in W) of each turbine at time sol.k.
%         sol.ct:    Thrust coefficient (-) of each turbine at time sol.k.
%         sol.x:     True system state (basically flow field excluding bcs).
%
%     - sys: this struct contains the system matrices at a certain timestep.
%         sys.A:     System matrix A in the grand picture: A*sol.x = b
%         sys.b:     System vector b in the grand picture: A*sol.x = b
%         sys.pRCM:  Reverse Cuthill-McKee algorithm for solving A*x=b faster.
%         sys.B1:    Important matrix in the boundary conditions.
%         sys.B2:    Important matrix in the boundary conditions.
%         sys.bc:    Important vector in the boundary conditions.

%%   Debugging and contributing:
%     - First, try to locate any errors by turning all possible outputs
%       on (printProgress, printConvergence, Animate, plotMesh).
%     - If you cannot solve your problems, reach out on the Github.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Define script settings
Wp.name      = '2turb_demo';    % Choose which scenario to simulate. See 'bin/core/meshing.m' for the full list.

% Model settings (recommended: leave default)
scriptOptions.Projection        = 0;        % Solve WFSim by projecting away the continuity equation (bool). Default: false.
scriptOptions.Linearversion     = 0;        % Calculate linear system matrices of WFSim (bool).              Default: false.
scriptOptions.exportLinearSol   = 0;        % Calculate linear solution of WFSim (bool).                     Default: false.
scriptOptions.Derivatives       = 0;        % Compute derivatives, useful for predictive control (bool).     Default: false.
scriptOptions.startUniform      = 1;        % Start with a uniform flowfield (true) or with fully developed wakes (false).
scriptOptions.exportPressures   = ~scriptOptions.Projection;   % Calculate pressure fields. Default: '~scriptOptions.Projection'

% Convergence settings (recommended: leave default)
scriptOptions.conv_eps          = 1e-6;     % Convergence threshold. Default: 1e-6.
scriptOptions.max_it_dyn        = 1;        % Maximum number of iterations for k > 1. Default: 1.
if scriptOptions.startUniform==1
    scriptOptions.max_it = 1;               % Maximum n.o. of iterations for k == 1, when startUniform = 1.
else
    scriptOptions.max_it = 50;              % Maximum n.o. of iterations for k == 1, when startUniform = 0.
end

% Display and visualization settings
scriptOptions.printProgress     = 1;    % Print progress in cmd window every timestep. Default: true.
scriptOptions.printConvergence  = 1;    % Print convergence values every timestep.     Default: false.
scriptOptions.Animate           = 10;   % Plot flow fields every [X] iterations (0: no plots). Default: 10.
scriptOptions.plotMesh          = 1;    % Plot mesh, turbine locations, and print grid offset values. Default: false.



%% Script core functions
run('WFSim_addpaths.m');                    % Add essential paths to MATLABs environment
[Wp,sol,sys] = InitWFSim(Wp,scriptOptions); % Initialize WFSim model

% Initialize variables and figure specific to this script
sol_array = {}; % Create empty array to save 'sol' to at each time instant
CPUTime   = zeros(Wp.sim.NN,1); % Create empty matrix to save CPU timings
if scriptOptions.Animate > 0 % Create empty figure if Animation is on
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
           [0 0 1 1],'ToolBar','none','visible', 'on');
end

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
