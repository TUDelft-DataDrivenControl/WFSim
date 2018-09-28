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
%         sys.B1:    Submatrix in the matrix sys.A.
%         sys.B2:    Submatrix in the matrix sys.A.
%         sys.bc:    Vector with boundary conditions of the continuity equation (part of sys.b).
%
%%   Debugging and contributing:
%     - First, try to locate any errors by turning all possible outputs
%       on (printProgress, printConvergence, Animate, plotMesh).
%     - If you cannot solve your problems, reach out on the Github.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Define script settings
addpath('layoutDefinitions')
Wp = palm_6turb_adm_turbl();    % Choose which scenario to simulate. See 'layoutDefinitions' folder for the full list.
NN = 1000; % Number of timesteps

% Model settings (recommended: leave default)
scriptOptions.Projection        = 0;        % Solve WFSim by projecting away the continuity equation (bool). Default: false.
scriptOptions.Linearversion     = 1;        % Calculate linear system matrices of WFSim (bool).              Default: false.
scriptOptions.exportLinearSol   = 0;        % Calculate linear solution of WFSim (bool).                     Default: false.
scriptOptions.Derivatives       = 0;        % Compute derivatives, useful for predictive control (bool).     Default: false.
scriptOptions.exportPressures   = ~scriptOptions.Projection;   % Calculate pressure fields. Default: '~scriptOptions.Projection'

% Convergence settings (recommended: leave default)
scriptOptions.conv_eps          = 1e-6;     % Convergence threshold. Default: 1e-6.
scriptOptions.max_it_dyn        = 1;        % Maximum number of iterations for k > 1. Default: 1.

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
CPUTime   = []; % Create empty matrix to save CPU timings
if scriptOptions.Animate > 0 % Create empty figure if Animation is on
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
           [0 0 1 1],'ToolBar','none','visible', 'on');
end
if Wp.sim.startUniform==1
    scriptOptions.max_it = 1;               % Maximum n.o. of iterations for k == 1, when startUniform = 1.
else
    scriptOptions.max_it = 50;              % Maximum n.o. of iterations for k == 1, when startUniform = 0.
end

% Performing forward time propagations
disp(['Performing ' num2str(NN) ' forward simulations..']);
while sol.k < NN
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