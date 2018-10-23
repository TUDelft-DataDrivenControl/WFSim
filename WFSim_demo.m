clear; clc; close all; %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%      WIND FARM SIMULATOR (WFSIM) by S. Boersma and B. Doekemeijer
%                 Delft University of Technology, 2018
%          Repo: https://github.com/TUDelft-DataDrivenControl/WFSim
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%
%%   Quick use:
%     1. Specify the wind farm you would like to simulate on line 73.
%        A list of all the wind farm scenarios can be found in
%        the 'layoutDefinitions' folder. You can also create your own 
%        wind farm scenario here.
%     2. Either load a predefined timeseries of control inputs using line
%        74, or alternatively add your turbine/farm controllers in lines
%        115-119.
%     3. Setup the model solver settings in line 75. The default selection
%        is 'solverSet_default(Wp)', as defined in 'solverDefintions'.
%     3. Setup the simulation settings in lines 78-81.
%     3. Press start.
%
%%   Relevant input/output variables
%     - modelOptions: this struct contains simulation settings
%     related to the wind farm itself (solution methodology, etc.)
%     - scriptOptions: this struct contains simulation settings, not
%     related to the wind farm itself (outputs, etc.)
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
%         sol.turbInput: turbine inputs at current time
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


%% Define simulation settings: layout, control inputs and simulation duration
addpath('layoutDefinitions') % Folder with predefined wind farm layouts
Wp = layoutSet_palm_4turb_adm_turbl(); % Choose which scenario to simulate. See 'layoutDefinitions' folder for the full list.
addpath('controlDefinitions') % Make use of a predefined timeseries of control inputs
turbInputSet = controlSet_palm_4turb_adm_turbl(Wp); % Choose control set 
addpath('solverDefinitions'); % Folder with model options, solver settings, etc.
modelOptions = solverSet_default(Wp); % Choose model solver options

% Simulation length, display and visualization settings
NN = floor(turbInputSet.t(end)/Wp.sim.h); % Number of timesteps in simulation
verboseOptions.printProgress = 1;    % Print progress in cmd window every timestep. Default: true.
verboseOptions.Animate       = 0;   % Plot flow fields every [X] iterations (0: no plots). Default: 10.
verboseOptions.plotMesh      = 0;    % Plot mesh, turbine locations, and print grid offset values. Default: false.


%% Script core functions
run('WFSim_addpaths.m');                    % Add essential paths to MATLABs environment
[Wp,sol,sys] = InitWFSim(Wp,modelOptions,verboseOptions.plotMesh); % Initialize WFSim model

% Initialize variables and figure specific to this script
CPUTime   = []; % Create empty matrix to save CPU timings
if verboseOptions.Animate > 0 % Create empty figure if Animation is on
    hfig = figure('color',[0 166/255 214/255],'units','normalized',...
        'outerposition',[0 0 1 1],'ToolBar','none','visible', 'on');
end

% Performing forward time propagations
disp(['Performing ' num2str(NN) ' forward simulations..']);
while sol.k < NN
    tic;                    % Start stopwatch
    
    % Determine control setting at current time by interpolation of time series
    turbInput = struct('t',sol.time);
    for i = 1:Wp.turbine.N
        
        %turbInput.CT_prime(i,1) = interp1(turbInputSet.t,turbInputSet.CT_prime(i,:),sol.time,turbInputSet.interpMethod);        
        turbInput.CT_prime(i,1)  = 2.0;
        
        if i==1
          Amsin         = 1.5;
          fsin          = 0.225*Wp.site.u_Inf/Wp.turbine.Drotor; % St*Uinf/D
          St            = fsin*Wp.turbine.Drotor/Wp.site.u_Inf;  % St_opt=0.225    
          CTprime       = (-4.0*sqrt( (1.0-0.89 )/4.0 ) + 2.0 )/( 0.5+sqrt( ( 1.0-0.89 )/4.0 ) )...
                            + Amsin * sin( 2.0*pi*St*turbInput.t*Wp.site.u_Inf/Wp.turbine.Drotor );
          
          turbInput.CT_prime(i,1) = 1.0 - 4.0*( CTprime ./( 4.0 + CTprime ) - 0.5 ).^2;
        end
        
        turbInput.phi(i,1)      = interp1(turbInputSet.t,turbInputSet.phi(i,:),     sol.time,turbInputSet.interpMethod);
    end
    
    % Propagate the WFSim model
    [sol,sys]      = WFSim_timestepping(sol,sys,Wp,turbInput,modelOptions); % forward timestep: x_k+1 = f(x_k)
    CPUTime(sol.k) = toc;   % Stop stopwatch
    
    % Save sol to cell array
    sol_array(sol.k) = sol; 
    
    % Print progress, if necessary
    if verboseOptions.printProgress
        disp(['Simulated t(' num2str(sol.k) ') = ' num2str(sol.time) ...
              ' s. CPU: ' num2str(CPUTime(sol.k)*1e3,3) ' ms.']);
    end
    
    % Plot animations, if necessary
    if verboseOptions.Animate > 0
        if ~rem(sol.k,verboseOptions.Animate)
            hfig = WFSim_animation(Wp,sol,hfig); 
        end
    end
end
disp(' ')
disp(['Completed ' num2str(NN) ' forward simulations. Average CPU time: ' num2str(mean(CPUTime)*10^3,3) ' ms.']);

%%
if Wp.turbine.N==4
    % plot apfelstrudel results
%     Powers_baseline =    1.0e+09*[
%    1.577097785855393
%    0.617718736241326
%    0.568546065420571
%    0.499913095972285];
    Powers_baseline =    1.0e+09*[
   1.659988378014872
   0.569277441672924
   0.562909939599069
   0.496278092984221];
    for kk=1:NN
        P(:,kk) =sol_array(kk).turbine.power;
    end
    Increase = (sum(sum(P,2))/sum(Powers_baseline)-1)*100;
    
    figure(2);clf
    stem([1,2,3,4],sum(P,2)./Powers_baseline);grid
    xlabel('Row','interpreter', 'latex')
    ylabel('$P/P_{baseline}$','interpreter', 'latex')
    set(gca,'xtick',0:5);xlim([0 5])
    title(['Total power increase = ',num2str(Increase,2),' \%'],'interpreter','latex')
end