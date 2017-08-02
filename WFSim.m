clear; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%          WIND FARM SIMULATOR (WFSIM) by S. Boersma
%               Delft University of Technology, 2017
%
%            URL: https://github.com/Bartdoekemeijer/WFSim
%
%  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  %  % 
%
%%   Quick use:
%     1. Specify the wind farm you would like to simulate on line 31.
%     2. Set up the scriptOptions settings in lines 33-60. Leave the
%        default if you are unfamiliar with the settings.
%     3.
%
% Types of systems when Projection=0:
% Linearversion = 0;
% sys.A * x_{k+1}  = sys.b(x_k,u_k)
% sys.A * x_{k+1}  = sys.M * x_k + sys.m(x_k,u_k)
% Linearversion = 1;
% sys.A * dx_{k+1} = sys.Al * dx_k + sys.bl(du_k)
% sys.A * dx_{k+1} = sys.Al * dx_k + sys.Bl * du_k
% sys.Ac * \dot{x} = sys.Alc * dx + sys.Bl * u



%% Define script settings
% Wind farm settings
Wp.name      = 'YawCase3_50x50_lin_OBS';    % See bin/core/meshing.m

% Model settings
scriptOptions.Projection        = 0;        % Use projection (true/false)
scriptOptions.Linearversion     = 0;        % Provide linear variant of WFSim (true/false)
scriptOptions.exportLinearSol   = 0;        % Calculate linear solution of WFSim
scriptOptions.Derivatives       = 0;        % Compute derivatives
scriptOptions.startUniform      = 1;        % Start from a uniform flowfield (true) or a steady-state solution (false)
scriptOptions.exportPressures   = ~scriptOptions.Projection;   % Calculate pressure fields

% Convergence settings
scriptOptions.conv_eps          = 1e-6;     % Convergence threshold
scriptOptions.max_it_dyn        = 1;        % Maximum number of iterations for k > 1
if scriptOptions.startUniform==1
    scriptOptions.max_it = 1; 
else
    scriptOptions.max_it = 50;
end

% Display and visualization settings
scriptOptions.printProgress     = 1;  % Print progress every timestep
scriptOptions.printConvergence  = 1;  % Print convergence parameters every timestep
scriptOptions.Animate           = 1;  % Show 2D flow fields every x iterations (0: no plots)
scriptOptions.plotMesh          = 0;  % Show meshing and turbine locations

%%%------------------------------------------------------------------------%%%%

%% Script core
WFSim_addpaths;                             % Add paths
[Wp,sol,sys] = InitWFSim(Wp,scriptOptions); % Initialize WFSim model

% Initialize variables and figure specific to this script
sol_array = {};
CPUTime   = zeros(Wp.sim.NN,1);
if scriptOptions.Animate > 0
    %scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
           [0 0 1 1],'ToolBar','none','visible', 'on');
end


% Performing timestepping until end
disp(['Performing ' num2str(Wp.sim.NN) ' forward simulations..']);
while sol.k < Wp.sim.NN
    tic;         % Intialize timer
    [sol,sys]      = WFSim_timestepping(sol,sys,Wp,scriptOptions); % forward timestep with WFSim
    CPUTime(sol.k) = toc; % Take time
    
    % Save sol to a cell array
    sol_array{sol.k} = sol;
    
    % Display progress and animations
    if scriptOptions.printProgress
        disp(['Simulated t(' num2str(sol.k) ') = ' num2str(sol.time) ' s. CPU: ' num2str(CPUTime(sol.k)*1e3,3) ' ms.']);
    end;
    if scriptOptions.Animate > 0
        if ~rem(sol.k,scriptOptions.Animate)
            hfig = WFSim_animation(Wp,sol,hfig); 
        end; 
    end; 
end;
disp(['Completed ' num2str(Wp.sim.NN) ' forward simulations. Average CPU time: ' num2str(mean(CPUTime)*10^3,3) ' ms.']);