function [modelOptions] = solverSet_minimal(Wp)
% This is the solverSet used for the EnKF and the UKF in WFObs

% Model settings (recommended: leave default)
modelOptions.Projection        = 0;        % Solve WFSim by projecting away the continuity equation (bool). Default: false.
modelOptions.Linearversion     = 0;        % Calculate linear system matrices of WFSim (bool).              Default: false.
modelOptions.exportLinearSol   = 0;        % Calculate linear solution of WFSim (bool).                     Default: false.
modelOptions.exportPressures   = 0;        % Calculate pressure fields. Default: 0

% Convergence settings (recommended: leave default)
modelOptions.printConvergence = 0;    % Print convergence values every timestep. Default: false.
modelOptions.conv_eps         = 1e-6; % Convergence threshold. Default: 1e-6.
modelOptions.max_it_dyn       = 1;    % Maximum number of iterations for k > 1. Default: 1.

if Wp.sim.startUniform==1
    modelOptions.max_it = 1;               % Maximum n.o. of iterations for k == 1, when startUniform = 1.
else
    modelOptions.max_it = 50;              % Maximum n.o. of iterations for k == 1, when startUniform = 0.
end
end

