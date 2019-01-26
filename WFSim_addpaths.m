%% Add necessary paths for WFSim
if ispc; slashSymbol = '\'; else; slashSymbol = '/'; end;
[WFSimFolder, ~, ~] = fileparts(which([mfilename '.m']));         % Get WFSim directory
addpath([WFSimFolder slashSymbol 'bin' slashSymbol 'core']);      % Add core files
addpath([WFSimFolder slashSymbol 'bin' slashSymbol 'analysis']);  % Add analysis files
addpath(genpath([WFSimFolder slashSymbol 'libraries']));          % Add external libraries
addpath([WFSimFolder slashSymbol 'bin' slashSymbol 'mpc']);       % Add mpc files

clear WFSimFolder slashSymbol