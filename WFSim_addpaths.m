%% Add necessary paths for WFSim
[WFSimFolder, ~, ~] = fileparts(which([mfilename '.m']));   % Get WFSim directory
addpath([WFSimFolder '\bin\core']);                         % Add core files
addpath([WFSimFolder '\bin\analysis']);                     % Add analysis files
addpath(genpath([WFSimFolder '\libraries']));               % Add external libraries

clear WFSimFolder