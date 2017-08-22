function [ Wp, input ] = meshing( Wp, plotMesh, PrintGridMismatch )
%MESHING2 Meshing for the WFSim code
%  Usage: [ Wp, input ] = meshing( Wp, plotMesh, PrintGridMismatch )
%   This script produces linear/exponential meshes for the WFSim flow
%   solver. Inputs are:
%   Wp:         structure containing Wp.name -- a string with the selected meshing.
%   Lingrid:    Linear meshing (true/false). Recommended: false.
%   plotMesh:   show resulting meshings (true/false).
%
%   run this code without any inputs for an example.
%
%   Date: October 18th, 2016.
%   Author: Bart

% Default settings
if nargin <= 0; error('Please specify a meshing case.'); end;
if nargin <= 1; plotMesh = true;                         end;
if nargin <= 2; PrintGridMismatch = true;                end;

% Some pre-processing to determine correct input file location
meshingloc  = which('meshing.m');
% next 4 lines I comment since WFSimfolder appeard empty with Windows OS
% I added two lines after from previous version
%if ispc; slashSymbol = '\'; else; slashSymbol = '/'; end;
%strbslash   = strfind(meshingloc,slashSymbol);
%WFSimfolder = meshingloc(1:strbslash(end-2));
%WFSimfolder = WFSimfolder(length(pwd)+2:end);
strbslash   = strfind(meshingloc,'\');
WFSimfolder = meshingloc(1:strbslash(end-2));

switch lower(Wp.name)
    % Wind farms for which PALM data is available
        case lower('1turb_adm')
        type   = 'lin';                     % Meshing type ('lin' or 'exp')
        Lx     = 7700-5500;                % Domain length in x-direction (m)
        Ly     = 1660-900;                % Domain length in y-direction (m)
        Nx     = 50;                       % Number of grid points in x-direction
        Ny     = 25;                        % Number of grid points in y-direction
        Crx    = 6000-5500;        % Turbine locations in x-direction (m)
        Cry    = 1280-900;        % Turbine locations in y-direction (m)
        
        %M = [Time   UR  Uinf  Ct_adm  a Yaw Thrust Power  WFPower]       
        %M = [Time   UR  Uinf  Ct_adm  a Yaw Thrust Power  WFPower]       
        for kk=1:length(Crx)
            M{kk} = load(['../../Data_PALM/' char(Wp.name) '/' char(Wp.name) '_matlab_turbine_parameters0' num2str(kk) '.txt']);
        end
                 
        % Correctly format inputs (temporary function)
        for j = 1:length(M{1})
            input{j}.t    = M{1}(j,1);
            input{j}.beta = M{1}(j,5)...
                            ./(1-M{1}(j,5));
            input{j}.phi  = M{1}(j,6);
            input{j}.CT   = M{1}(j,4);
        end;
        
        % Calculate delta inputs 
        for j = 1:length(M{1})-1
            input{j}.dbeta = input{j+1}.beta-input{j}.beta;
            input{j}.dphi  = input{j+1}.phi-input{j}.phi ;
        end;
        
        Drotor      = 126;    % Turbine rotor diameter in (m)
        powerscale  = 1;    % Turbine powerscaling
        forcescale  = 1.8;    % Turbine force scaling
        
        h        = 1;         % Sampling time (s)
        L        = floor(M{1}(end-1,1));% Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = mean(M{1}(1,3));   % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = .5;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 4;
        
    case lower('2turb_adm')
        type   = 'lin';                     % Meshing type ('lin' or 'exp')
        Lx     = 7500-5500;                % Domain length in x-direction (m)
        Ly     = 1460-900;                % Domain length in y-direction (m)
        Nx     = 50;                       % Number of grid points in x-direction
        Ny     = 25;                        % Number of grid points in y-direction
        Crx    = [5700, 6456]-5500;        % Turbine locations in x-direction (m)
        Cry    = [1175, 1175]-900;        % Turbine locations in y-direction (m)
        
       %M = [Time   UR  Uinf  Ct_adm  a Yaw Thrust Power  WFPower]       
        for kk=1:length(Crx)
            M{kk} = load(['../../Data_PALM/' char(Wp.name) '/' char(Wp.name) '_matlab_turbine_parameters0' num2str(kk) '.txt']);
        end
                 
        % Correctly format inputs (temporary function)
        for j = 1:length(M{1})
            input{j}.t    = M{1}(j,1);
            input{j}.beta = [M{1}(j,5);M{2}(j,5)]...
                            ./(1-[M{1}(j,5);M{2}(j,5)]);
            input{j}.phi  = [M{1}(j,6);M{2}(j,6)];
            input{j}.CT   = [M{1}(j,4);M{2}(j,4)];
        end;
        
        % Calculate delta inputs 
        for j = 1:length(M{1})-1
            input{j}.dbeta = input{j+1}.beta-input{j}.beta;
            input{j}.dphi  = input{j+1}.phi-input{j}.phi ;
        end;
        
        Drotor      = 126;    % Turbine rotor diameter in (m)
        powerscale  = 1;    % Turbine powerscaling
        forcescale  = 1.9;    % Turbine force scaling
        
        h        = 1;         % Sampling time (s)
        L        = floor(M{1}(end-1,1));% Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = mean(M{1}(:,3));   % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = .4;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 4;
  case lower('2turb_adm')
        type   = 'lin';                     % Meshing type ('lin' or 'exp')
        Lx     = 7500-5500;                % Domain length in x-direction (m)
        Ly     = 1460-900;                % Domain length in y-direction (m)
        Nx     = 50;                       % Number of grid points in x-direction
        Ny     = 25;                        % Number of grid points in y-direction
        Crx    = [5700, 6456]-5500;        % Turbine locations in x-direction (m)
        Cry    = [1175, 1175]-900;        % Turbine locations in y-direction (m)
        
       %M = [Time   UR  Uinf  Ct_adm  a Yaw Thrust Power  WFPower]       
        for kk=1:length(Crx)
            M{kk} = load(['../../Data_PALM/' char(Wp.name) '/' char(Wp.name) '_matlab_turbine_parameters0' num2str(kk) '.txt']);
        end
                 
        % Correctly format inputs (temporary function)
        for j = 1:length(M{1})
            input{j}.t    = M{1}(j,1);
            input{j}.beta = [M{1}(j,5);M{2}(j,5)]...
                            ./(1-[M{1}(j,5);M{2}(j,5)]);
            input{j}.phi  = [M{1}(j,6);M{2}(j,6)];
            input{j}.CT   = [M{1}(j,4);M{2}(j,4)];
        end;
        
        % Calculate delta inputs 
        for j = 1:length(M{1})-1
            input{j}.dbeta = input{j+1}.beta-input{j}.beta;
            input{j}.dphi  = input{j+1}.phi-input{j}.phi ;
        end;
        
        Drotor      = 126;    % Turbine rotor diameter in (m)
        powerscale  = 1;    % Turbine powerscaling
        forcescale  = 1.9;    % Turbine force scaling
        
        h        = 1;         % Sampling time (s)
        L        = floor(M{1}(end-1,1));% Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = mean(M{1}(:,3));   % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = .4;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 4;
        
     case lower('2turb')
        type   = 'lin';                     % Meshing type ('lin' or 'exp')
        Lx     = 1900-100;                      % Domain length in x-direction (m)
        Ly     = 1260-200;                      % Domain length in y-direction (m)
        Nx     = 50;                       % Number of grid points in x-direction
        Ny     = 25;                        % Number of grid points in y-direction
        Cry    = [900, 900]-100;
        Crx    = [500, 1258]-200;
        
        for kk=1:length(Crx)
            M{kk} = dlmread(['../../Data_PALM/' char(Wp.name) '/' char(Wp.name) '_matlab_turbine_parameters0' num2str(kk) '.txt'],'',1,0);
        end
                 
        % Correctly format inputs (temporary function)
        for j = 1:length(M{1})
            input{j}.t    = M{1}(j,1);
            input{j}.beta = [.3;.3];
            input{j}.phi  = [M{1}(j,11);M{2}(j,11)];
            input{j}.CT   = [8/9;8/9];
        end;
        
        % Calculate delta inputs 
        for j = 1:length(M{1})-1
            input{j}.dbeta = input{j+1}.beta-input{j}.beta;
            input{j}.dphi  = input{j+1}.phi-input{j}.phi ;
        end;
        
        Drotor      = 120;    % Turbine rotor diameter in (m)
        powerscale  = .9;    % Turbine powerscaling
        forcescale  = 1.9;    % Turbine force scaling
        
        h        = 1;         % Sampling time (s)
        L        = floor(M{1}(end-1,1));% Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = mean(M{1}(:,3));   % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = .5;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 4;
        
        case lower('9turb_adm_turb')
        type   = 'lin';                     % Meshing type ('lin' or 'exp')
        Lx     = 2800;                      % Domain length in x-direction (m)
        Ly     = 1360;                      % Domain length in y-direction (m)
        Nx     = 100;                       % Number of grid points in x-direction
        Ny     = 50;                        % Number of grid points in y-direction
        Cry    = [300.0, 300.0, 300.0, 680.0, 680.0, 680.0, 1060.0, 1060.0, 1060.0];% Turbine locations in x-direction (m)
        Crx    = [300.0, 930.0, 1560.0, 300.0, 930.0, 1560.0, 300.0, 930.0, 1560.0];% Turbine locations in y-direction (m)
        
       %M = [Time   UR  Uinf  Ct_adm  a Yaw Thrust Power  WFPower]       
        for kk=1:length(Crx)
            M{kk} = load(['../../Data_PALM/' char(Wp.name) '/' char(Wp.name) '_matlab_turbine_parameters0' num2str(kk) '.txt']);
        end
                 
        % Correctly format inputs (temporary function)
        for j = 1:length(M{1})
            input{j}.t    = M{1}(j,1);
            input{j}.beta = [M{1}(j,5);M{2}(j,5);M{3}(j,5);M{4}(j,5);M{5}(j,5);M{6}(j,5);M{7}(j,5);M{8}(j,5);M{9}(j,5)]...
                            ./(1-[M{1}(j,5);M{2}(j,5);M{3}(j,5);M{4}(j,5);M{5}(j,5);M{6}(j,5);M{7}(j,5);M{8}(j,5);M{9}(j,5)]);
            input{j}.phi  = [M{1}(j,6);M{2}(j,6);M{3}(j,6);M{4}(j,6);M{5}(j,6);M{6}(j,6);M{7}(j,6);M{8}(j,6);M{9}(j,6)];
            input{j}.CT   = [M{1}(j,4);M{2}(j,4);M{3}(j,4);M{4}(j,4);M{5}(j,4);M{6}(j,4);M{7}(j,4);M{8}(j,4);M{9}(j,4)];
        end;
        
        % Calculate delta inputs 
        for j = 1:length(M{1})-1
            input{j}.dbeta = input{j+1}.beta-input{j}.beta;
            input{j}.dphi  = input{j+1}.phi-input{j}.phi ;
        end;
        
        Drotor      = 120;    % Turbine rotor diameter in (m)
        powerscale  = .9;    % Turbine powerscaling
        forcescale  = 2.0;    % Turbine force scaling
        
        h        = 1;         % Sampling time (s)
        L        = floor(M{1}(end-1,1));% Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = mean(M{1}(:,3));   % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = .8;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 4;
        
        % Wind farms for which SOWFA data is available
    case lower('YawCase3_50x50_lin_OBS')
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = 2481.9702;      % Domain length in x-direction (m)
        Ly     = 1400;           % Domain length in y-direction (m)
        Nx     = 50;             % Number of grid points in x-direction
        Ny     = 25;             % Number of grid points in y-direction
        Crx    = [400, 1281.97]; % Turbine locations in x-direction (m)
        Cry    = [700, 700];     % Turbine locations in y-direction (m)
        
        loadedinput = load([WFSimfolder 'data_SOWFA/YawCase3/system_input.mat']); % load input settings
        Wp.Turbulencemodel  = 'WFSim3';
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = loadedinput.input.beta(j,:)';
            input{j}.phi  = loadedinput.input.phi(j,:)';
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = loadedinput.input.beta(j+1,:)'- loadedinput.input.beta(j,:)';
            input{j}.dphi  = loadedinput.input.phi(j+1,:)' - loadedinput.input.phi(j,:)' ;
        end;
        
        Drotor      = 126.4;  % Turbine rotor diameter in (m)
        powerscale  = 1.0;    % Turbine powerscaling
        forcescale  = 0.75;   % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 999;       % Simulation length (s)
        mu       = 10;        % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 8.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = 1;         % Mixing length in x-direction (m)
        turbul   = false;     % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 8;
        
    case lower('YawCase3_50x50_lin')
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = 2481.9702;      % Domain length in x-direction (m)
        Ly     = 1400;           % Domain length in y-direction (m)
        Nx     = 50;             % Number of grid points in x-direction
        Ny     = 25;             % Number of grid points in y-direction
        Crx    = [400, 1281.97]; % Turbine locations in x-direction (m)
        Cry    = [700, 700];     % Turbine locations in y-direction (m)
        
        loadedinput = load([WFSimfolder 'data_SOWFA/YawCase3/system_input.mat']); % load input settings
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = loadedinput.input.beta(j,:)';
            input{j}.phi  = loadedinput.input.phi(j,:)';
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = loadedinput.input.beta(j+1,:)'- loadedinput.input.beta(j,:)';
            input{j}.dphi  = loadedinput.input.phi(j+1,:)' - loadedinput.input.phi(j,:)' ;
        end;
        
        Drotor      = 126.4;  % Turbine rotor diameter in (m)
        powerscale  = 1.0;    % Turbine powerscaling
        forcescale  = 1.2;    % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 999;       % Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 8.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = 1;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 8;
        
    case lower('YawCase3_50x50_exp')
        type     = 'exp';      % Meshing type ('lin' or 'exp')
        Lx       = 2500;       % Domain length in x-direction (m)
        Ly       = 1400;        % Domain length in y-direction (m)
        cellSize = [50, 50];   % Approximate cell size for equidistant points (m)
        Crx      = [400, 1300];% Turbine locations in x-direction (m)
        Cry      = [700, 700]; % Turbine locations in y-direction (m)
        R_con    = [60, 60];   % Concentration radii (m) (default: equal to rotor radius)
        N_con    = [11, 11];   % Number of grid points inside concentration radii
        dx_min   = [3, 3];     % Minimal cell size inside concentration radii (m)
        
        loadedinput = load([WFSimfolder 'data_SOWFA/YawCase3/system_input.mat']); % load input settings
        loadedinput.input.phi = 0*loadedinput.input.phi;
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = loadedinput.input.beta(j,:)';
            input{j}.phi  = loadedinput.input.phi(j,:)';
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = loadedinput.input.beta(j+1,:)'- loadedinput.input.beta(j,:)';
            input{j}.dphi  = loadedinput.input.phi(j+1,:)' - loadedinput.input.phi(j,:)' ;
        end;
        
        Drotor      = 126.4;  % Turbine rotor diameter in (m)
        powerscale  = 1.0;    % Turbine powerscaling
        forcescale  = 1.0;    % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 250;         % Simulation length (s)
        mu       = 18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 8.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = 2;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 8;
        
    case lower('Yawcase1_2turb_100x50_lin')
        
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = 2481.9702;      % Domain length in x-direction (m)
        Ly     = 1400;           % Domain length in y-direction (m)
        Nx     = 100;             % Number of grid points in x-direction
        Ny     = 50;             % Number of grid points in y-direction
        Crx    = [400, 1281.97]; % Turbine locations in x-direction (m)
        Cry    = [700, 700];     % Turbine locations in y-direction (m)
        
        loadedinput = load([WFSimfolder 'data_SOWFA/YawCase1/system_input.mat']); % load input settings
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = loadedinput.input.beta(j,:)';
            input{j}.phi  = loadedinput.input.phi(j,:)';
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = loadedinput.input.beta(j+1,:)'- loadedinput.input.beta(j,:)';
            input{j}.dphi  = loadedinput.input.phi(j+1,:)' - loadedinput.input.phi(j,:)' ;
        end;
        
        Drotor      = 126.4;  % Turbine rotor diameter in (m)
        powerscale  = 1.75;    % Turbine powerscaling
        forcescale  = 1.2;    % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 1186;       % Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 8.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = 1;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 8;
        
    case lower('NoPrecursor_2turb_60x30_lin')
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = 2232.0623;
        Ly     = 1400;
        Nx     = 60;
        Ny     = 30;
        Crx    = [400 1032];
        Cry    = [700 700];
        
        loadedinput = load([WFSimfolder 'data_SOWFA/NoPrecursor/system_input.mat']); % load input settings
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = loadedinput.input.beta(j,:)';
            input{j}.phi  = loadedinput.input.phi(j,:)';
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = loadedinput.input.beta(j+1,:)'- loadedinput.input.beta(j,:)';
            input{j}.dphi  = loadedinput.input.phi(j+1,:)' - loadedinput.input.phi(j,:)' ;
        end;
        
        Drotor      = 126;  % Turbine rotor diameter in (m)
        powerscale  = 1.0;    % Turbine powerscaling
        forcescale  = 1.4;%1.2;    % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 10;      % Simulation length (s)
        mu       = 0*18e-5;   % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 8.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = 2;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 8;
        
    case lower('NoPrecursor_2turb_50x25_lin')
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = 2232.0623;
        Ly     = 1400;
        Nx     = 50;
        Ny     = 25;
        Crx    = [400 1032.062];
        Cry    = [700 700];
        
        loadedinput = load([WFSimfolder 'data_SOWFA/WithPrecursor/system_input.mat']); % load input settings
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = loadedinput.input.beta(j,:)';
            input{j}.phi  = loadedinput.input.phi(j,:)';
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = loadedinput.input.beta(j+1,:)'- loadedinput.input.beta(j,:)';
            input{j}.dphi  = loadedinput.input.phi(j+1,:)' - loadedinput.input.phi(j,:)' ;
        end;
        
        Drotor      = 126.3992;  % Turbine rotor diameter in (m)
        powerscale  = 1.0;       % Turbine powerscaling
        forcescale  = 1.2;       % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 750;      % Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 8.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        Wp.Turbulencemodel  = 'WFSim3';
        lmu      = 2;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 8;
        
    case lower('APC_3x3turb_noyaw_9turb_100x50_lin')
        %load([WFSimfolder 'data_SOWFA/APC/9turb_100x50_lin/workspace.mat']); % load input settings 
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = Wp.mesh.Lx;
        Ly     = Wp.mesh.Ly;
        Nx     = Wp.mesh.Nx;
        Ny     = Wp.mesh.Ny;
        Crx    = Wp.turbine.Crx;
        Cry    = Wp.turbine.Cry;
        
%         loadedinput = load([WFSimfolder 'data_SOWFA/APC/system_input.mat']); % load input settings
%         loadedinput.input.beta = [loadedinput.input.beta(:,3) loadedinput.input.beta(:,6) loadedinput.input.beta(:,1)...
%             loadedinput.input.beta(:,2) loadedinput.input.beta(:,4) loadedinput.input.beta(:,5)...
%             loadedinput.input.beta(:,7) loadedinput.input.beta(:,8) loadedinput.input.beta(:,9)];
                
        % Filter the input signals
        for j = 1:size(loadedinput.input.beta,2)
            loadedinput.input.beta(:,j)= lsim(ss(tf(1,[15 1])),loadedinput.input.beta(:,j),...
                loadedinput.input.t,loadedinput.input.beta(1,j));
        end
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = loadedinput.input.beta(j,:)';
            input{j}.phi  = loadedinput.input.phi(j,:)';
            input{j}.CT   = 4*input{j}.beta./(1+input{j}.beta).*(1-input{j}.beta./(1+input{j}.beta));
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = loadedinput.input.beta(j+1,:)'- loadedinput.input.beta(j,:)';
            input{j}.dphi  = loadedinput.input.phi(j+1,:)' - loadedinput.input.phi(j,:)' ;
        end;
        
        Drotor      = 126.3992;  % Turbine rotor diameter in (m)
        powerscale  = .45;%.55;    % Turbine powerscaling
        forcescale  = 2.5;%1.75    % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 999;       % Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 12.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
         
        Wp.Turbulencemodel  = 'WFSim3';
        lmu      = 1.5;      % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 3;
        
        % Wind farms for which no SOWFA data is available
    case lower('SingleTurbine_50x50_lin')
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = 2000;      % Domain length in x-direction (m)
        Ly     = 800;           % Domain length in y-direction (m)
        Nx     = 100;             % Number of grid points in x-direction
        Ny     = 50;             % Number of grid points in y-direction
        Crx    = 750;           % Turbine locations in x-direction (m)
        Cry    = 400;            % Turbine locations in y-direction (m)
        
        loadedinput = load([WFSimfolder 'data_SOWFA/YawCase3/system_input.mat']); % load input settings
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = 1/3;%loadedinput.input.beta(j,1)';%+.6;
            input{j}.phi  = 0*loadedinput.input.phi(j,1)';
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = loadedinput.input.beta(j+1,1)'- loadedinput.input.beta(j,1)';
            input{j}.dphi  = loadedinput.input.phi(j+1,1)' - loadedinput.input.phi(j,1)' ;
        end;
        
        Drotor      = 126.4;  % Turbine rotor diameter in (m)
        powerscale  = 1.0;    % Turbine powerscaling
        forcescale  = 2.5;    % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 50;         % Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 8.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        Wp.Turbulencemodel  = 'WFSim3';
        lmu      = 2;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 0;

    case lower('WP_CPUTime')
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = 2232.0623;
        Ly     = 1400;
        Nx     = sqrt(Wp.ll)*50*50;
        Ny     = sqrt(Wp.ll)*25*50;
        Crx    = [400 1032.062];
        Cry    = [700 700];
        
        loadedinput = load([WFSimfolder 'data_SOWFA/NoPrecursor/system_input.mat']); % load input settings
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = [.5;.5];%loadedinput.input.beta(j,:)';
            input{j}.phi  = loadedinput.input.phi(j,:)';
            input{j}.CT  = [.5;.5];
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = loadedinput.input.beta(j+1,:)'- loadedinput.input.beta(j,:)';
            input{j}.dphi  = loadedinput.input.phi(j+1,:)' - loadedinput.input.phi(j,:)' ;
        end;
        
        Drotor      = 126.3992;  % Turbine rotor diameter in (m)
        powerscale  = 1.0;    % Turbine powerscaling
        forcescale  = 1.2;    % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 5;      % Simulation length (s)
        mu       = 0*18e-5;   % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 8.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = 2;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 8;

        % Wind farms used to do MPC
    case lower('TwoTurbine_mpc')
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = 2481.9702;      % Domain length in x-direction (m)
        Ly     = 1400;           % Domain length in y-direction (m)
        Nx     = 50;             % Number of grid points in x-direction
        Ny     = 25;             % Number of grid points in y-direction
        Crx    = [400, 1281.97]; % Turbine locations in x-direction (m)
        Cry    = [700, 700];     % Turbine locations in y-direction (m)
        
        loadedinput = load([WFSimfolder 'data_SOWFA/YawCase3/system_input.mat']); % load input settings
        loadedinput.input.phi = 0*loadedinput.input.phi;
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = loadedinput.input.beta(j,:)';
            input{j}.phi  = loadedinput.input.phi(j,:)';
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = loadedinput.input.beta(j+1,:)'- loadedinput.input.beta(j,:)';
            input{j}.dphi  = loadedinput.input.phi(j+1,:)' - loadedinput.input.phi(j,:)' ;
        end;
        
        Drotor      = 126.4;  % Turbine rotor diameter in (m)
        powerscale  = 1.0;    % Turbine powerscaling
        forcescale  = 1.0;    % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 600;       % Simulation length (s)
        mu       = 18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 8.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = 2;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 8;
        
    case lower('ThreeTurbine_mpc')
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = 2000;           % Domain length in x-direction (m)
        Ly     = 600;           % Domain length in y-direction (m)
        Nx     = 50;             % Number of grid points in x-direction
        Ny     = 25;             % Number of grid points in y-direction
        Crx    = [400, 400+6*90, 400+6*90+6*90];     % Turbine locations in x-direction (m)
        Cry    = [300, 300, 300];       % Turbine locations in y-direction (m)
        
        loadedinput = load([WFSimfolder 'data_SOWFA/YawCase3/system_input.mat']); % load input settings
        loadedinput.input.phi = 0*loadedinput.input.phi;
        
        % Correctly format inputs (temporary function)
        for j = 1:length(loadedinput.input.t)
            input{j}.t    = loadedinput.input.t(j);
            input{j}.beta = [loadedinput.input.beta(j,:)';.3]+.2;
            input{j}.phi  = [loadedinput.input.phi(j,:)';0];
        end;
        
        % Calculate delta inputs
        for j = 1:length(loadedinput.input.t)-1
            input{j}.dbeta = [loadedinput.input.beta(j+1,:)'- loadedinput.input.beta(j,:)';0];
            input{j}.dphi  = [loadedinput.input.phi(j+1,:)' - loadedinput.input.phi(j,:)';0] ;
        end;
        
        Drotor      = 90;     % Turbine rotor diameter in (m)
        powerscale  = 1.0;    % Turbine powerscaling
        forcescale  = 1;    % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 350;       % Simulation length (s)
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 9.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = 1.75;         % Mixing length in x-direction (m)
        turbul   = true;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 6;
        
    case lower('ThreeTurbine_Ampc')
        type   = 'lin';          % Meshing type ('lin' or 'exp')
        Lx     = 2500;           % Domain length in x-direction (m)
        Ly     = 1000;           % Domain length in y-direction (m)
        Nx     = 50;             % Number of grid points in x-direction
        Ny     = 25;             % Number of grid points in y-direction
        Cry    =   [Ly/2; Ly/2; Ly/2];            % X-coordinate of rotor (center)
        Crx    =   [500; 500+5*90;  500+2*5*90];         % Y-coordinate of rotor (center)
        
        Drotor      = 90;     % Turbine rotor diameter in (m)
        powerscale  = 1.0;    % Turbine powerscaling
        forcescale  = 1;      % Turbine force scaling
        
        h        = 1.0;       % Sampling time (s)
        L        = 200;       % Simulation length (s)
        time     = 0:h:L;
        
        % Correctly format inputs (temporary function)
        for j = 1:length(time)
            input{j}.t    = time(j);
            input{j}.beta = [.5;.5;.5];
            input{j}.phi  = [0;0;0];
        end;
        
        % Calculate delta inputs
        for j = 1:length(time)-1
            input{j}.dbeta = input{j+1}.beta- input{j}.beta;
            input{j}.dphi  = input{j+1}.phi - input{j}.phi ;
        end;
        
        mu       = 0*18e-5;     % Dynamic flow viscosity
        Rho      = 1.20;      % Flow density (kg m-3)
        u_Inf    = 8.0;       % Freestream flow velocity x-direction (m/s)
        v_Inf    = 0.0;       % Freestream flow velocity y-direction (m/s)
        p_init   = 0.0;       % Initial values for pressure terms (Pa)
        
        lmu      = 5;         % Mixing length in x-direction (m)
        turbul   = false;      % Use mixing length turbulence model (true/false)
        n        = 2;
        m        = 4;
        
    otherwise
        error('No valid meshing specified. Please take a look at Wp.name.');
end;

%% Calculate time
time = 0:h:L;          % time vector for simulation
NN   = length(time)-1; % Total number of simulation steps

if strcmp(lower(type),'lin')
    % linear gridding
    ldx  = linspace(0,Lx,Nx);
    ldy  = linspace(0,Ly,Ny);
elseif strcmp(lower(type),'exp')
    % exponential gridding
    ldx  = expvec( Lx, cellSize(1), uniquetol(Crx,1e-2), R_con(1), N_con(1), dx_min(1) );
    ldy  = expvec( Ly, cellSize(2), uniquetol(Cry,1e-2), R_con(2), N_con(2), dx_min(2) );
else
    error('Wrong meshing type specified in "type".');
end;

Nx   = length(ldx);
Ny   = length(ldy);
ldxx = repmat(ldx',1,Ny);
ldyy = repmat(ldy,Nx,1);

% Create secondary grid from primary grid
ldx2  = 0.5*(ldx(1:end-1)+ldx(2:end));
ldx2  = [ldx2 2*ldx2(end)-ldx2(end-1)]; % add extra cells
ldy2  = 0.5*(ldy(1:end-1)+ldy(2:end));
ldy2  = [ldy2 2*ldy2(end)-ldy2(end-1)]; % add extra cells
ldxx2 = repmat(ldx2',1,Ny);
ldyy2 = repmat(ldy2,Nx,1);

% Calculate cell dimensions
dx   = diff(ldx);
dxx  = repmat([dx'; dx(end)],1,Ny);
dx2  = diff(ldx2);
dxx2 = repmat([dx2'; dx2(end)],1,Ny);
dy   = diff(ldy);
dyy  = repmat([dy, dy(end)],Nx,1);
dy2  = diff(ldy2);
dyy2 = repmat([dy2, dy2(end)],Nx,1);

% Calculate location of turbines in grid and grid mismatch
for i = 1:length(Crx)
    % Calculate cells relevant for turbine (x-dir) on primary grid
    [~,xline(i,1)] = min(abs(ldx-Crx(i))); % automatically picks earliest entry in vector
    
    % Calculate cells closest to turbines (y-dir) on both grids
    [ML_prim, L_prim ] = min(abs(ldy- (Cry(i)-Drotor/2)));
    [ML_sec , L_sec  ] = min(abs(ldy2-(Cry(i)-Drotor/2)));
    [MR_prim, R_prim ] = min(abs(ldy- (Cry(i)+Drotor/2)));
    [MR_sec , R_sec  ] = min(abs(ldy2-(Cry(i)+Drotor/2)));
    
    yline{i}  = L_prim:1: R_prim; % turbine cells for primary grid
    % ylinev{i} = L_sec :1: R_sec ; % turbine cells for secondary grid
    ylinev{i} = L_prim:1: R_prim+1; % JWs code -- Bart: I feel like this needs fixing
    
    if PrintGridMismatch
        % Calculate turbine-grid mismatch
        disp([' TURBINE ' num2str(i) ' GRID MISMATCH:']);
        disp(['                    Primary           Secondary ']);
        disp(['       center:   (' num2str(min(abs(Crx(i)-ldx)),'%10.2f/n') ',' num2str(min(abs(Cry(i)-ldy)),'%10.2f/n') ') m.      (' num2str(min(abs(Crx(i)-ldx2)),'%10.2f/n') ',' num2str(min(abs(Cry(i)-ldy2)),'%10.2f/n') ') m.']);
        disp(['   left blade:   (' num2str(min(abs(Crx(i)-ldx)),'%10.2f/n') ',' num2str(ML_prim,'%10.2f/n')             ') m.      (' num2str(min(abs(Crx(i)-ldx2)),'%10.2f/n') ',' num2str(ML_sec,'%10.2f/n')                ') m.']);
        disp(['  right blade:   (' num2str(min(abs(Crx(i)-ldx)),'%10.2f/n') ',' num2str(MR_prim,'%10.2f/n')             ') m.      (' num2str(min(abs(Crx(i)-ldx2)),'%10.2f/n') ',' num2str(MR_sec,'%10.2f/n')                ') m.']);
        disp(' ');
    end;
end;

% Calculate state sizes
Nu = (Nx-3)*(Ny-2);   % Number of u velocities in state vector
Nv = (Nx-2)*(Ny-3);   % Number of v velocities in state vector
Np = (Nx-2)*(Ny-2)-2; % Number of pressure terms in state vector

%% Display results
if plotMesh
    clf;
    Z1 = -2*ones(size(ldxx));
    mesh(ldyy,ldxx,Z1,'FaceAlpha',0,'EdgeColor','black','LineWidth',0.1)
    hold on;
    Z2 = -1*ones(size(ldxx2));
    mesh(ldyy2,ldxx2,Z2,'FaceAlpha',0,'EdgeColor','blue','LineWidth',0.1)
    hold on;
    plot([Cry-Drotor/2;Cry+Drotor/2],[Crx;Crx],'r-','LineWidth',2.0)
    axis equal
    xlim([-0.1*Ly 1.2*Ly]);
    ylim([-0.1*Lx 1.1*Lx]);
    xlabel('y (m)');
    ylabel('x (m)');
    view(0,90); % view from the top
    drawnow;
end;

%% Export to Wp and input
Wp         = struct('Nu',Nu,'Nv',Nv,'Np',Np,'name',Wp.name,'Turbulencemodel',Wp.Turbulencemodel);
Wp.sim     = struct('h',h,'time',time,'L',L,'NN',NN);
Wp.turbine = struct('Drotor',Drotor,'powerscale',powerscale,'forcescale',forcescale, ...
    'N',length(Crx),'Crx',Crx,'Cry',Cry);
Wp.site    = struct('mu',mu,'Rho',Rho,'u_Inf',u_Inf,'v_Inf',v_Inf,'p_init',p_init, ...
    'lmu',lmu,'turbul',turbul,'m',m,'n',n);
Wp.mesh    = struct('Lx',Lx,'Ly',Ly,'Nx',Nx,'Ny',Ny,'ldxx',ldxx,'ldyy',ldyy,'ldxx2',...
    ldxx2,'ldyy2',ldyy2,'dxx',dxx,'dyy',dyy,'dxx2',dxx2,'dyy2',dyy2,...
    'xline',xline,'type',type);Wp.mesh.yline = yline; Wp.mesh.ylinev = ylinev; % Do not support struct command


%% Construct mu if no turbulence
if turbul==0; ConstructMu; Wp.site(:).mu = mu; end

end