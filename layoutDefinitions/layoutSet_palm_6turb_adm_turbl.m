function [Wp] = layoutSet_palm_6turb_adm_turbl()

Wp = struct('description','6 NREL 5MW turbines case, turbulent inflow, based on a PALM ADM simulation');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[180 180 810 810 1440 1440],... % X-coordinates of turbines (m)
    'Cry',[185 563 185 563 185 563],... % Y-coordinates of turbines (m)
    'Drotor',120.0,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.95,... % Turbine power scaling
    'forcescale',1.50 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',7.705,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'turbul',true,... % Use mixing length turbulence model (true/false)
    'turbModel','WFSim3',...  % Turbulence model of choice
    'lmu',0.60,... % Mixing length in x-direction (m)
    'm',4,... % Turbulence model gridding property
    'n',2,... % Turbulence model gridding property
    'mu',0.0,... % Dynamic flow viscosity
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'gridType','lin',... % Grid type ('lin' the only supported one currently)
    'Lx',1740,... % Domain length in x-direction
    'Ly',748,... % Domain length in y-direction
    'Nx',80,... % Number of cells in x-direction
    'Ny',40 ... % Number of cells in y-direction
    );

%         d_lower  = 73.3;  % Turbulence model gridding property (WES: d')
%         d_upper  = 601.9; % Turbulence model gridding property (WES: d)
%         lm_slope = 0.068; % Turbulence model gridding property (WES: l_s)
        % Ranges: lmu= xxx, f = xxx, m = xxx, n = xxx  
end