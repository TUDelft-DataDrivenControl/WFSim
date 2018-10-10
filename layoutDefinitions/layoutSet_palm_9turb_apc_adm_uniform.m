function [Wp] = layoutSet_palm_9turb_apc_adm_uniform()

Wp = struct('description','9 NREL 5MW turbines case, uniform inflow, based on a PALM ADM simulation');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[repmat([300 930 1560],1,3)],... % X-coordinates of turbines (m)
    'Cry',[repmat(300,1,3) repmat(680,1,3) repmat(1060,1,3)],... % Y-coordinates of turbines (m)
    'Drotor',120.0,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.90,... % Turbine power scaling
    'forcescale',1.60 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',7.856,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'turbul',true,... % Use mixing length turbulence model (true/false)
    'turbModel','WFSim3',...  % Turbulence model of choice
    'lmu',0.30,... % Mixing length in x-direction (m)
    'm',4,... % Turbulence model gridding property
    'n',2,... % Turbulence model gridding property
    'mu',0.0,... % Dynamic flow viscosity
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'gridType','lin',... % Grid type ('lin' the only supported one currently)
    'Lx',2800,... % Domain length in x-direction
    'Ly',1360,... % Domain length in y-direction
    'Nx',100,... % Number of cells in x-direction
    'Ny',50 ... % Number of cells in y-direction
    );

% Tuning notes 'apc_9turb_adm_noturb' (Sep 10th, 2017):
% Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
end