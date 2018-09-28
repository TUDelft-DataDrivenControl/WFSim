function [Wp] = palm_2turb_yaw_adm_uniform()

Wp = struct('description','2 NREL 5MW turbines case, uniform inflow, based on a PALM ADM simulation in which turbine 1 is yawed');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[180 936],... % X-coordinates of turbines (m)
    'Cry',[250 250],... % Y-coordinates of turbines (m)
    'Drotor',120.0,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.95,... % Turbine power scaling
    'forcescale',1.60 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',7.923,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'turbul',true,... % Use mixing length turbulence model (true/false)
    'turbModel','WFSim3',...  % Turbulence model of choice
    'lmu',0.60,... % Mixing length in x-direction (m)
    'm',1,... % Turbulence model gridding property
    'n',4,... % Turbulence model gridding property
    'mu',0.0,... % Dynamic flow viscosity
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'gridType','lin',... % Grid type ('lin' the only supported one currently)
    'Lx',1686,... % Domain length in x-direction
    'Ly',500,... % Domain length in y-direction
    'Nx',50,... % Number of cells in x-direction
    'Ny',25 ... % Number of cells in y-direction
    );


% Tuning notes '2turb_yaw_adm_noturb' (Sep 7th, 2017):
% Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
end