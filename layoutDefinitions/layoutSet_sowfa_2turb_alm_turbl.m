function [Wp] = layoutSet_sowfa_2turb_alm_turbl()

Wp = struct('description','2 NREL 5MW turbines case, turbulent inflow, based on a SOWFA ALM simulation');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[0.4000    1.0321]*1e3,... % X-coordinates of turbines (m)
    'Cry',[400.0000  398.5747],... % Y-coordinates of turbines (m)
    'Drotor',126.4,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.95,... % Turbine power scaling
    'forcescale',1.40 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',8.0446,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'turbul',true,... % Use mixing length turbulence model (true/false)
    'lm_slope',0.1,... % Mixing length in x-direction (m)
    'd_lower',73.3,... % Turbulence model gridding property
    'd_upper',601.9,... % Turbulence model gridding property
    'mu',0.0,... % Dynamic flow viscosity
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'Lx',1882.1,... % Domain length in x-direction
    'Ly',800.0,... % Domain length in y-direction
    'Nx',50,... % Number of cells in x-direction
    'Ny',25 ... % Number of cells in y-direction
    );

% Tuning notes '2turb_alm_turb' (Sep 6th, 2017):
% Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
end