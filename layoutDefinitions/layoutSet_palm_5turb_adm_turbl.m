function [Wp] = layoutSet_palm_5turb_adm_turbl()

Wp = struct('description','5 NREL 5MW turbines case, turbulent inflow, based on a PALM ADM simulation');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[250,880,1510,2140,2770],... % X-coordinates of turbines (m)
    'Cry',[200,200,200,200,200],... % Y-coordinates of turbines (m)
    'Drotor',126,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.99,... % Turbine power scaling
    'forcescale',1.9 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',7.5,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'lm_slope',0.05,... % Mixing length in x-direction (m)
    'd_lower',140.0,... % Turbulence model gridding property
    'd_upper',1000.0,... % Turbulence model gridding property
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'Lx',3500,... % Domain length in x-direction
    'Ly',400,... % Domain length in y-direction
    'Nx',80,... % Number of cells in x-direction
    'Ny',40 ... % Number of cells in y-direction
    );

% Tuning notes '5turb_adm_turbl' (Jan 6th, 2019):
% Manually tuned

end