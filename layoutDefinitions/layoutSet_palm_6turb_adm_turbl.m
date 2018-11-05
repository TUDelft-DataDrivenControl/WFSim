function [Wp] = layoutSet_palm_6turb_adm_turbl()

Wp = struct('description','6 NREL 5MW turbines case, turbulent inflow, based on a PALM ADM simulation');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',false ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
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
    'lm_slope',0.05,... % Mixing length in x-direction (m)
    'd_lower',73.3,... % Turbulence model gridding property
    'd_upper',601.9,... % Turbulence model gridding property
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'Lx',1740,... % Domain length in x-direction
    'Ly',748,... % Domain length in y-direction
    'Nx',80,... % Number of cells in x-direction
    'Ny',40 ... % Number of cells in y-direction
    );

end