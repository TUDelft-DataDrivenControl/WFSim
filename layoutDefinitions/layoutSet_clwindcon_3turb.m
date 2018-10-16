function [Wp] = layoutSet_clwindcon_3turb()
Wp = struct('description','3 DTU 10MW turbines case for the CL-Windcon project');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[1. 6. 11.]*178.3,... % X-coordinates of turbines (m)
    'Cry',[2.5 2.5 3.]*178.3,... % Y-coordinates of turbines (m)
    'Drotor',178.3,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.95,... % Turbine power scaling
    'forcescale',1.50 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',8.0,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'lm_slope',0.05,... % Mixing length in x-direction (m)
    'd_lower',50,... % Turbulence model gridding property
    'd_upper',700,... % Turbulence model gridding property
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'Lx',15*178.3,... % Domain length in x-direction
    'Ly', 5*178.3,... % Domain length in y-direction
    'Nx',100,... % Number of cells in x-direction
    'Ny',40 ... % Number of cells in y-direction
    );
end