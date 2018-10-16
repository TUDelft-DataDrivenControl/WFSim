function [Wp] = layoutSet_clwindcon_windTunnel()
Wp = struct('description','2 TUM G1 turbines for the CL-Windcon project');

Wp.sim = struct(...
    'h',0.03,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[2.0 6.0],... % X-coordinates of turbines (m)
    'Cry',[2.0 2.5],... % Y-coordinates of turbines (m)
    'Drotor',1.0,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',1.0,... % Turbine power scaling
    'forcescale',1.50 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',5.7,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'lm_slope',0.10,... % Mixing length in x-direction (m)
    'd_lower',0.1,... % Turbulence model gridding property
    'd_upper',3.0,... % Turbulence model gridding property
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'Lx',10.0,... % Domain length in x-direction
    'Ly',5.0,... % Domain length in y-direction
    'Nx',40,... % Number of cells in x-direction
    'Ny',20 ... % Number of cells in y-direction
    );
end