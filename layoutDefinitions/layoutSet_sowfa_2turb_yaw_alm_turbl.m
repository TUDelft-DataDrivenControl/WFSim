function [Wp] = layoutSet_sowfa_2turb_yaw_alm_turbl()
error('This case has not yet been tuned/validated. Remove this message manually to continue.');

Wp = struct('description','2 NREL 5MW turbines case, turbulent inflow, based on a SOWFA ALM simulation where turbine 1 is yawed');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[0.4000    1.2820]*1e3,... % X-coordinates of turbines (m)
    'Cry',[400.0000  397.6982],... % Y-coordinates of turbines (m)
    'Drotor',126.4,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',1.0,... % Turbine power scaling
    'forcescale',1.20 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',8.0641,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'lm_slope',0.1*ones(1,size(Wp.turbine.Crx,2)),... % Mixing length in x-direction (m)
    'd_lower',73.3*ones(1,size(Wp.turbine.Crx,2)),... % Turbulence model gridding property
    'd_upper',601.9*ones(1,size(Wp.turbine.Crx,2)),... % Turbulence model gridding property
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'Lx',2132.0,... % Domain length in x-direction
    'Ly',800.0,... % Domain length in y-direction
    'Nx',50,... % Number of cells in x-direction
    'Ny',25 ... % Number of cells in y-direction
    );
end