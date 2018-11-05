function [Wp] = layoutSet_palm_2turb_adm_uniform()

Wp = struct('description','2 NREL 5MW turbines case, uniform inflow, based on a PALM ADM simulation');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',false ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[180 936],... % X-coordinates of turbines (m)
    'Cry',[255 255],... % Y-coordinates of turbines (m)
    'Drotor',110.0,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',1.0,... % Turbine power scaling
    'forcescale',1.7 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',7.986,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'lm_slope',0.01,... % Mixing length in x-direction (m)
    'd_lower',73.3,... % Turbulence model gridding property
    'd_upper',601.9,... % Turbulence model gridding property
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'Lx',1686,... % Domain length in x-direction
    'Ly',500,... % Domain length in y-direction
    'Nx',40,... % Number of cells in x-direction
    'Ny',20 ... % Number of cells in y-direction
    );


% Tuned by hand on November 5, 2018. There may be some discrepancies
% at the rotor/wake edge due to mesh mismatches. This is unavoidable and
% would be significantly less noticeable under turbulence and mixing. In
% general, the dynamics are very well captured.
end