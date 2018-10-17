function [Wp] = layoutSet_palm_2turb_adm_uniform()

Wp = struct('description','2 NREL 5MW turbines case, uniform inflow, based on a PALM ADM simulation');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[180 936],... % X-coordinates of turbines (m)
    'Cry',[250 250],... % Y-coordinates of turbines (m)
    'Drotor',120.0,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',1.0,... % Turbine power scaling
    'forcescale',1.7 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',7.986,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'lm_slope',0.1*ones(1,size(Wp.turbine.Crx,2)),... % Mixing length in x-direction (m)
    'd_lower',73.3*ones(1,size(Wp.turbine.Crx,2)),... % Turbulence model gridding property
    'd_upper',601.9*ones(1,size(Wp.turbine.Crx,2)),... % Turbulence model gridding property
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'Lx',1686,... % Domain length in x-direction
    'Ly',500,... % Domain length in y-direction
    'Nx',40,... % Number of cells in x-direction
    'Ny',20 ... % Number of cells in y-direction
    );


% Tuning notes '2turb_adm_turb', '2turb_adm_noturb' (Sep 6th, 2017):
% Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
% Note:   Gridsearched multiple ways, should be very reliable now.
% Results for both cases are so similar that they have been merged.
end