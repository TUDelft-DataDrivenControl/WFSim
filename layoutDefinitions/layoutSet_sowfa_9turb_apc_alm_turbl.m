function [Wp] = layoutSet_sowfa_9turb_apc_alm_turbl()

Wp = struct('description','9 NREL 5MW turbines case, based on a SOWFA ALM simulation');

Wp.sim = struct(...
    'h',1.0,... % timestep (s)
    'startUniform',true ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );

Wp.turbine = struct(...
    'Crx',[0.4048 0.4024 0.40 1.0368 1.0344 1.0320 1.6688 1.6663 1.6639]*1e3,... % X-coordinates of turbines (m)
    'Cry',[1.1584 0.7792 0.40 1.1543 0.7752 0.3960 1.1503 0.7711 0.3919]*1e3,... % Y-coordinates of turbines (m)
    'Drotor',126.4,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.99,... % Turbine power scaling
    'forcescale',1.9 ... % Turbine force scaling
    );

Wp.site = struct(...
    'u_Inf',12.0214,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'lm_slope',0.05,... % Mixing length in x-direction (m)
    'd_lower',140.0,... % Turbulence model gridding property
    'd_upper',1000.0,... % Turbulence model gridding property
    'Rho',1.20 ... % Air density
    );

Wp.mesh = struct(...
    'Lx',2518.8,... % Domain length in x-direction
    'Ly',1558.4,... % Domain length in y-direction
    'Nx',100,... % Number of cells in x-direction
    'Ny',42 ... % Number of cells in y-direction
    );

% Tuning notes 'apc_9turb_alm_turb' (Sep 11th, 2017):
% Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.5, m = 1:8, n = 1:4

% Tuning notes on Nov 5, 2018
% Optimal out of:
%     'forcescale', 1.5:0.1:2.5, ...
%     'lm_slope', 0.005:0.005:0.10, ...
%     'd_lower', 0.1:20:200.1,...
%     'd_upper', 300:50:1000,...
end