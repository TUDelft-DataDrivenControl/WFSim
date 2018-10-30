clear all; close all; clc;
% This function allows one to compare WFSim to time-averaged LES data for a single turbine simulation under various yaw angles.
% LES data is included for two cases: one with laminar inflow, the other with turbulent inflow, both generated using SOWFA.

% Load LES data
% load('uniformInflow_yawRange.mat');
load('turbInflow_yawRange.mat');

% Set up WFSim
Wp = layoutSet_sowfa_1turb(); % Choose which scenario to simulate. See 'layoutDefinitions' folder for the full list.
addpath('../solverDefinitions'); % Folder with model options, solver settings, etc.
modelOptions = solverSet_default(Wp); % Choose model solver options

% Create interpolation objects
F_u = scatteredInterpolant(cellCenters(:,1),cellCenters(:,2),cellData(1).cellData(:,1),'nearest');
F_v = scatteredInterpolant(cellCenters(:,1),cellCenters(:,2),cellData(1).cellData(:,2),'nearest');

for i = 1:length(cellData)
    % Control setting
    yawRange = [30:-10:-30];
    turbInput = struct('t',0,'CT_prime',2,'phi',yawRange(i));

    %% Script core functions
    run('../WFSim_addpaths.m'); % Add essential paths to MATLABs environment
    [Wp,sol,sys] = InitWFSim(Wp,modelOptions,0); % Initialize WFSim model

    % Propagate the WFSim model
    [sol,sys] = WFSim_timestepping(sol,sys,Wp,turbInput,modelOptions); % forward timestep: x_k+1 = f(x_k)

    % Compare the two
    F_u.Values = cellData(i).cellData(:,1);
    F_v.Values = cellData(i).cellData(:,2);
    u_LES = F_u(Wp.mesh.ldxx2,Wp.mesh.ldyy );
    v_LES = F_v(Wp.mesh.ldxx, Wp.mesh.ldyy2);

    figure('Position',[1 39.4000 1536 752]);
    subplot(2,3,1);
    contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',sol.u,(0:0.1:Wp.site.u_Inf*1.2),'Linecolor','none'); 
    caxis([min(min(sol.u))-2 Wp.site.u_Inf*1.04]); axis equal;
    xlim([1000 2000]); ylim([500 3000]);
    colorbar;
    subplot(2,3,2);
    contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',u_LES,(0:0.1:Wp.site.u_Inf*1.2),'Linecolor','none');  
    caxis([min(min(sol.u))-2 Wp.site.u_Inf*1.04]); axis equal;
    xlim([1000 2000]); ylim([500 3000]);
    colorbar;
    subplot(2,3,3);
    contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',abs(sol.u-u_LES),(0:0.1:3),'Linecolor','none'); 
    caxis([0 3]); axis equal;
    xlim([1000 2000]); ylim([500 3000]);
    colorbar;
    subplot(2,3,4);
    contourf(Wp.mesh.ldyy2(1,:),Wp.mesh.ldxx(:,1)',sol.v,(-3:0.1:3),'Linecolor','none');  
    caxis([min(min(sol.v))-1 max(max(sol.v))+1]); axis equal;
    xlim([1000 2000]); ylim([500 3000]);
    colorbar;
    subplot(2,3,5);
    contourf(Wp.mesh.ldyy2(1,:),Wp.mesh.ldxx(:,1)',v_LES,(-3:0.1:3),'Linecolor','none');  
    caxis([min(min(sol.v))-1 max(max(sol.v))+1]); axis equal;
    xlim([1000 2000]); ylim([500 3000]);
    colorbar;
    subplot(2,3,6);
    contourf(Wp.mesh.ldyy2(1,:),Wp.mesh.ldxx(:,1)',abs(sol.v-v_LES),(0:0.1:3),'Linecolor','none'); 
    caxis([0 3]); axis equal;
    xlim([1000 2000]); ylim([500 3000]);
    colorbar;
    
    drawnow;
end

function [Wp] = layoutSet_sowfa_1turb()
Wp = struct('description','1 NREL 5MW turbine case, uniform inflow, based on a SOWFA ALM simulation');
Wp.sim = struct(...
    'h',Inf,... % timestep (s)
    'startUniform',false ... % Start from a uniform flow field (T) or from a fully developed waked flow field (F).
    );
Wp.turbine = struct(...
    'Crx',[1000],... % X-coordinates of turbines (m)
    'Cry',[1500],... % Y-coordinates of turbines (m)
    'Drotor',126.4*0.9,... % Rotor diameter (m), note that WFSim only supports a uniform Drotor for now
    'powerscale',0.95*(1/0.9),... % Turbine power scaling
    'forcescale',1.40 ... % Turbine force scaling
    );
Wp.site = struct(...
    'u_Inf',8.0,... % Initial long. wind speed in m/s
    'v_Inf',0.0,... % Initial lat. wind speed in m/s
    'p_init',0.0,... % Initial values for pressure terms (Pa)
    'lm_slope',0.07,... % Mixing length in x-direction (m)
    'd_lower',73.3,... % Turbulence model gridding property
    'd_upper',601.9,... % Turbulence model gridding property
    'Rho',1.20 ... % Air density
    );
Wp.mesh = struct(...
    'Lx',3000,... % Domain length in x-direction
    'Ly',3000,... % Domain length in y-direction
    'Nx',100,... % Number of cells in x-direction
    'Ny',100 ... % Number of cells in y-direction
    );
end