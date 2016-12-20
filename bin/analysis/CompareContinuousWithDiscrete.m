%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Types of systems when Projection=0:
% Linearversion = 0;
% sys.A * x_{k+1}  = sys.b(x_k,u_k)
% sys.A * x_{k+1}  = sys.M * x_k + sys.m(x_k,u_k)
% Linearversion = 1;
% sys.A * dx_{k+1} = sys.Al * dx_k + sys.bl(du_k)
% sys.A * dx_{k+1} = sys.Al * dx_k + sys.Bl * du_k
% sys.Ac * \dot{x} = sys.Alc * dx + sys.Bl * u

clear; clc; close all;

%% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 1;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name       = 'TwoTurbine_Dev_lin';   % Meshing name (see "\bin\core\meshing.m")

Animate       = 0;                      % Show 2D flow fields every x iterations (0: no plots)
plotMesh      = 0;                      % Show meshing and turbine locations
conv_eps      = 1e-6;                   % Convergence threshold
max_it_dyn    = 1;                      % Maximum number of iterations for k > 1

if options.startUniform==1
    max_it = 1; 
else
    max_it = 50;
end

% WFSim general initialization script
[Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] ...
    = InitWFSim(Wp,options,plotMesh);

if Animate > 0
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
           [0 0 1 1],'ToolBar','none','visible', 'on');
end


%% Loop
for k=1:Wp.sim.NN
    tic
    it        = 0;
    eps       = 1e19;
    epss      = 1e20;
    
    while ( eps>conv_eps && it<max_it && eps<epss );
        it   = it+1;
        epss = eps;
        
        if k>1 
            max_it = max_it_dyn; 
        end
        
        [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k)] = ...
                    Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options); % Create system matrices
        [sol,sys] = Computesol(sys,input{k},sol,k,it,options);                   % Compute solution
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);         % Map solution to field
        
    end
    toc
    
    if Animate > 0
        if ~rem(k,Animate)
            Animation; 
        end; 
    end; 
end;

%%
C         = zeros(2,Wp.Nu+Wp.Nv+Wp.Np);
C(1,(Wp.mesh.xline(1)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{1}(1)-1:...
    (Wp.mesh.xline(1)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{1}(end)-1) = 1/2; % Extracts sum of rotor velocities first turbine
C(2,(Wp.mesh.xline(2)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{2}(1)-1:...
    (Wp.mesh.xline(2)-3)*(Wp.mesh.Ny-2)+Wp.mesh.yline{2}(end)-1) = 1/2; % Extracts sum of rotor velocities second turbine

dsys = dss(full(sys.Al),full(sys.Bl),C,0,full(sys.A),1);
csys = d2c(dsys,'tustin'); 
size(dsys)

nu = 1;
ny = 1;
figure(1);clf
step(dsys(ny,nu));hold on;
step(csys(ny,nu));


