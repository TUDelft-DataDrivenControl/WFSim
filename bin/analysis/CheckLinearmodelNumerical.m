%% Compare the linear model with the nonlinear model
% This script to simulates linear model and the nonlinear model
% with equivalent delta perturbations. Note: Not finished yet!!
clear; clc; close all;

%% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 1;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 1;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name       = 'Reinier_2turbine_20x12';       % Meshing name (see "\bin\core\meshing.m")

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
end

%% Check linearization numerical
% sys.A * x_{k+1}  = sys.M * x_k + sys.m(x_k,u_k) is nonlinear model and 
% sol.dx  = sys.A\sys.Al*sol.dx+sys.bl;
% the linear model we will check with sys.Al = dcdx+dbcdx+dSm.dx-Al 
% and sys.bl =[vec(Sm.dx');vec(Sm.dy');0.*bc];
       
x0     = sol.x;

[sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k)] = ...
                    Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options); % Create system matrices
                
E0     = sys.A;
B0     = sys.m;
A0     = sys.M;
As0    = E0\A0;
Bs0    = E0\B0;

deltax = .01;
I      = eye(size(sys.M,1));

for kk=1:size(sys.M,1)
    
    % Small perturbation on one state
    sol.x     = x0 + I(:,kk)*deltax;
    
    % Construct analytic linear model with perturbed state
    [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);    % Map solution to field
    [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k)] = ...
        Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options);              % Create system matrices
    [sol,sys] = Computesol(sys,input{k},sol,k,it,options);              % Compute solution
    
    El     = sys.A;
    Al     = sys.Al;
    Bl     = sys.bl;   
    Asl    = El\Al;
    Bsl    = El\Bl;
    
    
    % Update nonlinear model with perturbed state
    E3     = sys.A;
    B3     = sys.m;
    A3     = sys.M;
    
    % Numerical differentiate nonlinear model
    dEdx_numeric(:,kk)  = (E3-E0)*sol.x./deltax;       % Sould be equal to dEdx
    dBdx_numeric(:,kk)  = (B3-B0)./deltax;             % Sould be equal to dBdx     
    
    Al_numeric(:,kk)    = dEdx_numeric(:,kk)+dBdx_numeric(:,kk);
end

[dEdx_numeric dEdx]
[dBdx_numeric dBdx]
[As_numeric Asl]

% Perturbation on u
x0     = x;
u0     = u;
E0     = [x0(1) x0(2); x0(2) x0(2)];
B0     = [x0(1)^2*cos(u0); x0(2)^2*sin(u0)];
A0     = A;
As0    = E0\A0;
Bs0    = E0\B0;

deltau = .0001;
I      = eye(size(Bl,2));

for kk=1:size(Bl,2)
    
    % Small perturbation on input
    x1     = x0;
    u1     = u0 + I(:,kk)*deltau;
    
    % Construct analytic linear model with perturbed state
    dEdx   = [x1(1) x1(2);0 x1(1)+x1(2)];
    dBdx   = [2*x1(1)*cos(u1) 0;0 2*x1(2)*sin(u1)];
    El     = [x1(1) x1(2); x1(2) x1(2)];
    Al     = A+dBdx-dEdx;
    Bl     = [-x1(1)^2*sin(u1); x1(2)^2*cos(u1)];   % dBdu
    Asl    = El\Al;
    Bsl    = El\Bl; 
    
    
    % Update nonlinear model with perturbed state
    E1     = [x1(1) x1(2); x1(2) x1(2)];
    B1     = [x1(1)^2*cos(u1); x1(2)^2*sin(u1)];
    A1     = A0;
    As1    = E1\A1; 
    Bs1    = E1\B1;
    
    % Numerical differentiate nonlinear model
    dBdu_numeric(:,kk)  = (B1-B0)./deltau;        % Sould be equal to Bl = dBdu 
    Bs_numeric(:,kk)    = (Bs1-Bs0)./deltau;      % Sould be equal to Bsl
        
end

[dBdu_numeric Bl]
[Bs_numeric Bsl]


