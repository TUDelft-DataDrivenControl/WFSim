% This function checks by running different simulations if the power
% expression is correct.
% 1) Single turbine, power as a function of axial induction
% 2) Single turbine, power as a function of yaw angle
% 3) Two turbines partial overlap, power as a function of changing yaw angle turbine 1 
% 4) Comparison with SOWFA simulation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name             = 'SingleTurbine_50x50_lin';   % Meshing name (see "\bin\core\meshing.m")
Wp.Turbulencemodel  = 'WFSim3';

Animate       = 0;                      % Show 2D flow fields every x iterations (0: no plots)
plotMesh      = 0;                      % Show meshing and turbine locations
conv_eps      = 1e-6;                   % Convergence threshold
max_it_dyn    = 1;                      % Maximum number of iterations for k > 1

if options.startUniform==1
    max_it = 1;
else
    max_it = 30;
end

% WFSim general initialization script
[Wp,sol,sys,Power,CT,a,Ueffect,input,B1,B2,bc] ...
    = InitWFSim(Wp,options,plotMesh);

if Animate > 0
    scrsz = get(0,'ScreenSize');
    hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
end

Beta       = linspace(0,1,35); 
Power_ak   = zeros(1,length(Beta));

for l=1:length(Beta)
    
    input{1}.beta = Beta(l);
        
    %% Loop
    for k=1:Wp.sim.NN
        tic
        it        = 0;
        eps       = 1e19;
        epss      = 1e20;
        
        
        while ( eps>conv_eps && it<max_it && eps<epss );
            it   = it+1;
            epss = eps;
                    
            [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k)] = ...
                Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options);                  % Create system matrices
            [sol,sys] = Computesol(sys,input{k},sol,k,it,options);                  % Compute solution
            [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);        % Map solution to field
            
        end
        toc
        
        if Animate > 0
            if ~rem(k,Animate)
                Animation;
            end;
        end;
    end;
    
    Power_ak(l) = Power(:,k);
    ak(l)       = a(:,k);
    
end

disp('Completed simulations.');

figure(2);clf
plot(ak,Power_ak/max(Power_ak));grid
ylabel('$\overline{P}_1$','interpreter','latex');xlabel('$a$','interpreter','latex')

disp(' ')
disp('Hit a key for following test')
disp(' ')
pause
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name       = 'SingleTurbine_50x50_lin';   % Meshing name (see "\bin\core\meshing.m")
Wp.Turbulencemodel  = 'WFSim3';


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

Phi       = linspace(-40,40,50);
Power_Phi = zeros(1,length(Phi));

for l=1:length(Phi)
    
    input{1}.phi = Phi(l);
    
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
            
            %display(['k ',num2str(k,'%-1000.1f'),', It ',num2str(it,'%-1000.0f'),', Nv=', num2str(Normv{k}(it),'%10.2e'), ', Nu=', num2str(Normu{k}(it),'%10.2e'), ', TN=',num2str(eps,'%10.2e'),', Np=','Mean effective=',num2str(mean(Ueffect(1,k)),'%-1000.2f')]) ;
        end
        toc
        
        if Animate > 0
            if ~rem(k,Animate)
                Animation;
            end;
        end;
    end;
    
    Power_Phi(l) = Power(:,k);
end

disp('Completed simulations.');

figure(2);clf
plot(Phi,Power_Phi/max(Power_Phi));grid
ylabel('Power');xlabel('\phi')

disp(' ')
disp('Hit a key for following test')
disp(' ')
pause
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name       = 'TwoTurbinePartialOverlap_lin';   % Meshing name (see "\bin\core\meshing.m")
Wp.Turbulencemodel  = 'WFSim3';

Animate       = 50;                     % Show 2D flow fields every x iterations (0: no plots)
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

% Loop
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
        
        %display(['k ',num2str(k,'%-1000.1f'),', It ',num2str(it,'%-1000.0f'),', Nv=', num2str(Normv{k}(it),'%10.2e'), ', Nu=', num2str(Normu{k}(it),'%10.2e'), ', TN=',num2str(eps,'%10.2e'),', Np=','Mean effective=',num2str(mean(Ueffect(1,k)),'%-1000.2f')]) ;
    end
    toc
    
    if Animate > 0
        if ~rem(k,Animate)
            Animation;
        end;
    end;
end;

disp('Completed simulations.');

m = sum(Power(:,1));
figure(2);clf
plot(Wp.sim.time(1:end-1),Power(1,:)/m); hold on
plot(Wp.sim.time(1:end-1),Power(2,:)/m,'r');
plot(Wp.sim.time(1:end-1),sum(Power(:,:))/m,'b');grid;
ylabel('Power [W]'); xlabel('Time [s]')
legend('Power T_1','Power T_2','Total Power')

disp(' ')
disp('Hit a key for following test')
disp(' ')
pause
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; clc; close all;

%% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

Wp.name       = 'YawCase3_50x50_lin';   % Meshing name (see "\bin\core\meshing.m")
Wp.Turbulencemodel  = 'WFSim3';

Animate       = 0;                     % Show 2D flow fields every x iterations (0: no plots)
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
            Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options);                  % Create system matrices
        [sol,sys] = Computesol(sys,input{k},sol,k,it,options);                  % Compute solution
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);        % Map solution to field
        
        %display(['k ',num2str(k,'%-1000.1f'),', It ',num2str(it,'%-1000.0f'),', Nv=', num2str(Normv{k}(it),'%10.2e'), ', Nu=', num2str(Normu{k}(it),'%10.2e'), ', TN=',num2str(eps,'%10.2e'),', Np=','Mean effective=',num2str(mean(Ueffect(1,k)),'%-1000.2f')]) ;
    end
    toc
    
    if Animate > 0
        if ~rem(k,Animate)
            Animation;
        end;
    end;
end;

disp('Completed simulations.');

load('..\..\Data_SOWFA\YawCase3\Power_YawCase3.mat');
        
figure(2);clf
plot(Wp.sim.time(1:end-1),Power(1,:)); hold on
plot(Wp.sim.time(1:end-1),Power(2,:),'r');
plot(Wp.sim.time(1:end-1),Power_SOWFA(1,1:Wp.sim.NN),'b--');
plot(Wp.sim.time(1:end-1),Power_SOWFA(2,1:Wp.sim.NN),'r--');grid;
ylabel('Power [W]'); xlabel('Time [s]');
legend('T_1 (WFSim)','T_2 (WFSim)','T_1 (SOWFA)','T_2 (SOWFA)')

disp(' ')
disp('End of tests')
disp(' ')