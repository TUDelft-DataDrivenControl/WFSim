function [ Wp ] = meshing( scenarioName, plotMesh, PrintGridMismatch )
%MESHING Meshing and settings function for the WFSim code
% This code includes all the topology information, atmospheric
% information, turbine properties and turbine control settings for any
% simulation scenario. Basically, all wind farm-related information is
% defined here.

% Default settings
if nargin <= 0; error('Please specify a meshing case.'); end;
if nargin <= 1; plotMesh = true;                         end;
if nargin <= 2; PrintGridMismatch = true;                end;

% Some pre-processing to determine correct input file location
binloc  = which('meshing.m');
if ispc; slashSymbol = '\'; else; slashSymbol = '/'; end; % Compatible with both UNIX and PC
strbslash   = strfind(binloc ,[slashSymbol 'WFSim']);
WFSimfolder = binloc(1:strbslash(end)+6);

switch lower(scenarioName)    
    
    % Wind farms for which PALM data is available        
    case{lower('2turb_adm_turb'),lower('2turb_adm_noturb')}
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = 1.0;        % Turbine power scaling
        forcescale = 1.7;        % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = true;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice
        lmu        = 0.45;       % Mixing length in x-direction (m)
        mu         = 0*18e-5;    % Dynamic flow viscosity        
        m          = 1;          % Turbulence model gridding property        
        n          = 4;          % Turbulence model gridding property
                         
        % Tuning notes '2turb_adm_turb', '2turb_adm_noturb' (Sep 6th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
        % Note:   Gridsearched multiple ways, should be very reliable now.
        % Results for both cases are so similar that they have been merged.
        
    case lower('2turb_yaw_adm_noturb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = .95;        % Turbine power scaling
        forcescale = 1.6;        % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = true;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice   
        lmu        = 0.6;        % Mixing length in x-direction (m)
        mu         = 0.0;        % Dynamic flow viscosity
        m          = 1;          % Turbulence model gridding property        
        n          = 4;          % Turbulence model gridding property
        
        % Tuning notes '2turb_yaw_adm_noturb' (Sep 7th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
    
    case lower('6turb_adm_turb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = .95;        % Turbine power scaling
        forcescale = 1.5;        % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = false;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice   
        lmu        = 0.6;        % Mixing length in x-direction (m)
        mu         = 0.0;        % Dynamic flow viscosity
        m          = 4;          % Turbulence model gridding property        
        n          = 2;          % Turbulence model gridding property        
        
        % Tuning notes '6turb_adm_turb' (Oct 06th, 2017): 
        % Ranges: lmu= xxx, f = xxx, m = xxx, n = xxx
    
    case lower('apc_9turb_adm_noturb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = .90;        % Turbine power scaling
        forcescale = 1.6;        % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = true;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice   
        lmu        = 0.3;        % Mixing length in x-direction (m)
        mu         = 0.0;        % Dynamic flow viscosity
        m          = 4;          % Turbulence model gridding property        
        n          = 2;          % Turbulence model gridding property        
        
        % Tuning notes 'apc_9turb_adm_noturb' (Sep 10th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
        
    % Wind farms for which SOWFA data is available
    case lower('2turb_alm_noturb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = 0.95;       % Turbine power scaling
        forcescale = 1.40;       % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = true;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice   
        lmu        = 1.7;        % Mixing length in x-direction (m)
        mu         = 0*18e-5;    % Dynamic flow viscosity        
        m          = 1;          % Turbulence model gridding property        
        n          = 4;          % Turbulence model gridding property

        % Tuning notes '2turb_alm_noturb' (Sep 5th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
                  
    case lower('2turb_alm_turb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);            % Load the LES meshing file
        Drotor     = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale = 0.95;       % Turbine power scaling
        forcescale = 1.4;        % Turbine force scaling
        p_init     = 0.0;        % Initial values for pressure terms (Pa)
        turbul     = true;       % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';   % Turbulence model of choice   
        lmu        = 1.4;        % Mixing length in x-direction (m)
        mu         = 0*18e-5;    % Dynamic flow viscosity
        m          = 1;          % Turbulence model gridding property        
        n          = 3;          % Turbulence model gridding property

        % Tuning notes '2turb_alm_turb' (Sep 6th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.0, m = 1:8, n = 1:4
         
    case lower('apc_9turb_alm_turb')
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);             % Load the LES meshing file
        Drotor      = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale  = 0.97;        % Turbine power scaling
        forcescale  = 2.0;        % Turbine force scaling
        p_init   = 0.0;           % Initial values for pressure terms (Pa)
        turbul   = true;          % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';    % Turbulence model of choice   
        lmu      = 1.20;          % Mixing length in x-direction (m)
        mu       = 0*18e-5;       % Dynamic flow viscosity        
        m        = 7;             % Turbulence model gridding property          
        n        = 1;             % Turbulence model gridding property
        
        % Tuning notes 'apc_9turb_alm_turb' (Sep 11th, 2017): 
        % Ranges: lmu= 0.1:0.1:2.0, f = 0.8:0.1:2.5, m = 1:8, n = 1:4
        
    case lower('yaw_2turb_alm_noturb')
        error('This case has not yet been tuned/validated. Remove this message manually in meshing.m to continue.');
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);             % Load the LES meshing file
        Drotor      = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale  = 0.95;        % Turbine power scaling
        forcescale  = 1.4;        % Turbine force scaling
        p_init   = 0.0;           % Initial values for pressure terms (Pa)
        turbul   = true;          % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';    % Turbulence model of choice   
        lmu      = 1.0;           % Mixing length in x-direction (m)
        mu       = 0*18e-5;       % Dynamic flow viscosity        
        m        = 2;             % Turbulence model gridding property        
        n        = 1;             % Turbulence model gridding property
                
    case lower('yaw_2turb_alm_turb')
        error('This case has not yet been tuned/validated. Remove this message manually in meshing.m to continue.');
        [meshFn,measurementFn] = downloadLESdata( WFSimfolder, lower(scenarioName) ); % Download files
        load(meshFn);             % Load the LES meshing file
        Drotor      = Drotor(1);  % WFSim only supports a uniform Drotor for now
        powerscale  = 1.0;        % Turbine power scaling
        forcescale  = 1.2;        % Turbine force scaling
        p_init   = 0.0;           % Initial values for pressure terms (Pa)
        turbul   = true;          % Use mixing length turbulence model (true/false)        
        turbModel  = 'WFSim3';    % Turbulence model of choice   
        lmu      = 1.0;           % Mixing length in x-direction (m)
        mu       = 0*18e-5;       % Dynamic flow viscosity        
        m        = 8;             % Turbulence model gridding property        
        n        = 2;             % Turbulence model gridding property
       
    otherwise
        error('No valid meshing specified. Please take a look at Wp.name.');
end;

% calculate time
time = 0:h:L;          % time vector for simulation
NN   = length(time)-1; % Total number of simulation steps

% compute input to linear model
for kk=2:length(time)
    turbInput(kk-1).dCT_prime = turbInput(kk).CT_prime - turbInput(kk-1).CT_prime; 
end
turbInput(end).dCT_prime = turbInput(end-1).dCT_prime;

% construct grid
ldx  = linspace(0,Lx,Nx);
ldy  = linspace(0,Ly,Ny);

if strcmp(lower(gridType),'lin')
    % linear gridding
    ldx  = linspace(0,Lx,Nx);
    ldy  = linspace(0,Ly,Ny);
elseif strcmp(lower(gridType),'exp')
    error('Exponential gridding currently not supported.');
%     % exponential gridding
%     ldx  = expvec( Lx, cellSize(1), uniquetol(Crx,1e-2), R_con(1), N_con(1), dx_min(1) );
%     ldy  = expvec( Ly, cellSize(2), uniquetol(Cry,1e-2), R_con(2), N_con(2), dx_min(2) );
%     Nx   = length(ldx);
%     Ny   = length(ldy);
%    
else
    error('Wrong meshing type specified in "type".');
end;

ldxx = repmat(ldx',1,Ny);
ldyy = repmat(ldy,Nx,1);

% Create secondary grid from primary grid
ldx2  = 0.5*(ldx(1:end-1)+ldx(2:end));
ldx2  = [ldx2 2*ldx2(end)-ldx2(end-1)]; % add extra cells
ldy2  = 0.5*(ldy(1:end-1)+ldy(2:end));
ldy2  = [ldy2 2*ldy2(end)-ldy2(end-1)]; % add extra cells
ldxx2 = repmat(ldx2',1,Ny);
ldyy2 = repmat(ldy2,Nx,1);

% Calculate cell dimensions
dx   = diff(ldx);
dxx  = repmat([dx'; dx(end)],1,Ny);
dx2  = diff(ldx2);
dxx2 = repmat([dx2'; dx2(end)],1,Ny);
dy   = diff(ldy);
dyy  = repmat([dy, dy(end)],Nx,1);
dy2  = diff(ldy2);
dyy2 = repmat([dy2, dy2(end)],Nx,1);

% Calculate location of turbines in grid and grid mismatch
for i = 1:length(Crx)
    % Calculate cells relevant for turbine (x-dir) on primary grid
    [~,xline(i,1)] = min(abs(ldx-Crx(i))); % automatically picks earliest entry in vector
    
    % Calculate cells closest to turbines (y-dir) on both grids
    [ML_prim, L_prim ] = min(abs(ldy- (Cry(i)-Drotor/2)));
    [ML_sec , L_sec  ] = min(abs(ldy2-(Cry(i)-Drotor/2)));
    [MR_prim, R_prim ] = min(abs(ldy- (Cry(i)+Drotor/2)));
    [MR_sec , R_sec  ] = min(abs(ldy2-(Cry(i)+Drotor/2)));
    
    yline{i}  = L_prim:1: R_prim; % turbine cells for primary grid
    % ylinev{i} = L_sec :1: R_sec ; % turbine cells for secondary grid
    ylinev{i} = L_prim:1: R_prim+1; % JWs code -- Bart: I feel like this needs fixing
    
    if PrintGridMismatch
        % Calculate turbine-grid mismatch
        disp([' TURBINE ' num2str(i) ' GRID MISMATCH:']);
        disp(['                    Primary           Secondary ']);
        disp(['       center:   (' num2str(min(abs(Crx(i)-ldx)),'%10.2f/n') ',' num2str(min(abs(Cry(i)-ldy)),'%10.2f/n') ') m.      (' num2str(min(abs(Crx(i)-ldx2)),'%10.2f/n') ',' num2str(min(abs(Cry(i)-ldy2)),'%10.2f/n') ') m.']);
        disp(['   left blade:   (' num2str(min(abs(Crx(i)-ldx)),'%10.2f/n') ',' num2str(ML_prim,'%10.2f/n')             ') m.      (' num2str(min(abs(Crx(i)-ldx2)),'%10.2f/n') ',' num2str(ML_sec,'%10.2f/n')                ') m.']);
        disp(['  right blade:   (' num2str(min(abs(Crx(i)-ldx)),'%10.2f/n') ',' num2str(MR_prim,'%10.2f/n')             ') m.      (' num2str(min(abs(Crx(i)-ldx2)),'%10.2f/n') ',' num2str(MR_sec,'%10.2f/n')                ') m.']);
        disp(' ');
    end;
end;

% Calculate state sizes
Nu = (Nx-3)*(Ny-2);   % Number of u velocities in state vector
Nv = (Nx-2)*(Ny-3);   % Number of v velocities in state vector
Np = (Nx-2)*(Ny-2)-2; % Number of pressure terms in state vector

%% Display results
if plotMesh
    clf;
    Z1 = -2*ones(size(ldxx));
    mesh(ldyy,ldxx,Z1,'FaceAlpha',0,'EdgeColor','black','LineWidth',0.1,'DisplayName','Primary mesh')
    hold on;
    Z2 = -1*ones(size(ldxx2));
    mesh(ldyy2,ldxx2,Z2,'FaceAlpha',0,'EdgeColor','blue','LineWidth',0.1,'DisplayName','Secondary mesh')
    hold on;
    for j = 1:length(Cry)
        plot([Cry(j)-Drotor/2;Cry(j)+Drotor/2],[Crx(j);Crx(j)],'LineWidth',3.0,'DisplayName',['Turbine ' num2str(j)])
        hold on
        text(Cry(j),Crx(j),['T ' num2str(j)])
    end;
    axis equal
    xlim([-0.1*Ly 1.2*Ly]);
    ylim([-0.1*Lx 1.1*Lx]);
    xlabel('y (m)');
    ylabel('x (m)');
    view(0,90); % view from the top
    legend('-dynamicLegend');
    drawnow;
end;

%% Export to Wp and input
Wp         = struct('Nu',Nu,'Nv',Nv,'Np',Np,'name',scenarioName);
Wp.sim     = struct('h',h,'time',time,'L',L,'NN',NN);
Wp.turbine = struct('Drotor',Drotor,'powerscale',powerscale,'forcescale',forcescale, ...
    'N',length(Crx),'Crx',Crx,'Cry',Cry);
Wp.turbine.input = turbInput; % System inputs too
Wp.site    = struct('mu',mu,'Rho',Rho,'u_Inf',u_Inf,'v_Inf',v_Inf,'p_init',p_init, ...
    'lmu',lmu,'turbul',turbul,'m',m,'n',n,'Turbulencemodel',turbModel);
Wp.mesh    = struct('Lx',Lx,'Ly',Ly,'Nx',Nx,'Ny',Ny,'ldxx',ldxx,'ldyy',ldyy,'ldxx2',...
    ldxx2,'ldyy2',ldyy2,'dxx',dxx,'dyy',dyy,'dxx2',dxx2,'dyy2',dyy2,...
    'xline',xline,'type',gridType);
Wp.mesh.yline = yline; Wp.mesh.ylinev = ylinev; % Do not support struct command

if exist('measurementFn') == 1
    Wp.sim.measurementFile = measurementFn;
end

%% Construct mu if no turbulence
if turbul==0; mu = ConstructMu(Wp); Wp.site(:).mu = mu; end

end