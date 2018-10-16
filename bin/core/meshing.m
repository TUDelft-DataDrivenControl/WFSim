function [ Wp ] = meshing( Wp, plotMesh, PrintGridMismatch )
%MESHING Meshing and settings function for the WFSim code
% This code calculates the meshing and prepares the Wp struct of a
% specific wind farm case for WFSim simulations

% Default settings
if nargin <= 0; error('Please specify a meshing case.'); end;
if nargin <= 1; plotMesh = true;                         end;
if nargin <= 2; PrintGridMismatch = true;                end;

% construct grid
ldx  = linspace(0,Wp.mesh.Lx,Wp.mesh.Nx);
ldy  = linspace(0,Wp.mesh.Ly,Wp.mesh.Ny);

% linear gridding
ldx  = linspace(0,Wp.mesh.Lx,Wp.mesh.Nx);
ldy  = linspace(0,Wp.mesh.Ly,Wp.mesh.Ny);

ldxx = repmat(ldx',1,Wp.mesh.Ny);
ldyy = repmat(ldy,Wp.mesh.Nx,1);

% Create secondary grid from primary grid
ldx2  = 0.5*(ldx(1:end-1)+ldx(2:end));
ldx2  = [ldx2 2*ldx2(end)-ldx2(end-1)]; % add extra cells
ldy2  = 0.5*(ldy(1:end-1)+ldy(2:end));
ldy2  = [ldy2 2*ldy2(end)-ldy2(end-1)]; % add extra cells
ldxx2 = repmat(ldx2',1,Wp.mesh.Ny);
ldyy2 = repmat(ldy2,Wp.mesh.Nx,1);

% Calculate cell dimensions
dx   = diff(ldx);
dxx  = repmat([dx'; dx(end)],1,Wp.mesh.Ny);
dx2  = diff(ldx2);
dxx2 = repmat([dx2'; dx2(end)],1,Wp.mesh.Ny);
dy   = diff(ldy);
dyy  = repmat([dy, dy(end)],Wp.mesh.Nx,1);
dy2  = diff(ldy2);
dyy2 = repmat([dy2, dy2(end)],Wp.mesh.Nx,1);

% Calculate location of turbines in grid and grid mismatch
for i = 1:length(Wp.turbine.Crx)
    % Calculate cells relevant for turbine (x-dir) on primary grid
    [~,xline(i,1)] = min(abs(ldx-Wp.turbine.Crx(i))); % automatically picks earliest entry in vector
    
    % Calculate cells closest to turbines (y-dir) on both grids
    [ML_prim, L_prim ] = min(abs(ldy- (Wp.turbine.Cry(i)-Wp.turbine.Drotor/2)));
    [ML_sec , L_sec  ] = min(abs(ldy2-(Wp.turbine.Cry(i)-Wp.turbine.Drotor/2)));
    [MR_prim, R_prim ] = min(abs(ldy- (Wp.turbine.Cry(i)+Wp.turbine.Drotor/2)));
    [MR_sec , R_sec  ] = min(abs(ldy2-(Wp.turbine.Cry(i)+Wp.turbine.Drotor/2)));
    
    yline{i}  = L_prim:1: R_prim; % turbine cells for primary grid
    % ylinev{i} = L_sec :1: R_sec ; % turbine cells for secondary grid
    ylinev{i} = L_prim:1: R_prim+1; % JWs code -- Bart: I feel like this needs fixing
    
    if PrintGridMismatch
        % Calculate turbine-grid mismatch
        disp([' TURBINE ' num2str(i) ' GRID MISMATCH:']);
        disp(['                    Primary           Secondary ']);
        disp(['       center:   (' num2str(min(abs(Wp.turbine.Crx(i)-ldx)),'%10.2f/n') ',' num2str(min(abs(Wp.turbine.Cry(i)-ldy)),'%10.2f/n') ') m.      (' num2str(min(abs(Wp.turbine.Crx(i)-ldx2)),'%10.2f/n') ',' num2str(min(abs(Wp.turbine.Cry(i)-ldy2)),'%10.2f/n') ') m.']);
        disp(['   left blade:   (' num2str(min(abs(Wp.turbine.Crx(i)-ldx)),'%10.2f/n') ',' num2str(ML_prim,'%10.2f/n')             ') m.      (' num2str(min(abs(Wp.turbine.Crx(i)-ldx2)),'%10.2f/n') ',' num2str(ML_sec,'%10.2f/n')                ') m.']);
        disp(['  right blade:   (' num2str(min(abs(Wp.turbine.Crx(i)-ldx)),'%10.2f/n') ',' num2str(MR_prim,'%10.2f/n')             ') m.      (' num2str(min(abs(Wp.turbine.Crx(i)-ldx2)),'%10.2f/n') ',' num2str(MR_sec,'%10.2f/n')                ') m.']);
        disp(' ');
    end;
end;


%% Display results
if plotMesh
    clf;
    Z1 = -2*ones(size(ldxx));
    mesh(ldyy,ldxx,Z1,'FaceAlpha',0,'EdgeColor','black','LineWidth',0.1,'DisplayName','Primary mesh')
    hold on;
    Z2 = -1*ones(size(ldxx2));
    mesh(ldyy2,ldxx2,Z2,'FaceAlpha',0,'EdgeColor','blue','LineWidth',0.1,'DisplayName','Secondary mesh')
    hold on;
    for j = 1:length(Wp.turbine.Cry)
        plot([Wp.turbine.Cry(j)-Wp.turbine.Drotor/2;
              Wp.turbine.Cry(j)+Wp.turbine.Drotor/2],...
              [Wp.turbine.Crx(j);Wp.turbine.Crx(j)],...
              'LineWidth',3.0,'DisplayName',['Turbine ' num2str(j)])
        hold on
        text(Wp.turbine.Cry(j),Wp.turbine.Crx(j),['T ' num2str(j)])
    end;
    axis equal
    xlim([-0.1*Wp.mesh.Ly 1.2*Wp.mesh.Ly]);
    ylim([-0.1*Wp.mesh.Lx 1.1*Wp.mesh.Lx]);
    xlabel('y (m)');
    ylabel('x (m)');
    view(0,90); % view from the top
    legend('-dynamicLegend');
    drawnow;
end;

%% Export to Wp and input
Wp.Nu = (Wp.mesh.Nx-3)*(Wp.mesh.Ny-2);   % Number of u velocities in state vector
Wp.Nv = (Wp.mesh.Nx-2)*(Wp.mesh.Ny-3);   % Number of v velocities in state vector
Wp.Np = (Wp.mesh.Nx-2)*(Wp.mesh.Ny-2)-2; % Number of pressure terms in state vector

Wp.turbine.N = length(Wp.turbine.Crx);

% Write meshing
Wp.mesh.ldxx = ldxx;
Wp.mesh.ldyy = ldyy;
Wp.mesh.ldxx2= ldxx2;
Wp.mesh.ldyy2= ldyy2;
Wp.mesh.dxx = dxx;
Wp.mesh.dyy = dyy;
Wp.mesh.dxx2 = dxx2;
Wp.mesh.dyy2 = dyy2;
Wp.mesh.xline = xline;
Wp.mesh.yline = yline;
Wp.mesh.ylinev = ylinev;


%% Construct mu if no turbulence
if Wp.site.turbul==0; mu = ConstructMu(Wp); Wp.site(:).mu = mu; end

end