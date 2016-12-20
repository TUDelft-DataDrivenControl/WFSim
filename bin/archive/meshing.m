function Wp = meshing(plotting,Wp,MeshingMethod)
% ddx   = \Delta x_{I,I+1}
% ddy   = \Delta y_{J,J+1}
% ddx2  = \Delta x_{i,i+1}
% ddy2  = \Delta y_{j,j+1}
% dx    = \Delta x_{I-1,I}
% dy    = \Delta y_{J-1,J}
% dx2   = \Delta x_{i-1,i}
% dy2   = \Delta y_{j-1,j}
% ldxx  = I
% ldyy  = Jinput.beta(k,:)
% ldxx2 = i
% ldyy2 = j

switch (Wp.name)
    case {'apc_3x3turb_noyaw_9turb_100x50_lin'}
        Lx     = 2863.9514;
        Ly     = 2158.4055;
        Nx     = 100;
        Ny     = 50;
        Crx    = [400;1031.976;399.98;399.99;1031.956;1031.966;1663.931;1663.941;1663.951];
        Cry    = [700;700;1458.405;1079.203;1458.405;1079.203;1458.405;1079.203;700];
        turbine.Drotor = 126.3992;
        turbine.MaxCp  = 0.4866;
        Crxs   = [1;1;1;1;1;1;1;1;1];   % Obsolete option
        Qx     = 1;        % Irrelevant for linear grids
        Qy     = 1;        % Irrelevant for linear grids
        sigmax = 1;        % Irrelevant for linear grids
        sigmay = 1;        % Irrelevant for linear grids
   
    case {'apc_3x3turb_noyaw_9turb_50x25_lin'}
        Lx     = 2863.9514;
        Ly     = 2158.4055;
        Nx     = 50;
        Ny     = 25;
        Crx    = [400;1031.976;399.98;399.99;1031.956;1031.966;1663.931;1663.941;1663.951];
        Cry    = [700;700;1458.405;1079.203;1458.405;1079.203;1458.405;1079.203;700];
        turbine.Drotor = 126.3992;
        turbine.MaxCp  = 0.4866;
        Crxs   = [1;1;1;1;1;1;1;1;1];   % Obsolete option
        Qx     = 1;        % Irrelevant for linear grids
        Qy     = 1;        % Irrelevant for linear grids
        sigmax = 1;        % Irrelevant for linear grids
        sigmay = 1;        % Irrelevant for linear grids
        
    case {'noprecursor_2turb_60x30_lin'}
        Lx     = 2232.0623;
        Ly     = 1400;
        Nx     = 60;
        Ny     = 30;
        Crx    = [400; 1032.062];
        Cry    = [700; 700];
        turbine.Drotor = 126.3992;
        %turbine.MaxCp  = 0.4866;
        Crxs   = [1;1];    % Obsolete option
        Qx     = 1;        % Irrelevant for linear grids
        Qy     = 1;        % Irrelevant for linear grids
        sigmax = 1;        % Irrelevant for linear grids
        sigmay = 1;        % Irrelevant for linear grids
        
    case {'WithPrecursor_2turb_50x25_lin'}
        Lx     = 2232.0623;
        Ly     = 1400;
        Nx     = 50;
        Ny     = 25;
        Crx    = [400; 1032.062];
        Cry    = [700; 700];
        turbine.Drotor = 126.3992;
        turbine.MaxCp  = 0.4866;
        Crxs   = [1;1];    % Obsolete option
        Qx     = 1;        % Irrelevant for linear grids
        Qy     = 1;        % Irrelevant for linear grids
        sigmax = 1;        % Irrelevant for linear grids
        sigmay = 1;        % Irrelevant for linear grids
        
    case {'yaw_2turb_50x25_lin'}
        Lx     = 2481.9702;
        Ly     = 1400;
        Nx     = 50;
        Ny     = 25;
        Crx    = [400;1281.97];
        Cry    = [700;700];
        turbine.Drotor = 126.3992;
        turbine.MaxCp  = 0.4866;
        Crxs   = [1;1];   % Obsolete option
        Qx     = 1;        % Irrelevant for linear grids
        Qy     = 1;        % Irrelevant for linear grids
        sigmax = 1;        % Irrelevant for linear grids
        sigmay = 1;        % Irrelevant for linear grids
        Wp.powerscale = 1.1; 
        
    case {'yaw_2turb_100x50_lin'}
        Lx     = 2481.9702;
        Ly     = 1400;
        Nx     = 100;
        Ny     = 50;
        Crx    = [400;1281.97];
        Cry    = [700;700];
        turbine.Drotor = 126.3992;
        turbine.MaxCp  = 0.4866;
        Crxs   = [1;1];   % Obsolete option
        Qx     = 1;        % Irrelevant for linear grids
        Qy     = 1;        % Irrelevant for linear grids
        sigmax = 1;        % Irrelevant for linear grids
        sigmay = 1;        % Irrelevant for linear grids
        
    case {'sowfa_coarse_yaw'}
        turbine.MaxCp  = 0.4866;              % Maximup coefficient
        turbine.Drotor = 126.3992;            % Rotor diameter
        Lx          =   2000;                    % Length of the grid in x (N-S direction)
        Ly          =   1000;                    % Length of the grid in y (O-W direction)
        Cry         =   [Ly/2; Ly/2];            % X-coordinate of rotor (center)
        Crx         =   [400; 400+126.3992*7];         % Y-coordinate of rotor (center)
        Crxs        =   Crx-126.4;               % X-coordinate of rotor (center)
        Qx          =   1;                       % change of grid size (x-direction)
        Qy          =   1;                       % change of grid size (y-direction)
        Nx          =   50;                     % Number of grid points (x-direction)
        Ny          =   25;                      % Number of grid points (y-direction)
        sigmax      =   40;                      % Number of grid points (x-direction)
        sigmay      =   80;                      % Number of grid points (y-direction)

        %     case {'benchmark6_ieee'}
        %         load V90_data
        %         Lx      =   5000;            % Length of the grid in x (N-S direction)
        %         Ly      =   3000;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   90;             % Rotor diameter
        %         Cry     =   [1000; 1000+8*turbine.Drotor; 1000; 1000+8*turbine.Drotor; 1000; 1000+8*turbine.Drotor];
        %         Crx     =   [1800; 1800; 1800+10*turbine.Drotor; 1800+10*turbine.Drotor; 1800+20*turbine.Drotor; 1800+20*turbine.Drotor];             % Y-coordinate of rotor (center)
        %         Crxs     =   Crx-360;            % X-coordinate of rotor (center)
        %         Qx          =  1;% 40;   % change of grid size (x-direction)
        %         Qy          =   1;%p;20;   % change of grid size (y-direction)
        %         Nx          =   100;   % Number of grid points (x-direction)
        %         Ny          =   50;   % Number of grid points (y-direction)
        %         sigmax      =   40;   % Number of grid points (x-direction)
        %         sigmay      =   80;   % Number of grid points (y-direction)\
        %     case {'benchmark6_yaw'}
        %         load V90_data
        %         Lx      =   3800;            % Length of the grid in x (N-S direction)
        %         Ly      =   2700;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   90;             % Rotor diameter
        %         Cry     =   [1080; 1080+6*turbine.Drotor; 1080-0*turbine.Drotor; 1080+6*turbine.Drotor+0*turbine.Drotor; 1080; 1080+6*turbine.Drotor];
        %         Crx     =   [1260; 1260; 1260+7*turbine.Drotor; 1260+7*turbine.Drotor; 1260+14*turbine.Drotor; 1260+14*turbine.Drotor];             % Y-coordinate of rotor (center)
        %         Crxs     =   Crx-360;            % X-coordinate of rotor (center)
        %         Qx          =   1;%40;   % change of grid size (x-direction)
        %         Qy          =   1;%20;   % change of grid size (y-direction)
        %         Nx          =   150;   % Number of grid points (x-direction)
        %         Ny          =   150;   % Number of grid points (y-direction)
        %         sigmax      =   40;   % Number of grid points (x-direction)
        %         sigmay      =   80;   % Number of grid points (y-direction)
        %     case {'benchmark2_yaw'}
        %         load V90_data
        %         Lx      =   2000;            % Length of the grid in x (N-S direction)
        %         Ly      =   600;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   90;             % Rotor diameter
        %         Cry     =   [300;300-45];             % X-coordinate of rotor (center)
        %         Crx     =   [500; 500+6*turbine.Drotor];             % Y-coordinate of rotor (center)
        %         Crxs     =   Crx-360;%[55; 80];             % X-coordinate of rotor (center)
        %         Qx          =  1;% 40;   % change of grid size (x-direction)
        %         Qy          =   1;%p;20;   % change of grid size (y-direction)
        %         Nx          =   100;   % Number of grid points (x-direction)
        %         Ny          =   50;   % Number of grid points (y-direction)
        %         sigmax      =   40;   % Number of grid points (x-direction)
        %         sigmay      =   80;   % Number of grid points (y-direction)

        %     case {'benchmark2_ieee'}
        %         load V90_data
        %         Lx      =   6000;            % Length of the grid in x (N-S direction)
        %         Ly      =   2000;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   90;             % Rotor diameter
        %         Cry     =   [1000; 1000];             % X-coordinate of rotor (center)
        %         Crx     =   [1800; 3800];             % Y-coordinate of rotor (center)
        %         Crxs     =   Crx-360;%[55; 80];             % X-coordinate of rotor (center)
        %         Qx          =  1;% 40;   % change of grid size (x-direction)
        %         Qy          =   1;%p;20;   % change of grid size (y-direction)
        %         Nx          =   100;   % Number of grid points (x-direction)
        %         Ny          =   50;   % Number of grid points (y-direction)
        %         sigmax      =   40;   % Number of grid points (x-direction)
        %         sigmay      =   80;   % Number of grid points (y-direction)
        %     case {'benchmark2_coarse'}
        %         load V90_data
        %         Lx      =   7000;            % Length of the grid in x (N-S direction)
        %         Ly      =   2400;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   90;             % Rotor diameter
        %         Cry         =   [1200; 1200];             % X-coordinate of rotor (center)
        %         Crx         =   [3800; 4800];             % Y-coordinate of rotor (center)
        %         Crxs        =   Crx-360;%[55; 80];             % X-coordinate of rotor (center)
        %         Qx          =  40;   % change of grid size (x-direction)
        %         Qy          =  20;   % change of grid size (y-direction)
        %         Nx          =   60;   % Number of grid points (x-direction)
        %         Ny          =   30;   % Number of grid points (y-direction)
        %         sigmax      =   40;   % Number of grid points (x-direction)
        %         sigmay      =   80;   % Number of grid points (y-direction)
        %     case {'benchmark2_coarser'}
        %         load V90_data
        %         Lx      =   4000;            % Length of the grid in x (N-S direction)
        %         Ly      =   1500;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   90;             % Rotor diameter
        %         Cry         =   [750; 750];             % X-coordinate of rotor (center)
        %         Crx         =   [1500; 1500+7*turbine.Drotor];             % Y-coordinate of rotor (center)
        %         Crxs        =   Crx-360;%[55; 80];             % X-coordinate of rotor (center)
        %         Qx          =  40;   % change of grid size (x-direction)
        %         Qy          =  20;   % change of grid size (y-direction)
        %         Nx          =   20;   % Number of grid points (x-direction)
        %         Ny          =   10;   % Number of grid points (y-direction)
        %         sigmax      =   40;   % Number of grid points (x-direction)
        %         sigmay      =   80;   % Number of grid points (y-direction)
        %     case {'benchmark1_ieee'}
        %         load V90_data
        %         Lx          =   4400;            % Length of the grid in x (N-S direction)
        %         Ly          =   2200;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   90;             % Rotor diameter
        %         Crx         =   1800;             % X-coordinate of rotor (center)
        %         Crxs        =   Crx-0.2*turbine.Drotor;%[55; 80];             % X-coordinate of rotor (center)
        %         Cry         =   1100;             % Y-coordinate of rotor (center)
        %         Qx          =   20;   % change of grid size (x-direction)
        %         Qy          =   40;   % change of grid size (y-direction)
        %         Nx          =   120;   % Number of grid points (x-direction)
        %         Ny          =   100;   % Number of grid points (y-direction)
        %         sigmax      =   350;   % Number of grid points (x-direction)
        %         sigmay      =   400;   % Number of grid points (y-direction)
        %     case {'benchmark1_close'}
        %         load V90_data
        %         Lx          =   1500;            % Length of the grid in x (N-S direction)
        %         Ly          =   500;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   90;             % Rotor diameter
        %         Crx         =   400;             % X-coordinate of rotor (center)
        %         Crxs        =   Crx-0.2*turbine.Drotor;%[55; 80];             % X-coordinate of rotor (center)
        %         Cry         =   250;             % Y-coordinate of rotor (center)
        %         Qx          =   20;   % change of grid size (x-direction)
        %         Qy          =   40;   % change of grid size (y-direction)
        %         Nx          =   120;   % Number of grid points (x-direction)
        %         Ny          =   100;   % Number of grid points (y-direction)
        %         sigmax      =   350;   % Number of grid points (x-direction)
        %         sigmay      =   400;   % Number of grid points (y-direction)
        %     case {'benchmark1_coarse'}
        %         load V90_data
        %         Lx      =   1500;            % Length of the grid in x (N-S direction)
        %         Ly      =   500;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   90;             % Rotor diameter
        %         Crx         =   400;             % X-coordinate of rotor (center)
        %         Crxs        =   Crx-0.2*turbine.Drotor;%[55; 80];             % X-coordinate of rotor (center)
        %         Cry         =   250;             % Y-coordinate of rotor (center)
        %         Qx          =  20;   % change of grid size (x-direction)
        %         Qy          =  40;   % change of grid size (y-direction)
        %         Nx          =   120;   % Number of grid points (x-direction)
        %         Ny          =   100;   % Number of grid points (y-direction)
        %         sigmax      =   350;   % Number of grid points (x-direction)
        %         sigmay      =   400;   % Number of grid points (y-direction)
        %     case {'benchmark1_minimal'}
        %         load V90_data
        %         Lx          =   5;            % Length of the grid in x (N-S direction)
        %         Ly          =   6;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   1;        % Rotor diameter
        %         Cry         =   3;            % X-coordinate of rotor (center) ?
        %         Crx         =   3;            % Y-coordinate of rotor (center) ?
        %         Crxs        =   Crx-0.4*turbine.Drotor;%[55; 80]; % X-coordinate of rotor (center)
        %         Qx          =   1;   % change of grid size (x-direction)
        %         Qy          =   1;   % change of grid size (y-direction)
        %         Nx          =   5;   % Number of grid points (x-direction)
        %         Ny          =   6;   % Number of grid points (y-direction)
        %         sigmax      =   5;   % Number of grid points (x-direction) exp grid
        %         sigmay      =   6;   % Number of grid points (y-direction) exp grid
        %         case {'benchmark1_minimal2'}
        %         load V90_data
        %         Lx          =   5;            % Length of the grid in x (N-S direction)
        %         Ly          =   5;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   1;        % Rotor diameter
        %         Cry         =   4;            % X-coordinate of rotor (center) ?
        %         Crx         =   3;            % Y-coordinate of rotor (center) ?
        %         Crxs        =   Crx-0.4*turbine.Drotor;%[55; 80]; % X-coordinate of rotor (center)
        %         Qx          =   1;   % change of grid size (x-direction)
        %         Qy          =   1;   % change of grid size (y-direction)
        %         Nx          =   5;   % Number of grid points (x-direction)
        %         Ny          =   5;   % Number of grid points (y-direction)
        %         sigmax      =   5;   % Number of grid points (x-direction) exp grid
        %         sigmay      =   6;   % Number of grid points (y-direction) exp grid
        %         turbine.Drotor=1;
        %
        %                 case {'benchmark1_minimal3'}
        %         load V90_data
        %         Lx          =   5;            % Length of the grid in x (N-S direction)
        %         Ly          =   5;            % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   1;        % Rotor diameter
        %         Cry         =   4;            % X-coordinate of rotor (center) ?
        %         Crx         =   3;            % Y-coordinate of rotor (center) ?
        %         Crxs        =   Crx-0.4*turbine.Drotor;%[55; 80]; % X-coordinate of rotor (center)
        %         Qx          =   1;   % change of grid size (x-direction)
        %         Qy          =   1;   % change of grid size (y-direction)
        %         Nx          =   6;   % Number of grid points (x-direction)
        %         Ny          =   5;   % Number of grid points (y-direction)
        %         sigmax      =   5;   % Number of grid points (x-direction) exp grid
        %         sigmay      =   6;   % Number of grid points (y-direction) exp grid
        %         turbine.Drotor=1;
        %     case {'benchmark1_simple'}
        %         load V90_data
        %         ss          =   20;
        %         s           =   1;
        %         Lx          =   ss;                 % Length of the grid in x (N-S direction)
        %         Ly          =   ss;                 % Length of the grid in y (O-W direction)
        %         turbine.Drotor  =   2;              % Rotor diameter
        %         Cry         =   ss/2;               % X-coordinate of rotor (center) ?
        %         Crx         =   ss/2;               % Y-coordinate of rotor (center) ?
        %         Crxs        =   Crx-0.4*turbine.Drotor;%[55; 80]; % X-coordinate of rotor (center)
        %         Qx          =   s;      % change of grid size (x-direction)
        %         Qy          =   s;      % change of grid size (y-direction)
        %         Nx          =   ss;     % Number of grid points (x-direction)
        %         Ny          =   ss;     % Number of grid points (y-direction)
        %         sigmax      =   ss;     % Number of grid points (x-direction) exp grid
        %         sigmay      =   ss;     % Number of grid points (y-direction) exp grid
        
    otherwise
        disp('Wind Farm not in list')
end
N=size(Crx,1);

%% For X
% initial grid
Xgrid=0:Lx/5000:Lx;
Xres=zeros(1,length(Xgrid));

for j=1:1:N
    Xres=exp(-(Xgrid-Crx(j)).^2/2/sigmax^2)+Xres;
end
switch lower(MeshingMethod)
    case {'lin'}
        Xres=(sign(Xres+1e-12)+1)./2;
    otherwise
        disp('exponential meshing')
end

XRES=((1+Xres*Qx));

q1 = 0;
for i=1:1:length(XRES)
    q1(i+1)=q1(i)+XRES(i);
end

%plot(q1./q1(end)*Lx,(0:1:length(XRES))/length((XRES))*Lx)
Dist_dxx = interp1(q1./q1(end)*Lx,(0:1:length(XRES))/length(XRES).*Lx,0:Lx/(Nx-1):Lx);
%Dist_dxx = linspace(0,Lx + Lx/(Nx-1), Nx+1);
dvdxx    = diff(Dist_dxx);
dxx      = [dvdxx dvdxx(end)];

dxx2(1)=dxx(1)/2;
for i=1:1:length(dxx)
    if i==length(dxx); else
        dxx2(i+1)=dxx(i)/2+dxx(i+1)/2;
    end
end
%% For Y
% initial grid
Ygrid=0:Ly/5000:Ly;
Yres=zeros(1,length(Ygrid));

for j=1:1:N
    Yres=exp(-(Ygrid-Cry(j)).^2/2/sigmay^2)+Yres;
end
switch lower(MeshingMethod)
    case {'lin'}
        Yres=(sign(Yres+1e-12)+1)./2;
    otherwise
        disp(' ')
end
YRES=((1+Yres*Qy));

q1=[0];
for i=1:1:length(YRES)
    q1(i+1)=q1(i)+YRES(i);
end

%plot(q1,(0:1:length(XRES))/length((XRES)*Lx)
Dist_dyy = interp1(q1./q1(end)*Ly,(0:1:length(YRES))/length(YRES).*Ly,0:Ly/(Ny-1):Ly);
%Dist_dyy = linspace(0,Ly + Ly/(Ny-1), Ny+1);
dvdyy    = diff(Dist_dyy);
dyy      = [dvdyy dvdyy(end)];

dyy2(1)=dyy(1)/2;
for i=1:1:length(dyy)
    if i==length(dyy); else
        dyy2(i+1)=dyy(i)/2+dyy(i+1)/2;
    end
end

%%
dxx=dxx'*ones(1,Ny);
dxx2=dxx2'*ones(1,Ny);
dyy=ones(Nx,1)*dyy;
dyy2=ones(Nx,1)*dyy2;

% ldxx=Dist_dxx(1:end-1)'*ones(1,Ny);
% ldyy=ones(Nx,1)*Dist_dyy(1:end-1);
%
% Dist_dxx2=Dist_dxx(1:end-1)+diff(Dist_dxx)-0.5*Dist_dxx(1,2);
% Dist_dyy2=Dist_dyy(1:end-1)+diff(Dist_dyy)-0.5*Dist_dyy(1,2);

ldxx=Dist_dxx'*ones(1,Ny);
ldyy=ones(Nx,1)*Dist_dyy;

dvdxx = diff(Dist_dxx);
dvdyy = diff(Dist_dyy);

Dist_dxx2=Dist_dxx+[dvdxx dvdxx(end)]-0.5*Dist_dxx(1,2);
Dist_dyy2=Dist_dyy+[dvdyy dvdyy(end)]-0.5*Dist_dyy(1,2);

ldxx2=Dist_dxx2(1:end)'*ones(1,Ny);
ldyy2=ones(Nx,1)*Dist_dyy2(1:end);

if plotting==1;
    figure(20)
    axis([0 Ly 0 Lx])
    vline(Dist_dyy,'k')
    hline(Dist_dxx,'k')
    xlabel('y');ylabel('x');
end

for i=1:1:N
    [MINNIE,xline(i,:)]   =   min((Crx(i)-ldxx2(:,1)).^2);%floor(Crx/Dx);  % Y grid number of the turbine
end
for i=1:1:N
    [MINNIE,xlines(i,:)]   =   min((Crxs(i)-(ldxx2(:,1))).^2);%floor(Crx/Dx);  % Y grid number of the turbine
end

for i=1:1:N
    [Minnie,Q1]=min(((Cry(i)-0.5*turbine.Drotor-ldyy(1,:)')).^2);
    [Minnie,Q2]=min(((Cry(i)+0.5*turbine.Drotor-ldyy(1,:)')).^2);
    
    yline{i}   =   floor(Q1):1:floor(Q2);
    ylinev{i}   =   floor(Q1):1:floor(Q2)+1; %added JW
    %  ylineu{i}   =   floor(Q1)-1:1:floor(Q2)+1; %added JW
end

Wp.dxx     = dxx;
Wp.dyy     = dyy;
Wp.dxx2    = dxx2;
Wp.dyy2    = dyy2;
Wp.ldxx    = ldxx;
Wp.ldyy    = ldyy;
Wp.ldxx2   = ldxx2;
Wp.ldyy2   = ldyy2;
Wp.Nx      = Nx;
Wp.Ny      = Ny;
Wp.Lx      = Lx;
Wp.Ly      = Ly;
Wp.xline   = xline;
Wp.yline   = yline;
Wp.ylinev   = ylinev;
%  Wp.ylineu   = ylineu;
Wp.xlines  = xlines;
Wp.N       = N;
Wp.turbine = turbine;
Wp.Crx     = Crx;
Wp.Cry     = Cry;
Wp.Crxs    = Crxs;

% meshgrid(Wp.ldyy(1,1:end),Wp.ldxx(1:end,1))

%% some modifications to match sowfa

switch lower(Wp.name)
    case {'sowfa_coarse_yaw'}
        Wp.yline{1}=Wp.yline{1}(2:end-1);
        Wp.yline{2}=Wp.yline{2}(2:end-1);
        Wp.ylinev{1}=Wp.ylinev{1}(2:end-1);
        Wp.ylinev{2}=Wp.ylinev{2}(2:end-1);
        Wp.powerscale=turbine.Drotor/(ldyy(1,Wp.yline{1}(end))-ldyy(1,Wp.yline{1}(1)));
    case {'sowfa_coarse'}
        Wp.yline{1}=Wp.yline{1}(2:end-1);
        Wp.yline{2}=Wp.yline{2}(2:end-1);
        Wp.ylinev{1}=Wp.ylinev{1}(2:end-1);
        Wp.ylinev{2}=Wp.ylinev{2}(2:end-1);
        Wp.powerscale=turbine.Drotor/(ldyy(1,Wp.yline{1}(end))-ldyy(1,Wp.yline{1}(1)));
    case {'sowfa_yaw'}
        Wp.yline{1}=Wp.yline{1}(2:end-1);
        Wp.yline{2}=Wp.yline{2}(2:end-1);
        Wp.ylinev{1}=Wp.ylinev{1}(2:end-1);
        Wp.ylinev{2}=Wp.ylinev{2}(2:end-1);
        Wp.powerscale=turbine.Drotor/(ldyy(1,Wp.yline{1}(end))-ldyy(1,Wp.yline{1}(1)))*1.2061;
    case {'sowfa'}
        Wp.yline{1}=Wp.yline{1}(2:end-1);
        Wp.yline{2}=Wp.yline{2}(2:end-1);
        Wp.ylinev{1}=Wp.ylinev{1}(2:end-1);
        Wp.ylinev{2}=Wp.ylinev{2}(2:end-1);
        Wp.powerscale=turbine.Drotor/(ldyy(1,Wp.yline{1}(end))-ldyy(1,Wp.yline{1}(1)))*1.2061;
        
    case {'sowfa_power'}
        Wp.yline{1}=Wp.yline{1}(2:end);
        Wp.yline{4}=Wp.yline{4}(2:end);
        Wp.yline{7}=Wp.yline{7}(2:end);
        
        Wp.yline{3}=Wp.yline{3}(1:end-1);
        Wp.yline{6}=Wp.yline{6}(1:end-1);
        Wp.yline{9}=Wp.yline{9}(1:end-1);
        
        Wp.ylinev{1}=Wp.ylinev{1}(2:end);
        Wp.ylinev{4}=Wp.ylinev{4}(2:end);
        Wp.ylinev{7}=Wp.ylinev{7}(2:end);
        
        Wp.ylinev{3}=Wp.ylinev{3}(1:end-1);
        Wp.ylinev{6}=Wp.ylinev{6}(1:end-1);
        Wp.ylinev{9}=Wp.ylinev{9}(1:end-1);
        
        Wp.powerscale=turbine.Drotor/(ldyy(1,Wp.yline{1}(end))-ldyy(1,Wp.yline{1}(1)))*1.2061;
        
    otherwise
        Wp.powerscale=turbine.Drotor/(ldyy(1,Wp.yline{1}(end))-ldyy(1,Wp.yline{1}(1)));
end