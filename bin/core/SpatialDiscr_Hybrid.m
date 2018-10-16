function [output] = SpatialDiscr_Hybrid(Wp,sol,Linearversion)
% dxx   = \Delta x_{I,I+1}
% dyy   = \Delta y_{J,J+1}
% dxx2  = \Delta x_{i,i+1}
% dyy2  = \Delta y_{j,j+1}
% ldxx  = I
% ldyy  = J
% ldxx2 = i
% ldyy2 = j

Nx     = Wp.mesh.Nx;
Ny     = Wp.mesh.Ny;
dxx    = Wp.mesh.dxx;
dyy    = Wp.mesh.dyy;
dxx2   = Wp.mesh.dxx2;
dyy2   = Wp.mesh.dyy2;

Rho    = Wp.site.Rho;

u      = sol.u;
v      = sol.v;

% Init
[ax.aE,ax.aW,ax.aS,ax.aN,ax.aP]         = deal(zeros(Nx,Ny));
[Fex,Fwx,Fsx,Fnx,dFex,dFwx,dFnx,dFsx]   = deal(zeros(Nx,Ny));
[ay.aE,ax.aW,ay.aS,ay.aN,ay.aP]         = deal(zeros(Nx,Ny));
[Fey,Fwy,Fsy,Fny,dFey,dFwy,dFny,dFsy]   = deal(zeros(Nx,Ny));

%%  Setting the coefficients according to the hybrid differencing scheme %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% x-direction

% Define convection coefficients and its derivatives
% Fex = c ( u_{i,J} + u_{i+1,J} )
dFex(1:Nx-1,1:Ny) = Rho*0.5*dyy2(1:Nx-1,1:Ny);
Fex(1:Nx-1,1:Ny)  = dFex(1:Nx-1,1:Ny).*( u(2:Nx,1:Ny) + u(1:Nx-1,1:Ny) );
% Few = c ( u_{i,J} + u_{i-1,J} )
dFwx(2:Nx,1:Ny)   = Rho*0.5*dyy2(2:Nx,1:Ny);
Fwx(2:Nx,1:Ny)    = dFwx(2:Nx,1:Ny).*( u(2:Nx,1:Ny) + u(1:Nx-1,1:Ny) ); 
% Fnx = c ( v_{I-1,j+1} + v_{I,j+1} )
dFnx(2:Nx,1:Ny-1) = Rho*0.5*dxx(2:Nx,1:Ny-1);
Fnx(2:Nx,1:Ny-1)  = dFnx(2:Nx,1:Ny-1).*( v(2:Nx,2:Ny) + v(1:Nx-1,2:Ny) );
% Fsx = c ( v_{I-1,j} + v_{I,j} )
dFsx(2:Nx,1:Ny)   = Rho*0.5*dxx(2:Nx,1:Ny);
Fsx(2:Nx,1:Ny)    = dFsx(2:Nx,1:Ny).*( v(2:Nx,1:Ny) + v(1:Nx-1,1:Ny));

ax.aE             = max(max(-Fex,-0.5*Fex),zeros(Nx,Ny));
ax.aW             = max(max(Fwx,0.5.*Fwx),zeros(Nx,Ny));
ax.aN             = max(max(-Fnx,-0.5*Fnx),zeros(Nx,Ny));
ax.aS             = max(max(Fsx,0.5*Fsx),zeros(Nx,Ny));
ax.aP             = ax.aW + ax.aE + ax.aS + ax.aN + Fex - Fwx + Fnx - Fsx;

if Linearversion
    [dax.E,dax.W,dax.N,dax.S,dax.P]         = deal(zeros(Nx,Ny));
    [dax.SW,dax.NW,dax.SE,dax.NE]           = deal(zeros(Nx,Ny));
    [day.E,day.W,day.N,day.S,day.P]         = deal(zeros(Nx,Ny));
    [day.SW,day.NW,day.SE,day.NE]           = deal(zeros(Nx,Ny));
    
    % daxe/du_(i,J) = daxe/du_(i+1,J)
    dax.aE            = (-Fex>=(-Fex/2)).*(-Fex>zeros(Nx,Ny)).*-dFex + ((-Fex/2)>-Fex).*((-Fex/2)>zeros(Nx,Ny)).*-dFex/2; 
    % daxw/du_(i,J) = daxe/du_(i-1,J)
    dax.aW            = (Fwx>=(Fwx/2)).*(Fwx>zeros(Nx,Ny)).*dFwx + ((Fwx/2)>Fwx).*((Fwx/2)>zeros(Nx,Ny)).*dFwx/2;
    % daxn/dv_(I,j+1) = daxn/dv_(I-1,j+1)
    dax.aN            = (-Fnx>=(-Fnx/2)).*(-Fnx>zeros(Nx,Ny)).*-dFnx + ((-Fnx/2)>-Fnx).*((-Fnx/2)>zeros(Nx,Ny)).*-dFnx/2;
    % daxs/dv_(I,j) = daxs/dv_(I-1,j)
    dax.aS            = (Fsx>=(Fsx/2)).*(Fsx>zeros(Nx,Ny)).*dFsx + ((Fsx/2)>Fsx).*((Fsx/2)>zeros(Nx,Ny)).*dFsx/2;
    
    % daPx/du_{i+1,J}
    dax.aPE           = dax.aE + dFex;
    % daPx/du_{i-1,J}
    dax.aPW           = dax.aW - dFwx;
    % daPx/dv_{I,j+1}
    dax.aPN           = dax.aN + dFnx;
    % daPx/dv_{I,j}             
    dax.aPS           = dax.aS - dFsx;
    % daPx/du_{i,J}
    dax.aPP           = dax.aW + dax.aE - dFwx + dFex; 
    
    % Define derivatives for linearized model with ax = -aPx u_{i,J} + aEx u_{i+1,J} + aWx u_{i-1,J} + aNx u_{i,J+1} + aSx u_{i,J-1}
    
    % dax/du_(i-1,J) = aWx + daWx/du_{i-1,J} u_{i-1,J} - daPx/du_{i-1,J} u_{i,J}
    dax.W(2:Nx,1:Ny)    = ax.aW(2:Nx,1:Ny) + dax.aW(2:Nx,1:Ny).*u(1:Nx-1,1:Ny) - dax.aPW(2:Nx,1:Ny).*u(2:Nx,1:Ny);
    
    % dax/du_(i,J-1) = aSx
    dax.S               = ax.aS;
    
    % dax/du_(i,J)   = -aPx + daWx/du_{i,J} u_{i-1,J} + daEx/du_{i,J} u_{i+1,J} - daPx/du_{i,J} u_{i,J}
    dax.P(2:Nx-1,1:Ny)  = ax.aP(2:Nx-1,1:Ny) - dax.aW(2:Nx-1,1:Ny).*u(1:Nx-2,1:Ny)-dax.aE(2:Nx-1,1:Ny).*u(3:Nx,1:Ny)...
        + dax.aPP(2:Nx-1,1:Ny).*u(2:Nx-1,1:Ny); 
    
    % dax/du_(i,J+1)   = aNx
    dax.N               = ax.aN;
    
    % dax/du_(i+1,J) = daEx/du_{i+1,J} u_{i+1,J} + aEx - daPx/du_{i+1,J} u_{i,J}
    dax.E(1:Nx-1,1:Ny)  = ax.aE(1:Nx-1,1:Ny)  + dax.aE(1:Nx-1,1:Ny).*u(2:Nx,1:Ny)  - dax.aPE(1:Nx-1,1:Ny).*u(1:Nx-1,1:Ny);
    
    % dax/dv_(I-1,j) = daSx/dv_{I-1,J} u_{i,J-1} - daPx/dv_{I-1,J} u_{i,J}
    dax.SW(2:Nx,2:Ny)   = dax.aS(2:Nx,2:Ny).*u(2:Nx,1:Ny-1) - dax.aPS(2:Nx,2:Ny).*u(2:Nx,2:Ny); 
    
    % dax/dv_(I-1,j+1) = daNx/dv_{I-1,j+1} u_{i,J+1} - daPx/dv_{I-1,j+1} u_{i,J}
    dax.NW(1:Nx-1,1:Ny-1) = dax.aN(1:Nx-1,1:Ny-1).*u(1:Nx-1,2:Ny)  - dax.aPN(1:Nx-1,1:Ny-1).*u(1:Nx-1,1:Ny-1);
    
    % dax/dv_(I,j)   = daSx/dv_{I,j} u_{I,j-1} - daPx/dv_{I,j} u_{i,J}
    dax.SE(1:Nx,2:Ny)   = dax.aS(1:Nx,2:Ny).*u(1:Nx,1:Ny-1)  - dax.aPS(1:Nx,2:Ny).*u(1:Nx,2:Ny);
    
    % dax/dv_(I,j+1) = daNx/dv_{I,j+1} u_{i,J+1} - daPx/dv_{I,j+1} u_{i,J}
    dax.NE(1:Nx,1:Ny-1) = dax.aN(1:Nx,1:Ny-1).*u(1:Nx,2:Ny)  - dax.aPN(1:Nx,1:Ny-1).*u(1:Nx,1:Ny-1);
end;

%% y-direction

% Define convection coefficients and its derivatives

% Fey = c ( u_{i+1,J} + u_{i+1,J-1} )
dFey(1:Nx-1,2:Ny) = Rho*0.5*dyy(1:Nx-1,2:Ny);
Fey(1:Nx-1,2:Ny)  = dFey(1:Nx-1,2:Ny).*( u(2:Nx,2:Ny) + u(2:Nx,1:Ny-1) );
% Fwy = c ( u_{i,J} + u_{i,J-1} )
dFwy(1:Nx,2:Ny)   = Rho*0.5*dyy(1:Nx,2:Ny);
Fwy(1:Nx,2:Ny)    = dFwy(1:Nx,2:Ny).*( u(1:Nx,2:Ny) + u(1:Nx,1:Ny-1) );
% Fny = c ( v_{I,j+1} + v_{I,j} )
dFny(1:Nx,1:Ny-1) = Rho*0.5*dxx2(1:Nx,1:Ny-1);
Fny(1:Nx,1:Ny-1)  = dFny(1:Nx,1:Ny-1).*( v(1:Nx,1:Ny-1) + v(1:Nx,2:Ny) );
% Fsy = c ( v_{I,j-1} + v_{I,j} )
dFsy(1:Nx,2:Ny)   = Rho*0.5*dxx2(1:Nx,2:Ny);
Fsy(1:Nx,2:Ny)    = dFsy(1:Nx,2:Ny).*( v(1:Nx,1:Ny-1) + v(1:Nx,2:Ny) );

ay.aE             = max(max(-Fey,-0.5*Fey),zeros(Nx,Ny));
ay.aW             = max(max(Fwy,0.5.*Fwy),zeros(Nx,Ny));
ay.aN             = max(max(-Fny,-0.5*Fny),zeros(Nx,Ny));
ay.aS             = max(max(Fsy,0.5*Fsy),zeros(Nx,Ny));
ay.aP             = ay.aW + ay.aE + ay.aS + ay.aN + Fey - Fwy + Fny - Fsy;

if Linearversion
    day.aE            = (-Fey>=(-Fey/2)).*(-Fey>zeros(Nx,Ny)).*-dFey + ((-Fey/2)>-Fey).*((-Fey/2)>zeros(Nx,Ny)).*-dFey/2;
    day.aW            = (Fwy>=(Fwy/2)).*(Fwy>zeros(Nx,Ny)).*dFwy + ((Fwy/2)>Fwy).*((Fwy/2)>zeros(Nx,Ny)).*dFwy/2;
    day.aN            = (-Fny>=(-Fny/2)).*(-Fny>zeros(Nx,Ny)).*-dFny + ((-Fny/2)>-Fny).*((-Fny/2)>zeros(Nx,Ny)).*-dFny/2;
    day.aS            = (Fsy>=(Fsy/2)).*(Fsy>zeros(Nx,Ny)).*dFsy + ((Fsy/2)>Fsy).*((Fsy/2)>zeros(Nx,Ny)).*dFsy/2;
    
    % daPy/du_(i+1,J-1)
    day.aPE          = day.aE + dFey;
    % daPy/du_{i,J-1}
    day.aPW          = day.aW - dFwy;
    % daPy/dv_{I,j-1}
    day.aPS          = day.aS - dFsy;
    % daPy/dv_{I,j+1}
    day.aPN          = day.aN + dFny;
    % daPy/dv_{I,j}
    day.aPP          = day.aN + dFny + day.aS - dFsy;
    
    % Define derivatives for linearized model with ay = -aPy v_{I,j} + aEy v_{I+1,J} + aWy v_{I-1,J} + aNy v_{I,j+1} + aSy v_{I,j-1}
    
    % day/du_(i,J-1) = daWy/du_{i,J-1} v_{I-1,j} - daPy/du_{i,J-1} v_{I,j}
    day.SW(2:Nx,1:Ny)  = day.aW(2:Nx,1:Ny).*v(1:Nx-1,1:Ny) - day.aPW(2:Nx,1:Ny).*v(2:Nx,1:Ny);
    
    % day/du_(i,J) = daWy/du_{i,J} v_{I-1,j} - daPy/du_{i,J} v_{I,j}
    day.NW(2:Nx,1:Ny)  = day.aW(2:Nx,1:Ny).*v(1:Nx-1,1:Ny) - day.aPW(2:Nx,1:Ny).*v(2:Nx,1:Ny);
    
    % day/du_(i+1,J-1) = daEy/du_{i+1,J-1} v_{I+1,j} - daPy/du_{i+1,J-1} v_{I,j}
    day.SE(1:Nx-1,1:Ny)=  day.aE(1:Nx-1,1:Ny).*v(2:Nx,1:Ny)  - day.aPE(1:Nx-1,1:Ny).*v(1:Nx-1,1:Ny);
    
    % day/du_(i+1,J) = daEy/du_{i+1,J} v_{I+1,j} - daPy/du_{i+1,J} v_{I,j}
    day.NE(1:Nx-1,1:Ny)=  day.aE(1:Nx-1,1:Ny).*v(2:Nx,1:Ny)  - day.aPE(1:Nx-1,1:Ny).*v(1:Nx-1,1:Ny);
    
    % day/dv_(I-1,j) = aWy
    day.W              = ay.aW;
    
    % day/dv_(I,j-1) = aSy + daSy/dv_{I,j-1} v_{I,j-1} - daPy/dv_{I,j-1} v_{I,j}
    day.S(1:Nx,2:Ny) = ay.aS(1:Nx,2:Ny) + day.aS(1:Nx,2:Ny).*v(1:Nx,1:Ny-1) - day.aPS(1:Nx,2:Ny).*v(1:Nx,2:Ny);
    
    % day/dv_(I,j) = daNy/dv_{I,j} v_{I,j+1} + daSy/dv_{I,j} v_{I,j-1} - daPy/dv_{I,j} v_{I,j} - aPy
    day.P(1:Nx,2:Ny-1) = -day.aN(1:Nx,2:Ny-1).*v(1:Nx,3:Ny) - day.aS(1:Nx,2:Ny-1).*v(1:Nx,1:Ny-2) + ...
        day.aPP(1:Nx,2:Ny-1).*v(1:Nx,2:Ny-1) + ay.aP(1:Nx,2:Ny-1);

    
    % day/dv_(I,j+1) = daNy/dv_{I,j+1} v_{I,j+1} + aNy - daPy/dv_{I,j+1} v_{I,j}
    day.N(1:Nx,1:Ny-1)   = ay.aN(1:Nx,1:Ny-1)   + day.aN(1:Nx,1:Ny-1).*v(1:Nx,2:Ny) - day.aPN(1:Nx,1:Ny-1).*v(1:Nx,1:Ny-1);
    
    % day/dv_(I+1,j) = aEy
    day.E              = ay.aE;
end;

%% Turbulence model
Turbulence    


output.ax = ax;
output.ay = ay;

if Linearversion
    output.dax = dax;
    output.day = day;
end;