function [StrucDiscretization,StrucBCs] = BoundaryConditions(Wp,StrucDiscretization,sol,Linearversion)
% Import variables
Nx = Wp.mesh.Nx;
Ny = Wp.mesh.Ny;
u  = sol.u;
v  = sol.v;

% Execute script
ax = StrucDiscretization.ax;
ay = StrucDiscretization.ay;

bbx = sparse((Nx-3)*(Ny-2),(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2));
bby = sparse((Nx-2)*(Ny-3),(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2));

% Zero gradient outflow
% for u-direction
ax.aP((Nx-1),(2:Ny-1))  = ax.aP((Nx-1),(2:Ny-1)) - ax.aE((Nx-1),(2:Ny-1)); %NORTH
ax.aP((1:Nx-1),(Ny-1))  = ax.aP((1:Nx-1),(Ny-1)) - ax.aN((1:Nx-1),(Ny-1)); %EAST
ax.aP((1:Nx-1),(2))     = ax.aP((1:Nx-1),2)      - ax.aS((1:Nx-1),2);

if Linearversion
    dax = StrucDiscretization.dax;
    day = StrucDiscretization.day;
    
    dax.P((Nx-1),(2:Ny-1))  = dax.P((Nx-1),(2:Ny-1)) - dax.E((Nx-1),(2:Ny-1)); %NORTH
    dax.P((1:Nx-1),(Ny-1))  = dax.P((1:Nx-1),(Ny-1)) - dax.N((1:Nx-1),(Ny-1)); 
    dax.P((1:Nx-1),(2))     = dax.P((1:Nx-1),2)      - dax.S((1:Nx-1),2);

    %% inflow conditions (no zero gradient this seems to work)
    dax.P(2:3,:)=dax.P(2:3,:)+dax.aW(2:3,1:Ny).*u(1:2,1:Ny); % y momentum west side

    %% zero gradient
    dax.NW(:,2) =dax.NW(:,2)-dax.NW(:,1);%+  dax.aPN(:,2).*u(:,2); % a bit random (i can not explain)
    dax.NE(1:Nx,2) = dax.NE(1:Nx,2) -dax.NE(1:Nx,1) ; % a bit random (i can not explain)
    dax.SW(:,end-1) =dax.SW(:,end-1)-dax.SW(:,end);%-  dax.aS(1:Nx,end-1).*u(1:Nx,end-1) ; % a bit random (i can not explain)
    dax.SE(1:Nx,end-1) = dax.SE(1:Nx,end-1)-dax.SE(1:Nx,end);% - dax.aS(1:Nx,end-1).*u(1:Nx,end-1); % a bit random (i can not explain)
    day.SW(end-1,:) =day.SW(end-1,:)-day.SW(end,:);%  day.aW(end-1,:).*v(end-1,:) ; % a bit random (i can not explain)
    day.NW(end-1,:) =day.NW(end-1,:)-day.NW(end,:);%  day.aW(end-1,:).*v(end-1,:) ; % a bit random (i can not explain)

    %dax.NW(2:3,1:Ny-1) = dax.NW(2:3,1:Ny-1)-dax.aS(2:3,1:Ny-1).*u(2:3,2:Ny) ; %JW
end;

%% extra zero gradient derivatives
% dax.NW(1:Nx-1,2) =dax.NW(1:Nx-1,2)+dax.aS(1:Nx-1,2).*u(2:Nx,2);% dax.aN(1:Nx-1,1:Ny-1).*u(1:Nx-1,2:Ny)  - dax.aPN(1:Nx-1,1:Ny-1).*u(1:Nx-1,1:Ny-1); %JW
% dax.NE(1:Nx-1,2) = dax.NE(1:Nx-1,2)+ dax.aS(1:Nx-1,2).*u(2:Nx,2) ;

%% for v-direction
ay.aP((Nx-1),(1:Ny))  = ay.aP((Nx-1),(1:Ny)) - ay.aE((Nx-1),(1:Ny));
ay.aP((1:Nx),(Ny-1))  = ay.aP((1:Nx),(Ny-1)) - ay.aN((1:Nx),(Ny-1));
ay.aP((1:Nx),(3))     = ay.aP((1:Nx),3)      - ay.aS((1:Nx),3); % changed to 3 3 2 instead of 2 2 1

if Linearversion
    day.P((Nx-1),(1:Ny))  = day.P((Nx-1),(1:Ny)) - day.E((Nx-1),(1:Ny));
    day.P((1:Nx),(Ny-1))  = day.P((1:Nx),(Ny-1)) - day.N((1:Nx),(Ny-1));
    day.P((1:Nx),(3))     = day.P((1:Nx),3)      - day.S((1:Nx),3);% changed to 3 3 2 instead of2 2 1
end;

% Inflow boundary for non linear model
bx      = kron([1;zeros(Nx-4,1)],(ax.aW(3,2:end-1).*u(2,2:end-1))'); %changed to 3: 2 instead of 2:2
by      = [v(1,3:Ny-1)'.*ay.aW(2,3:Ny-1)';zeros((Nx-3)*(Ny-3),1)]; %changed to 2:3 inst


if Linearversion
% Inflow boundary for linear model
bbx(1:Ny-2,1:Ny-2) = diag((dax.aW(3,2:end-1).*u(2,2:end-1))');
bby(1:Ny-2,1:Ny-2) = 0.*diag((day.aW(2,2:end-1).*v(1,2:end-1))'); % I dont understand this but is works :)_

% Write to output
StrucBCs.dbcdx = [bbx;bby;sparse((Nx-2)*(Ny-2),(Nx-3)*(Ny-2)+(Nx-2)*(Ny-3)+(Nx-2)*(Ny-2))];
StrucDiscretization.dax = dax;
StrucDiscretization.day = day;
end;

% Write nonlinear to outputs
StrucDiscretization.ax  = ax;
StrucDiscretization.ay  = ay;
StrucBCs.bx             = bx;
StrucBCs.by             = by;
end