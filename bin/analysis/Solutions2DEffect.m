%% Simulate the qlpv model obtained without projection

clear; clc; close all;

Method = 1;    % Method 1 (scaling factors) and Method 2 (only x- and y-momentum equations) 

%% Initialize script
options.Projection    = 0;                      % Use projection (true/false)
options.Linearversion = 0;                      % Provide linear variant of WFSim (true/false)
options.exportLinearSol= 0;                     % Calculate linear solution of WFSim
options.Derivatives   = 0;                      % Compute derivatives
options.startUniform  = 0;                      % Start from a uniform flowfield (true) or a steady-state solution (false)
options.exportPressures= ~options.Projection;   % Calculate pressure fields

%Wp.name       = 'SingleTurbine_50x50_lin';   % Meshing name (see "\bin\core\meshing.m")
Wp.name       = 'YawCase3_50x50_lin';   % Meshing name (see "\bin\core\meshing.m")

Animate       = 5;                      % Show 2D flow fields every x iterations (0: no plots)
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

% For the 2D effect
if Method ==1
    B2                    = 2*B2;
    Wp.turbine.forcescale = 1.5;
end

% Initialize variables and figure specific to this script
uk    = Wp.site.u_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
vk    = Wp.site.v_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
pk    = Wp.site.p_init*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
ubark = Wp.site.u_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);
vbark = Wp.site.v_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny,Wp.sim.NN);

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
    
    % Write flow field solutions to a 3D matrix
    uk(:,:,k) = sol.u;
    vk(:,:,k) = sol.v;
    pk(:,:,k) = sol.p;
    
    while ( eps>conv_eps && it<max_it && eps<epss );
        it   = it+1;
        epss = eps;
        
        if k>1
            max_it = max_it_dyn;
        end
        
        [sys,Power(:,k),Ueffect(:,k),a(:,k),CT(:,k)] = ...
            Make_Ax_b(Wp,sys,sol,input{k},B1,B2,bc,k,options); % Create system matrices
        
        
        if (k==1 && it==1)
            sol.x = [Wp.site.u_Inf*ones(Wp.Nu,1);Wp.site.v_Inf*ones(Wp.Nv,1);Wp.site.p_init*ones(Wp.Np,1)];
        end
        F                 = sys.M(sys.pRCM,sys.pRCM)*sol.x(sys.pRCM,1)+sys.m(sys.pRCM);
        sol.x(sys.pRCM,1) = sys.A(sys.pRCM,sys.pRCM)\F;
        
        % For 2D effect
        if Method==2
            if (k==1 && it==1)
                xbar   = [Wp.site.u_Inf*ones(Wp.Nu,1);Wp.site.v_Inf*ones(Wp.Nv,1)];
                [u,uu] = deal(Wp.site.u_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny));
                [v,vv] = deal(Wp.site.v_Inf*ones(Wp.mesh.Nx,Wp.mesh.Ny));
            end
            ubark(:,:,k)    = u;
            vbark(:,:,k)    = v;
            
            Fbar    = sys.M(1:Wp.Nu+Wp.Nv,1:Wp.Nu+Wp.Nv)*xbar+sys.m(1:Wp.Nu+Wp.Nv);
            xbar    = sys.A(1:Wp.Nu+Wp.Nv,1:Wp.Nu+Wp.Nv)\Fbar;
            
            uu(3:end-1,2:end-1)     = reshape(xbar(1:(Wp.mesh.Nx-3)*(Wp.mesh.Ny-2)),Wp.mesh.Ny-2,Wp.mesh.Nx-3)';
            vv(2:end-1,3:end-1)     = reshape(xbar((Wp.mesh.Nx-3)*(Wp.mesh.Ny-2)+1:(Wp.mesh.Nx-3)*(Wp.mesh.Ny-2)...
                +(Wp.mesh.Nx-2)*(Wp.mesh.Ny-3)),Wp.mesh.Ny-3,Wp.mesh.Nx-2)';
            
            if k==1; aleph      = min(1-.9^it,1); else aleph=1; end;
            u(3:end-1,2:end-1)  = (1-aleph)*u(3:end-1,2:end-1)+(aleph)*uu(3:end-1,2:end-1);
            v(2:end-1,3:end-1)  = (1-aleph)*v(2:end-1,3:end-1)+(aleph)*vv(2:end-1,3:end-1);
            
            % Update velocities for next iteration and boundary conditions
            % Three zero gradient boundary conditions
            u(:,1)          =  u(:,2);                  % u_{i,1}  = u_{i+1,2}   for i = 1,..Nx
            u(:,Wp.mesh.Ny) =  u(:,Wp.mesh.Ny-1);       % u_{i,Ny} = u_{i,Ny-1}  for i = 1,..Nx
            u(Wp.mesh.Nx,:) =  u(Wp.mesh.Nx-1,:);       % u_{Nx,J} = u_{Nx-1,J}  for J = 1,..Ny (hence u_Inf comes via first row in field)
            
            v(:,2)          =  v(:,3);
            v(:,1)          =  v(:,2);                  % v_{I,1}  = v_{I,2}    for I = 1,..Nx
            v(:,Wp.mesh.Ny) =  v(:,Wp.mesh.Ny-1);       % v_{I,Ny} = v_{I,Ny-1} for I = 1,..Nx
            v(Wp.mesh.Nx,:) =  v(Wp.mesh.Nx-1,:);       % v_{Nx,j} = v_{Nx,j}   for j = 1,..Ny (hence v_Inf comes via row in field)
        end
        
        [sol,eps] = MapSolution(Wp.mesh.Nx,Wp.mesh.Ny,sol,k,it,options);         % Map solution to field
        
    end
    toc
    
    if Animate > 0
        if ~rem(k,Animate)
            Animation;
        end;
    end;
end;
disp('Completed simulations.');


% Plot for 2D effect
for kk=1:1:Wp.sim.NN
    
    if Method==1
        figure(2);clf;
        plot(Wp.mesh.ldyy(1,:),uk(Wp.mesh.xline(end)+13,:,kk),'linewidth',2)        % Cross section wake at approximately 13*50 of 13*25
        grid; xlabel('y [m]');ylabel('u [m/s]')
        title('Cross section of the wake with scaling factors')
    end
    
    if Method==2
        figure(2);clf;
        subplot(1,2,1)
        plot(Wp.mesh.ldyy(1,:),uk(Wp.mesh.xline(end)+13,:,kk),'linewidth',2)        % Cross section wake at approximately 13*50 of 13*25
        grid; xlabel('y [m]');ylabel('u [m/s]')
        title('Cross section of the wake')
        subplot(1,2,2)
        plot(Wp.mesh.ldyy(1,:),ubark(Wp.mesh.xline(end)+13,:,kk),'linewidth',2)        % Cross section wake at approximately 13*50 of 13*25
        grid; xlabel('y [m]');ylabel('u [m/s]')
        title('Only momentum equations')
        
        figure(3);clf
        subplot(2,2,[1 3])
        contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',uk(:,:,kk),(0:0.1:Wp.site.u_Inf*1.2),'Linecolor','none');  colormap(hot); caxis([min(min(uk(:,:,kk)))-2 Wp.site.u_Inf*1.04]);  hold all; colorbar;
        axis equal; axis tight;
        xlabel('y [m]')
        ylabel('x [m]');
        title('u [m/s]');
        subplot(2,2,[2 4])
        contourf(Wp.mesh.ldyy(1,:),Wp.mesh.ldxx2(:,1)',ubark(:,:,kk),(0:0.1:Wp.site.u_Inf*1.2),'Linecolor','none');  colormap(hot); caxis([min(min(ubark(:,:,kk)))-2 Wp.site.u_Inf*1.04]);  hold all; colorbar;
        axis equal; axis tight;
        xlabel('y [m]')
        title('u [m/s]');
        drawnow
        
        figure(4);clf
        D_ind    = yline{1};
        cent_ind = ceil(0.5*(D_ind(1)+D_ind(end)));
        up(:,k)  = mean(ubark(:,D_ind,kk),2);
        
        plot(ldxx2(:,1)',up(:,k),'k','Linewidth',2);
        xlabel('x [m]');ylabel('U_c [m/s]');grid;
        ylim([0 u_Inf+1]);xlim([ldxx2(1,1) ldxx2(end,1)]);
        for kk=1:N
            vline(Crx(kk),'r--')
            vline(Crx(kk)+6*Dr,'b:','6D')
        end
        hold off;
        drawnow;
    end
    
end
