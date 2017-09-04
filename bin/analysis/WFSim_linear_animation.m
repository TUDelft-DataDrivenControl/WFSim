function [hfig] = WFSim_linear_animation( Wp,sol,hfig )
    % Import local variables from large structs
    Dr     = Wp.turbine.Drotor;
    Nx     = Wp.mesh.Nx;
    Ny     = Wp.mesh.Ny;
    dxx    = Wp.mesh.dxx;
    dyy    = Wp.mesh.dyy;
    dxx2   = Wp.mesh.dxx2;
    dyy2   = Wp.mesh.dyy2;
    ldxx   = Wp.mesh.ldxx;
    ldyy   = Wp.mesh.ldyy;
    ldxx2  = Wp.mesh.ldxx2;
    ldyy2  = Wp.mesh.ldyy2;
    xline  = Wp.mesh.xline;
    yline  = Wp.mesh.yline;

    Rho    = Wp.site.Rho;
    mu     = Wp.site.mu;
    lmu    = Wp.site.lmu;
    Turb   = Wp.site.turbul;
    u_Inf  = Wp.site.u_Inf;
    v_Inf  = Wp.site.v_Inf;

    Drotor = Wp.turbine.Drotor;
    N      = Wp.turbine.N;
    Cry    = Wp.turbine.Cry;
    Crx    = Wp.turbine.Crx;
    input  = Wp.turbine.input(sol.k);

    time   = Wp.sim.time;
    k      = sol.k;
    
    turb_coord = .5*Dr*exp(1i*input.phi*pi/180);  % Yaw angles

    %% Create figure window, if not yet available
    if nargin <= 2
        scrsz = get(0,'ScreenSize');
        hfig = figure('color',[0 166/255 214/255],'units','normalized','outerposition',...
        [0 0 1 1],'ToolBar','none','visible', 'on');
    end;
    
    set(0,'CurrentFigure',hfig);

    %% Plot u velocity flow component from nonlinear model
    subplot(2,2,[1 3]);
    contourf(ldyy(1,:),ldxx2(:,1)',sol.u,(0:0.1:u_Inf*1.2),'Linecolor','none');  colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.04]);  hold all; colorbar;
    axis equal; axis tight;

    % Plot the turbines in the field
    for kk=1:N
        Qy     = (Cry(kk)-real(turb_coord(kk))):1:(Cry(kk)+real(turb_coord(kk)));
        Qx     = linspace(Crx(kk)-imag(turb_coord(kk)),Crx(kk)+imag(turb_coord(kk)),length(Qy));
        plot(Qy,Qx,'k','linewidth',1)
    end
    text(0,ldxx2(end,end)+80,['Time ', num2str(time(k),'%.1f'), 's']);
    xlabel('y [m]')
    ylabel('x [m]');
    title('u [m/s]');
    hold off;

    %% Plot the u velocity flow component from linear model
    subplot(2,2,2);
    contourf(ldyy(1,:),ldxx2(:,1)',sol.ul,(0:0.1:u_Inf*1.2),'Linecolor','none');  colormap(hot); caxis([min(min(sol.u))-2 u_Inf*1.04]);  hold all; colorbar;
    colorbar;
    for kk=1:N
        Qy     = (Cry(kk)-real(turb_coord(kk))):1:(Cry(kk)+real(turb_coord(kk)));
        Qx     = linspace(Crx(kk)-imag(turb_coord(kk)),Crx(kk)+imag(turb_coord(kk)),length(Qy));
        plot(Qy,Qx,'k','linewidth',1)
    end
    axis equal; axis tight
    xlabel('y [m]')
    ylabel('x [m]');
    title('u_l [m/s]')
    hold off;

    %% Plot the error between the linear and nonlinear model
    e               = sol.u-sol.ul;
    [maxe,maxeloc]  = max(abs(vec(e)));

    set(0,'CurrentFigure',hfig);
    subplot(2,2,4)
    contourf(ldyy(1,:),ldxx2(:,1)',e,(-2:0.1:2),'Linecolor','none');  colormap(hot); caxis([-2 2]);  hold all; colorbar;
    colorbar;
    for kk=1:N
        Qy     = (Cry(kk)-real(turb_coord(kk))):1:(Cry(kk)+real(turb_coord(kk)));
        Qx     = linspace(Crx(kk)-imag(turb_coord(kk)),Crx(kk)+imag(turb_coord(kk)),length(Qy));
        plot(Qy,Qx,'k','linewidth',1)
    end
    axis equal; axis tight
    xlabel('y [m]')
    ylabel('x [m]');
    title(['e [m/s] with max(e) = ', num2str(maxe,'%.1f'), '[m/s]'])
    hold off;
    drawnow;
end