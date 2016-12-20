Ar         = Wp.turbine.Drotor;
yaw_angles = .5*Ar*exp(1i*Phi(:,k)*pi/180);  % Yaw angles

if (Dynamic==0 && Wp.N==1)
    % Compute center line wake
    kd = .15;
    ad = -4.5;
    bd = -0.01;
    
    for kk=1:Wp.N
        X        = Wp.ldxx(Wp.xline(kk),1);
        x        = Wp.ldxx(Wp.xline(kk):end,1);
        temp     = 2*kd*(x-X)/Ar +1;
        xinit    = .5*cos(Phi(kk)*pi/180)^2*sin(Phi(kk)*pi/180)*CT(kk);
        xw(:,kk) = (xinit*15*temp.^4 + xinit^2) ./(30*kd/Ar*temp.^5)-...
            xinit*Ar*(15+xinit^2)/(30*kd);%+ad+bd*(x-X);
    end
    
    % Compute theoretical velocity deficit for one turbine
    ke = .065;
    MU = 5.5;   % For mixing region 5.5, near wake region .5 and far wale region 1
    aU = 5;
    bU = 1.66;
    
    for kk=1:Wp.N
        X     = Wp.ldxx(Wp.xline(kk),1);
        x     = Wp.ldxx(Wp.xline(kk):end,1);
        mU    = MU/(cos(aU+bU*Phi(kk)*pi/180));
        c     = (Ar./(Ar+2*ke*mU*(x-X))).^2;
        %Vwake = Ueffect(1)*(1-2*a(kk)*c);
        
        kf            = .075;
        Vwake(:,kk)   = Ueffect(kk)*(1-(1-sqrt(1-CT(kk)))./(1+2*kf*(x-X)./Ar).^2);
        % Compute theoretical wake diameter
        Dwm(:,kk) = max(Ar+2*ke*1*(x-X),0); % mixing region
        Dwf(:,kk) = max(Ar+2*ke*.22*(x-X),0); % far wake
        Dwn(:,kk) = max(Ar+2*ke*-.5*(x-X),0); % near wake
    end
    
end

%% Check Velocity profile behind one turbine using AD theory
if (Dynamic==0 && Wp.N==1)
    
    D_ind    = Wp.yline{1};
    cent_ind = ceil(0.5*(D_ind(1)+D_ind(end)));
    up       = mean(u(:,D_ind),2);
    
    figure
    plot(Wp.ldxx2(:,1)',up,'k','Linewidth',2); hold on
    
    % Axial induction factor formula based on the actuator disc (Look at Torres paper)
    x = -(Wp.ldxx2(:,cent_ind)-Wp.ldxx2(Wp.xline,cent_ind));
    V = u_Inf*(1-(2/pi)*a*((pi/2)+atan(2*x/Ar)));
    plot(x+Wp.ldxx2(Wp.xline)*ones(length(x),1),V,'Color',[0 0 0]+0.4,'Linewidth',2);
    plot(Wp.ldxx(Wp.xline(1):end,1),Vwake(:,1),'r--','Linewidth',2);
    xlabel('distance [m]'); ylabel('velocity [m/s]');
    legend('Model','AD theory','Jensen');
    axis([0 Wp.ldxx2(end,1) ceil(u_Inf/3)-1 u_Inf]);
    
    % Location of the rotor disc in the grid (for demonstration in plot)
    l = 2:1:11;
    plot(Wp.ldxx2(Wp.xline)*ones(1,length(l)),l,'k--','Linewidth',2); grid on;
    
    % Partitioning x-direction by rotor diameter (for demonstration in plot)
    for dd=2:3:11
        plot((Wp.ldxx2(Wp.xline)+dd*Ar)*ones(1,length(l)),l,'k--');
        annotation('textbox',...
            [0.44+0.031*(Ar/90)*((dd-1)/2) 0.735 0.03 0.05],...
            'String',{sprintf('%dD', dd)},...
            'FontSize',8,...
            'FontName','Arial',...
            'LineStyle','--',...
            'EdgeColor',[1 1 1],...
            'LineWidth',2,...
            'BackgroundColor',[1 1 1],...
            'Color',[0 0 0]);
    end
    
    % AD information in the plot
    X = 0:10:5000;
    plot(X,u_Inf*ones(length(X),1),'k--','Linewidth',1);
    annotation('textbox',...
        [0.150 0.85 0.1 0.05],...
        'String',{'U_\infty'},...
        'FontSize',12,...
        'FontName','Arial',...
        'LineStyle','--',...
        'EdgeColor',[1 1 1],...
        'LineWidth',2,...
        'BackgroundColor',[1 1 1],...
        'Color',[0 0 0]);
    plot(X,u_Inf*(1-a)*ones(length(X),1),'k--','Linewidth',1);
    annotation('textbox',...
        [0.150 0.55 0.1 0.05],...
        'String',{'(1-a)U_\infty'},...
        'FontSize',12,...
        'FontName','Arial',...
        'LineStyle','--',...
        'EdgeColor',[1 1 1],...
        'LineWidth',2,...
        'BackgroundColor',[1 1 1],...
        'Color',[0 0 0]);
    plot(X,u_Inf*(1-2*a)*ones(length(X),1),'k--','Linewidth',1);
    annotation('textbox',...
        [0.150 0.16 0.1 0.05],...
        'String',{'(1-2a)U_\infty'},...
        'FontSize',12,...
        'FontName','Arial',...
        'LineStyle','--',...
        'EdgeColor',[1 1 1],...
        'LineWidth',2,...
        'BackgroundColor',[1 1 1],...
        'Color',[0 0 0]);
end

%% Energy
if (Dynamic==0 && Wp.N==1)
    Mom = 0.5*Rho*(u.*u+v.*v)+pp;
    figure
    contourf(Wp.ldyy2(1,:),Wp.ldxx2(:,1)',min(Mom,max(max(Mom))*1.2),[0:5:max(max(Mom))*1.2],'Linecolor','none');
    colorbar; hold on
    
    Qy     = (Wp.Cry(1)-real(yaw_angles(1))):1:(Wp.Cry(1)+real(yaw_angles(1)));
    Qx     = linspace(Wp.Crx(1)-imag(yaw_angles(1)),Wp.Crx(1)+imag(yaw_angles(1)),length(Qy));
    plot(Qy,Qx,'k','linewidth',1);
    plot(Wp.Cry*ones(length(Wp.ldxx2(:,1)'),1)+2*Ar,Wp.ldxx2(:,1)','k--','linewidth',1);
    plot(Wp.Cry*ones(length(Wp.ldxx2(:,1)'),1)-2*Ar,Wp.ldxx2(:,1)','k--','linewidth',1);
    plot(Wp.Cry*ones(length(Wp.ldxx2(:,1)'),1)+.5*Ar,Wp.ldxx2(:,1)','k--','linewidth',1);
    plot(Wp.Cry*ones(length(Wp.ldxx2(:,1)'),1)-.5*Ar,Wp.ldxx2(:,1)','k--','linewidth',1);
    plot(Wp.Cry(1)-xw(:,1)+Dwn(:,1)/2,Wp.ldxx(Wp.xline(1):end,1),'r--');
    plot(Wp.Cry(1)-xw(:,1)-Dwn(:,1)/2,Wp.ldxx(Wp.xline(1):end,1),'r--');
    plot(Wp.Cry(1)-xw(:,1)+Dwf(:,1)/2,Wp.ldxx(Wp.xline(1):end,1),'r--');
    plot(Wp.Cry(1)-xw(:,1)-Dwf(:,1)/2,Wp.ldxx(Wp.xline(1):end,1),'r--');
    plot(Wp.Cry(1)-xw(:,1)+Dwm(:,1)/2,Wp.ldxx(Wp.xline(1):end,1),'r--');
    plot(Wp.Cry(1)-xw(:,1)-Dwm(:,1)/2,Wp.ldxx(Wp.xline(1):end,1),'r--');
    xlabel('y-direction [m]'); ylabel('x-direction [m]'); title('Flow energy');
    
end

if Dynamic==0
    figure
    subplot(2,2,1)
    plot3(Wp.ldxx,Wp.ldyy,pp);zlabel('Pressure');xlabel('x');ylabel('y');grid;axis tight
    subplot(2,2,2)
    plot3(Wp.ldxx,Wp.ldyy,u);zlabel('u-velocity');xlabel('x');ylabel('y');grid;axis tight
    subplot(2,2,3)
    plot3(Wp.ldxx,Wp.ldyy,v);zlabel('v-velocity');xlabel('x');ylabel('y');grid;axis tight
end

%% Vector plot velocities
if (Dynamic==0 && Wp.N<=2)
    
    figure();
    [~,ind] = min(abs(Wp.ldxx(Wp.xline(1))+Ar-Wp.ldxx(:,1)));
    subplot(1,7,1);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{1})),Wp.ldyy2(1,Wp.yline{1}),'r-','linewidth',3);
    grid;ylabel('y [m]');xlabel('u [m/s] at 1D');
    [~,ind] = min(abs(Wp.ldxx(Wp.xline(1))+2*Ar-Wp.ldxx(:,1)));
    subplot(1,7,2);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{1})),Wp.ldyy2(1,Wp.yline{1}),'r-','linewidth',3);
    grid;xlabel('u [m/s] at 2D');
    [~,ind] = min(abs(Wp.ldxx(Wp.xline(1))+3*Ar-Wp.ldxx(:,1)));
    subplot(1,7,3);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{1})),Wp.ldyy2(1,Wp.yline{1}),'r-','linewidth',3);
    grid;xlabel('u [m/s] at 3D');
    [~,ind] = min(abs(Wp.ldxx(Wp.xline(1))+4*Ar-Wp.ldxx(:,1)));
    subplot(1,7,4);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{1})),Wp.ldyy2(1,Wp.yline{1}),'r-','linewidth',3);
    grid;xlabel('u [m/s] at 4D');title('u-velocity profile turbine 1')
    [~,ind] = min(abs(Wp.ldxx(Wp.xline(1))+5*Ar-Wp.ldxx(:,1)));
    subplot(1,7,5);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{1})),Wp.ldyy2(1,Wp.yline{1}),'r-','linewidth',3);
    grid;xlabel('u [m/s] at 5D');
    [~,ind] = min(abs(Wp.ldxx(Wp.xline(1))+6*Ar-Wp.ldxx(:,1)));
    subplot(1,7,6);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{1})),Wp.ldyy2(1,Wp.yline{1}),'r-','linewidth',3);
    grid;xlabel('u [m/s] at 6D');
    [~,ind] = min(abs(Wp.ldxx(Wp.xline(1))+9*Ar-Wp.ldxx(:,1)));
    subplot(1,7,7);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{1})),Wp.ldyy2(1,Wp.yline{1}),'r-','linewidth',3);
    grid;xlabel('u [m/s] at 9D');
    
    if Wp.N==2
        figure();
        [~,ind] = min(abs(Wp.ldxx(Wp.xline(2))+Ar-Wp.ldxx(:,1)));
        subplot(1,7,1);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{2})),Wp.ldyy2(1,Wp.yline{2}),'r-','linewidth',3);
        grid;ylabel('y [m]');xlabel('u [m/s] at 1D');
        [~,ind] = min(abs(Wp.ldxx(Wp.xline(2))+2*Ar-Wp.ldxx(:,1)));
        subplot(1,7,2);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{2})),Wp.ldyy2(1,Wp.yline{2}),'r-','linewidth',3);
        grid;xlabel('u [m/s] at 2D');
        [~,ind] = min(abs(Wp.ldxx(Wp.xline(2))+3*Ar-Wp.ldxx(:,1)));
        subplot(1,7,3);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{2})),Wp.ldyy2(1,Wp.yline{2}),'r-','linewidth',3);
        grid;xlabel('u [m/s] at 3D');
        [~,ind] = min(abs(Wp.ldxx(Wp.xline(2))+4*Ar-Wp.ldxx(:,1)));
        subplot(1,7,4);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{2})),Wp.ldyy2(1,Wp.yline{2}),'r-','linewidth',3);
        grid;xlabel('u [m/s] at 4D');title('u-velocity profile turbine 2')
        [~,ind] = min(abs(Wp.ldxx(Wp.xline(2))+5*Ar-Wp.ldxx(:,1)));
        subplot(1,7,5);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{2})),Wp.ldyy2(1,Wp.yline{2}),'r-','linewidth',3);
        grid;xlabel('u [m/s] at 5D');
        [~,ind] = min(abs(Wp.ldxx(Wp.xline(2))+6*Ar-Wp.ldxx(:,1)));
        subplot(1,7,6);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{2})),Wp.ldyy2(1,Wp.yline{2}),'r-','linewidth',3);
        grid;xlabel('u [m/s] at 6D');
        [~,ind] = min(abs(Wp.ldxx(Wp.xline(2))+9*Ar-Wp.ldxx(:,1)));
        subplot(1,7,7);plot(u(ind,:),Wp.ldyy2(1,:));hold on;plot(zeros(1,length(Wp.yline{2})),Wp.ldyy2(1,Wp.yline{2}),'r-','linewidth',3);
        grid;xlabel('u [m/s] at 9D');
    end
    
    figure()
    for kk=1:Wp.N
        Qy     = (Wp.Cry(kk)-real(yaw_angles(kk))):1:(Wp.Cry(kk)+real(yaw_angles(kk)));
        Qx     = linspace(Wp.Crx(kk)-imag(yaw_angles(kk)),Wp.Crx(kk)+imag(yaw_angles(kk)),length(Qy));
        plot(Qy,Qx,'k','linewidth',1);hold on;
    end
    quiver(Wp.ldyy2,Wp.ldxx2,v,u);grid;
    xlabel('y-direction [m]');ylabel('x-direction [m]'); axis tight
end