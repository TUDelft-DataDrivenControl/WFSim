function [output,Ueffect,a,Power,CT] = Actuator_Joeri(Wp,input,sol,options)

Nx          = Wp.mesh.Nx;
Ny          = Wp.mesh.Ny;
dyy2        = Wp.mesh.dyy2;
xline       = Wp.mesh.xline;
yline       = Wp.mesh.yline;
ylinev      = Wp.mesh.ylinev;

Rho        = Wp.site.Rho;

Drotor     = Wp.turbine.Drotor;
N          = Wp.turbine.N;
F          = Wp.turbine.forcescale; % http://www.nrel.gov/docs/fy05osti/36834.pdf
powerscale = Wp.turbine.powerscale;

u       = sol.u;
v       = sol.v;
Phi     = input.phi;
beta    = input.beta;
dPhi    = input.dphi;
dbeta   = input.dbeta;

Sm.x     = sparse(Nx-3,Ny-2);
Sm.dx    = sparse(Nx-3,Ny-2);

Sm.xx    = sparse((Nx-3)*(Ny-2),N);
Sm.dxx   = sparse((Nx-3)*(Ny-2),N);
Sm.y     = sparse(Nx-2,Ny-3);
Sm.dy    = sparse(Nx-2,Ny-3);
tempy    = zeros(Nx-2,Ny-3);
tempdy   = zeros(Nx-2,Ny-3);
Sm.yy    = sparse((Nx-2)*(Ny-3),N);
Sm.dyy   = sparse((Nx-2)*(Ny-3),N);

Smdu     = sparse(Nx-3,Ny-2);
Smdv     = sparse(Nx-2,Ny-3);

dSm.xdu     = sparse((Nx-3)*(Ny-2),(Nx-3)*(Ny-2));
dSm.xdv     = sparse((Nx-3)*(Ny-2),(Nx-2)*(Ny-3) );
dSm.ydu     = sparse((Nx-2)*(Ny-3),(Nx-3)*(Ny-2));
dSm.ydv     = sparse((Nx-2)*(Ny-3),(Nx-2)*(Ny-3));

dSmxdbetakk = sparse(Nx-3,Ny-2);
dSmydbetakk = sparse(Nx-2,Ny-3);

dSmxdphikk  = sparse(Nx-3,Ny-2);
dSmydphikk  = sparse(Nx-2,Ny-3);

JXXu  = sparse(Nx-3,Ny-2);
JXXv  = sparse(Nx-2,Ny-3);
JXXp  = sparse(Nx-2,Ny-2);

for kk=1:N
    %% Oppassen met multiple turbine dSmxdbetakk  dSmxdphikk
    clear dSmxdbetakk  dSmxdphikk dSmydbetakk  dSmydphikk tempx tempdx
    tempx    = zeros(Nx-3,Ny-2);
    tempdx   = zeros(Nx-3,Ny-2);
    
    dSmxdbetakk = sparse(Nx-3,Ny-2);
    dSmydbetakk = sparse(Nx-2,Ny-3);
    
    dSmxdphikk  = sparse(Nx-3,Ny-2);
    dSmydphikk  = sparse(Nx-2,Ny-3);
    
    a(kk)         =     beta(kk)/(beta(kk)+1);                      % Axial induction factor at turbine
    vv            =     diff(v(xline(kk,:),ylinev{kk}))/2+v(xline(kk,:),ylinev{kk}(1:end-1));
    U(kk,:)       =     sqrt(u(xline(kk,:),yline{kk}).^2+vv.^2);
    
    
    dvvdv         = 1/2*diag(sign(ylinev{kk}))+1/2*diag(sign(ylinev{kk}(1:end-1)),1); dvvdv(end,:)=[];
    dUdv(kk,:)          =  vv./(U(kk,1:length(yline{kk})));
    dUdU(kk,1:length(yline{kk}))          =  u(xline(kk,:),yline{kk})./(U(kk,1:length(yline{kk})));
    phi(kk)       = 0;%-180/pi*atan(mean(v(xline(kk,:)-5,yline{kk}))/mean(u(xline(kk,:)-5,yline{kk}))); %% added an offset here
    Ue(kk,:)       =cos(atan(vv./(u(xline(kk,:),yline{kk})))+Phi(kk)/180*pi).*U(kk,1:length(yline{kk}));
    
    uu=u(xline(kk,:),yline{kk});
    dUedv(kk,:)=-sin(atan(vv./(u(xline(kk,:),yline{kk}))) + Phi(kk)/180*pi)./((u(xline(kk,:),yline{kk})).*(vv.^2./(u(xline(kk,:),yline{kk})).^2 + 1)).*U(kk,1:length(yline{kk}))+cos(atan(vv./(u(xline(kk,:),yline{kk})))+Phi(kk)/180*pi).*dUdv(kk,:); %% diff Ue wrt v
    dUedv(kk,:)=-uu.*sin(atan(vv./uu) + Phi(kk)/180*pi)./(vv.^2+uu.^2).*U(kk,1:length(yline{kk}))+cos(atan(vv./(u(xline(kk,:),yline{kk})))+Phi(kk)/180*pi).*dUdv(kk,:); %% diff Ue wrt v
    dUedU(kk,:)=(vv.*sin(atan(vv./(u(xline(kk,:),yline{kk}))) + Phi(kk)/180*pi))./((vv.^2 + (u(xline(kk,:),yline{kk})).^2)).*U(kk,1:length(yline{kk}))+cos(atan(vv./(u(xline(kk,:),yline{kk})))+Phi(kk)/180*pi).*dUdU(kk,1:length(yline{kk})); %% diff Ue wrt u
    dUedPhi(kk,:)=-1/180*pi*sin(atan(vv./(u(xline(kk,:),yline{kk}))) + Phi(kk)/180*pi).*U(kk,1:length(yline{kk})); %% diff wrt to Ue
    
    % http://www.nrel.gov/docs/fy05osti/36834.pdf
    if a(kk)> 0.4
        CT(kk)    = (8/9+(4*F-40/9)*a(kk)+(50/9-4*F)*a(kk)^2);%*cos(Phi(kk)*pi/180); %not a function of phi
        diff1 = 1/(beta(kk)+1).^2; % JW derivative of beta(kk)/(beta(kk)+1)
        diff2 = 2*beta(kk)/(beta(kk)+1)^3;  % JW derivative of (beta(kk)/(beta(kk)+1))^2)
        dCTdbeta(kk)    =((4*F-40/9)*diff1+(50/9-4*F)*diff2);%*cos(Phi(kk)*pi/180); % NEW JW
    else
        CT(kk)    =  4*a(kk)*F*(1-a(kk));%*cos(Phi(kk)*pi/180);
        diff1 = -(beta(kk)-1)/(beta(kk)+1)^3;% JW beta/(beta(kk)+1).^2; % derivative of beta(kk)/(beta(kk)+1)^2
        dCTdbeta(kk)    =  4*diff1*F;%*cos(Phi(kk)*pi/180); %  JW NEW JW
    end
    
    
    %% Gradient of the thrust with respect to beta
    Fthrust             = 1/2*Rho*CT(kk)*(Ue(kk,:)*(beta(kk)+1)).^2; % similar as 1/2*Rho*CT(kk)*(U/(1-a(kk))).^2; %somethin to try 2*60*ones(1,size(U(kk,1:length(yline{kk})),2));%
    dFthrustdbeta       = 1/2*Rho*CT(kk)*Ue(kk,:).^2*(2*beta(kk)+2)+1/2*Rho*dCTdbeta(kk)*Ue(kk,:).^2*(beta(kk)+1).^2; %diff to beta 1/2*Rho*CT(kk)*U^2*(beta(kk)+1).^2;
    dFxdbeta            = dFthrustdbeta*cos((Phi(kk)+phi(kk))*pi/180);
    dFydbeta            = dFthrustdbeta*sin((Phi(kk)+phi(kk))*pi/180);
    
    dSmxdbetakk(xline(kk,:)-2,yline{kk}-1) =  -(dFxdbeta.*dyy2(1,yline{kk}))';  %JW changed definition of force and location (need correction for yaw
    scale=3;
    dSmydbetakk(xline(kk,:)-1,yline{kk}(2:end)-2) = scale*(dFydbeta(2:end).*dyy2(1,yline{kk}(2:end)))'; %% JW made this 0
    
    dSm.dbeta(:,kk)   = [vec(dSmxdbetakk');vec(dSmydbetakk')];
    
    %% Gradient of the thrust with respect to yaw angle [to be checked]
    dFthrustdPhi                = Rho*CT(kk)*Ue(kk,:).*((beta(kk)+1)).^2.*dUedPhi(kk,:);%%1/2*Rho*dCTdPhi(kk)*(Ue(kk,:)*(beta(kk)+1)).^2;
    dFxdPhi                     = dFthrustdPhi*cos((Phi(kk)+phi(kk))*pi/180)-pi/180*Fthrust*sin((Phi(kk)+phi(kk))*pi/180);
    dFydPhi                     = +dFthrustdPhi*sin((Phi(kk)+phi(kk))*pi/180)+pi/180*Fthrust*cos((Phi(kk)+phi(kk))*pi/180);
    
    dSmxdphikk(xline(kk,:)-2,yline{kk}-1) = -dFxdPhi'.*dyy2(1,yline{kk})';
    dSmydphikk(xline(kk,:)-1,yline{kk}(2:end)-2) = scale*dFydPhi(2:end)'.*dyy2(1,yline{kk}(2:end))';
    
    dSm.dphi(:,kk)   = [vec(dSmxdphikk');vec(dSmydphikk')];
    %% Gradient of the Thrust with respect to the state
    dFxdU = Rho*CT(kk)*Ue(kk,:)*((beta(kk)+1)).^2*cos((phi(kk)+Phi(kk))*pi/180);
    dFydU = Rho*CT(kk)*Ue(kk,:)*((beta(kk)+1)).^2*sin((phi(kk)+Phi(kk))*pi/180);
    dSm.xdu((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-3)*(Ny-2)+yline{kk}-1) = diag((-dFxdU.*dUedU(kk,1:length(yline{kk})).*dyy2(1,yline{kk})));
    %dSm.xdv((xline(kk,:)-3)*(Ny-2)+ylinev{kk}-1,(xline(kk,:)-2)*(Ny-3)+ylinev{kk}-2) = diag((-dFxdU.*dUedv(kk,:).*dyy2(1,yline{kk})*dvvdv));
    dSm.xdv((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-2)*(Ny-3)+yline{kk}-2) = 1/2*diag((-dFxdU.*dUedv(kk,:).*dyy2(1,yline{kk})));
    dSm.xdv((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-2)*(Ny-3)+yline{kk}-1) = dSm.xdv((xline(kk,:)-3)*(Ny-2)+yline{kk}-1,(xline(kk,:)-2)*(Ny-3)+yline{kk}-1)+ 1/2*diag((-dFxdU.*dUedv(kk,:).*dyy2(1,yline{kk})));
    dSm.ydu((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-3)*(Ny-2)+yline{kk}(2:end)-1) = diag(scale*(dFydU(2:end).*dUedU(kk,2:end).*dyy2(1,yline{kk}(2:end))));
    %dSm.ydv((xline(kk,:)-2)*(Ny-3)+yline{kk}-2,(xline(kk,:)-2)*(Ny-3)+yline{kk}-2) = diag(scale*(dFydU(2:end).*dUedv(kk,2:end).*dyy2(1,yline{kk}(2:end))*dvvdv(2:end,2:end)));
    dSm.ydv((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2) = 1/2*diag(scale*(dFydU(2:end).*dUedv(kk,2:end).*dyy2(1,yline{kk}(2:end))));
    dSm.ydv((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-1) =   dSm.ydv((xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-2,(xline(kk,:)-2)*(Ny-3)+yline{kk}(2:end)-1)+ 1/2*diag(scale*(dFydU(2:end).*dUedv(kk,2:end).*dyy2(1,yline{kk}(2:end))));
    
    %%
    dvvdv         = 1/2*diag(sign(ylinev{kk}))+1/2*diag(sign(ylinev{kk}(1:end-1)),1); dvvdv(end,:)=[];
    dUdv(kk,:)    =  vv./U(kk,:);
    dUdu(kk,:)    =  uu./U(kk,:);
    
    dUedv(kk,:)=-sin(atan(vv./(uu)) + input.phi(kk)/180*pi)./((uu).*...
        (vv.^2./(uu).^2 + 1)).*U(kk,:)+cos(atan(vv./(uu))+input.phi(kk)/180*pi).*dUdv(kk,:); %% diff Ue wrt v
    dUedv(kk,:)=-uu.*sin(atan(vv./uu) + input.phi(kk)/180*pi)./(vv.^2+uu.^2).*U(kk,:)+cos(atan(vv./(uu))+...
        input.phi(kk)/180*pi).*dUdv(kk,:); %% diff Ue wrt v
    dUedu(kk,:)=(vv.*sin(atan(vv./(uu)) + input.phi(kk)/180*pi))./((vv.^2 + (uu).^2)).*...
        U(kk,:)+cos(atan(vv./(uu))+input.phi(kk)/180*pi).*dUdu(kk,:); %% diff Ue wrt u
    dSm.dbetadPower_in(kk) = 1/(1/F*powerscale*2*Rho*pi*(0.5*Drotor)^2*mean(Ue(kk,:)).^3);
    dSm.dJdPower_in(kk)    = -1/F*powerscale*pi*(0.5*Drotor)^2*(2*Rho).*mean(Ue(kk,:)).^3*dSm.dbetadPower_in(kk);
    dSm.dJdPower_in(kk)   = -1/F*powerscale*pi*(0.5*Drotor)^2*(2*Rho).*mean(Ue(kk,:)).^3*dSm.dbetadPower_in(kk);
    
    % Derivative cost function with respect to state and beta
    JXXu(xline(kk,:)-2,yline{kk}-1) =  -1/F*pi*(0.5*Drotor)^2*6*Rho*input.beta(kk).*mean(Ue(kk,:)).^2.*dUedu(kk,:)./(length(yline{kk}));%-3*mean(urotor).^2*pi*(0.5*turbine.Drotor)^2*2*input.beta(kk)/(length(solind)/2);
    JXXv(xline(kk,:)-1,ylinev{kk}-2) = -1/F*pi*(0.5*Drotor)^2*6*Rho*input.beta(kk).*mean(Ue(kk,:)).^2.*dUedv(kk,:)./((length(yline{kk})))*dvvdv;%-3*mean(urotor).^2*pi*(0.5*turbine.Drotor)^2*2*input.beta(kk)/(length(solind)/2);
    
    % two lines above should be checked!!
    dSm.dJdx        = [vec(JXXu');vec(JXXv');vec(JXXp')];
    dSm.dJdbeta(kk) = -1/F*pi*(0.5*Drotor)^2*(2*Rho).*mean(Ue(kk,:)).^3;
    dSm.dJdPhi(kk)  = -1/F*6*Rho*pi*(0.5*Drotor)^2*input.beta(kk)*mean(Ue(kk,:)).^2.*mean(dUedPhi(kk,:));
    %%
    
    %     dSm.xdu(xline(kk,:)-2,yline{kk}-1)    = (-dFxdU.*dUedu.*dyy2(1,yline{kk}));
    %     dSm.xdv(xline(kk,:)-1,ylinev{kk}-2)    = (-dFxdU.*dUedv.*dyy2(1,yline{kk})*dvvdv)';
    %     dSm.ydu(xline(kk,:)-2,yline{kk}(2:end)-1)    =scale*(dFydU(2:end).*dUedu(2:end).*dyy2(1,yline{kk}(2:end)))';
    %     dSm.ydv(xline(kk,:)-1,yline{kk}(1:end)-2)    =scale*(dFydU(2:end).*dUedv(2:end).*dyy2(1,yline{kk}(2:end))*dvvdv(2:end,2:end))';
    %
    %     dSm.xxdu    = diag(vec(dSm.xdu')');
    %     dSm.xxdv    = diag(vec(dSm.xdv')');
    %     dSm.yydu    = diag(vec(dSm.ydu')');
    %     dSm.yydv    = diag(vec(dSm.ydv')');
    %%
    %  Fthrust       = 1/2*Rho*CT(kk)*(U/(1-a(kk))).^2;    % dFhrust has to change when using this Fthrust
    
    dFthrust      = 2*Rho*Ue(kk,:).^2*cos((Phi(kk)+phi(kk))*pi/180); %% JW what is this???
    
    Fx            = Fthrust*cos((phi(kk)+Phi(kk))*pi/180);
    dFx           = dFthrust*cos((phi(kk)+Phi(kk))*pi/180);
    Fy            = Fthrust*sin((phi(kk)+Phi(kk))*pi/180);
    dFy           = dFthrust*sin((phi(kk)+Phi(kk))*pi/180);
    
    Fxmean(kk)    = mean(-Fx'.*dyy2(1,yline{kk})');
    Fymean(kk)    = mean(Fy'.*dyy2(1,yline{kk})');
    
    %%
    CP(kk)        = 4*beta(kk)/(beta(kk)+1)^3*0.768*cos((Phi(kk)-phi(kk))*pi/180)^(1.88);
    Power2(kk)     = 1/F*mean(2*Rho*pi*(0.5*Drotor)^2*beta(kk)*Ue(kk,:).^3); %% JW ??
    Power(kk)     = 1/F*2*Rho*pi*(0.5*Drotor)^2*beta(kk)*mean(Ue(kk,:)).^3; %% JW ??
    
    %% Derivative cost function with respect to state and beta
    JXXu(xline(kk,:)-2,yline{kk}-1) =  -1/F*pi*(0.5*Drotor)^2*6*Rho*beta(kk).*mean(Ue(kk,:)).^2.*dUedU(kk,1:length(yline{kk}))./(length(yline{kk}));%-3*mean(urotor).^2*pi*(0.5*Drotor)^2*2*beta(kk)/(length(solind)/2);
    % JXXv(xline(kk,:)-1,yline{kk}-2) = -1/F*pi*(0.5*Drotor)^2*6*Rho*beta(kk).*(Ue(kk,:).^2.*dUedv./((length(yline{kk}))));%-3*mean(urotor).^2*pi*(0.5*Drotor)^2*2*beta(kk)/(length(solind)/2);
    JXXv(xline(kk,:)-1,ylinev{kk}-2) = -1/F*pi*(0.5*Drotor)^2*6*Rho*beta(kk).*mean(Ue(kk,:)).^2.*dUedv(kk,:)./((length(yline{kk})))*dvvdv;%-3*mean(urotor).^2*pi*(0.5*Drotor)^2*2*beta(kk)/(length(solind)/2);
    
    %% two lines above should be checked!!
    dSm.dJdx  = [vec(JXXu');vec(JXXv');vec(JXXp')];
    dSm.dJdbeta(kk)   = -1/F*pi*(0.5*Drotor)^2*(2*Rho).*mean(Ue(kk,:)).^3;
    dSm.dJdPhi(kk)  = -1/F*6*Rho*pi*(0.5*Drotor)^2*beta(kk)*mean(Ue(kk,:)).^2.*mean(dUedPhi(kk,:));

    dSmxdbetakk = sparse(Nx-3,Ny-2);
    dSmydbetakk = sparse(Nx-2,Ny-3);
    dSmxdphikk  = sparse(Nx-3,Ny-2);
    dSmydphikk  = sparse(Nx-2,Ny-3);
    
    dSmxdbetakk(xline(kk,:)-2,yline{kk}-1)        =  -(dFxdbeta.*dyy2(1,yline{kk}))';  %JW changed definition of force and location (need correction for yaw
    dSmydbetakk(xline(kk,:)-1,yline{kk}(2:end)-2) = scale*(dFydbeta(2:end).*dyy2(1,yline{kk}(2:end)))'; %% JW made this 0
    dSm.dbeta(:,kk)                               = [vec(dSmxdbetakk');vec(dSmydbetakk')];
    
    dSmxdphikk(xline(kk,:)-2,yline{kk}-1)         = -dFxdPhi'.*dyy2(1,yline{kk})';
    dSmydphikk(xline(kk,:)-1,yline{kk}(2:end)-2)  = scale*dFydPhi(2:end)'.*dyy2(1,yline{kk}(2:end))';
    dSm.dphi(:,kk)                                = [vec(dSmxdphikk');vec(dSmydphikk')];
    
    dSm.dPower_in(:,kk)                           = [vec(dSmxdbetakk');vec(dSmydbetakk')].*dSm.dbetadPower_in(kk);

    
    %% Input to Ax=b
    Sm.x(xline(kk,:)-2,yline{kk}-1)      = -Fx'.*dyy2(1,yline{kk})'*beta(kk)/beta(kk);
    Sm.dx(xline(kk,:)-2,yline{kk}-1)     = -dFx'.*dyy2(1,yline{kk})'*dbeta(kk);   %dSm.dbeta * dbeta
    
    Sm.y(xline(kk,:)-1,yline{kk}(2:end)-2)      = scale*Fy(2:end)'.*dyy2(1,yline{kk}(2:end))'*beta(kk)/beta(kk);
    Sm.dy(xline(kk,:)-1,yline{kk}-2)     = dFy'.*dyy2(1,yline{kk})'*dbeta(kk);
    
    
    dFu                                   = 4*beta(kk)*Rho*u(xline(kk,:),yline{kk})*cos(Phi(kk)*pi/180);
    dFv                                   = 4*beta(kk)*Rho*v(xline(kk,:),yline{kk})*sin(Phi(kk)*pi/180);
    Smdu(xline(kk,:)-2,yline{kk}-1) = -dFu'.*dyy2(1,yline{kk})';
    Smdv(xline(kk,:)-1,yline{kk}-2) = dFv'.*dyy2(1,yline{kk})';
    
    dSm.dx = blkdiag(diag(vec(Smdu')'),diag(vec(Smdv')'),sparse((Ny-2)*(Nx-2),(Ny-2)*(Nx-2)));
    
    %% Input to qLPV
    tempx(xline(kk,:)-2,yline{kk}-1,kk) = -Fx'.*dyy2(1,yline{kk})';
    tempy(xline(kk,:)-1,yline{kk}(2:end)-2,kk) = scale*Fy(2:end)'.*dyy2(1,yline{kk}(2:end))';
    Sm.xx(:,kk)                               = vec(tempx(:,:,kk)');
    Sm.yy(:,kk)                               = vec(tempy(:,:,kk)');
    
    tempdx(xline(kk,:)-2,yline{kk}-1,kk) = -dFx'.*dyy2(1,yline{kk})';
    tempdy(xline(kk,:)-1,yline{kk}-2,kk) = dFy'.*dyy2(1,yline{kk})';
    Sm.dxx(:,kk)                               = vec(tempdx(:,:,kk)');
    Sm.dyy(:,kk)                               = vec(tempdy(:,:,kk)');
    
    %% checks
    Sm.UU=U;
    Sm.UU2=Ue;
    Sm.Uv=v(xline(kk,:),ylinev{kk});
    Sm.Uu=u(xline(kk,:),yline{kk});
    Sm.Power2=Power2;
    
    %% don't really need them
    Umean(kk)     = mean(u(xline(kk,:),yline{kk}));
    Vmean(kk)     = mean(v(xline(kk,:),yline{kk}));
    Vmean(kk)     = mean(vv);
    Umatrix       = sqrt(u(xline(kk,:):xline(kk,:),yline{kk}).^2+v(xline(kk,:):xline(kk,:),yline{kk}).^2);
    Ur(kk,:)      = mean(Ue(kk,:));
    Ueffect(kk)   = Ur(kk)/(1-a(kk));
    
end

Sm.dxx  = [Sm.dxx zeros(size(Sm.dxx,1),size(dPhi,1))];
Sm.dyy  = [Sm.dyy zeros(size(Sm.dyy,1),size(dPhi,1))];
Sm.xx   = [Sm.xx zeros(size(Sm.xx,1),size(Phi,1))];
Sm.yy   = [Sm.yy zeros(size(Sm.yy,1),size(Phi,1))];
output.Sm  = Sm;
output.dSm = dSm;
