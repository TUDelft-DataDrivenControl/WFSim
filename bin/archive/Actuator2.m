function [Sm,dSm,Ueffect,Ur,a,Power,CT,CP,phi,Fxmean,Fymean,Umean,Vmean,Umatrix] = Actuator(F,Wp,u,v,Phi,beta,dbeta,Rho)

Sm.x     = sparse(Wp.Nx-3,Wp.Ny-2);
Sm.dx    = sparse(Wp.Nx-3,Wp.Ny-2);

Sm.xx    = sparse((Wp.Nx-3)*(Wp.Ny-2),Wp.N);
Sm.dxx   = sparse((Wp.Nx-3)*(Wp.Ny-2),Wp.N);
Sm.y     = sparse(Wp.Nx-2,Wp.Ny-3);
Sm.dy    = sparse(Wp.Nx-2,Wp.Ny-3);
tempy    = zeros(Wp.Nx-2,Wp.Ny-3);
tempdy   = zeros(Wp.Nx-2,Wp.Ny-3);
Sm.yy    = sparse((Wp.Nx-2)*(Wp.Ny-3),Wp.N);
Sm.dyy   = sparse((Wp.Nx-2)*(Wp.Ny-3),Wp.N);

Smdu     = sparse(Wp.Nx-3,Wp.Ny-2);
Smdv     = sparse(Wp.Nx-2,Wp.Ny-3);

dSm.xdu     = sparse((Wp.Nx-3)*(Wp.Ny-2),(Wp.Nx-3)*(Wp.Ny-2));
dSm.xdv     = sparse((Wp.Nx-3)*(Wp.Ny-2),(Wp.Nx-2)*(Wp.Ny-3) );
dSm.ydu     = sparse((Wp.Nx-2)*(Wp.Ny-3),(Wp.Nx-3)*(Wp.Ny-2));
dSm.ydv     = sparse((Wp.Nx-2)*(Wp.Ny-3),(Wp.Nx-2)*(Wp.Ny-3));

dSmxdbetakk = sparse(Wp.Nx-3,Wp.Ny-2);
dSmydbetakk = sparse(Wp.Nx-2,Wp.Ny-3);

dSmxdphikk  = sparse(Wp.Nx-3,Wp.Ny-2);
dSmydphikk  = sparse(Wp.Nx-2,Wp.Ny-3);

JXXu  = sparse(Wp.Nx-3,Wp.Ny-2);
JXXv  = sparse(Wp.Nx-2,Wp.Ny-3);
JXXp  = sparse(Wp.Nx-2,Wp.Ny-2);

for kk=1:Wp.N
    %% Oppassen met multiple turbine dSmxdbetakk  dSmxdphikk
    clear dSmxdbetakk  dSmxdphikk dSmydbetakk  dSmydphikk tempx tempdx
    tempx    = zeros(Wp.Nx-3,Wp.Ny-2);
    tempdx   = zeros(Wp.Nx-3,Wp.Ny-2);
    
    dSmxdbetakk = sparse(Wp.Nx-3,Wp.Ny-2);
    dSmydbetakk = sparse(Wp.Nx-2,Wp.Ny-3);
    
    dSmxdphikk  = sparse(Wp.Nx-3,Wp.Ny-2);
    dSmydphikk  = sparse(Wp.Nx-2,Wp.Ny-3);
    
    a(kk)         =     beta(kk)/(beta(kk)+1);                      % Axial induction factor at turbine
    vv            =     diff(v(Wp.xline(kk,:),Wp.ylinev{kk}))/2+v(Wp.xline(kk,:),Wp.ylinev{kk}(1:end-1));
    U(kk,:)       =     sqrt(u(Wp.xline(kk,:),Wp.yline{kk}).^2+vv.^2);
    uu            =     u(Wp.xline(kk,:),Wp.yline{kk});
    
    dvvdv         =     1/2*diag(sign(Wp.ylinev{kk}))+1/2*diag(sign(Wp.ylinev{kk}(1:end-1)),1); dvvdv(end,:)=[];
    dUdv(kk,:)    =     vv./(U(kk,:));
    dUdu(kk,:)    =     u(Wp.xline(kk,:),Wp.yline{kk})./(U(kk,:));
    phi(kk)       =     0;%-180/pi*atan(mean(v(Wp.xline(kk,:)-5,Wp.yline{kk}))/mean(u(Wp.xline(kk,:)-5,Wp.yline{kk}))); %% added an offset here
    Ue(kk,:)      =     cos(atan(vv./(u(Wp.xline(kk,:),Wp.yline{kk})))+Phi(kk)/180*pi).*U(kk,:);
  
    
    
    uu=u(Wp.xline(kk,:),Wp.yline{kk});
    dUedv(kk,:)=-sin(atan(vv./(u(Wp.xline(kk,:),Wp.yline{kk}))) + Phi(kk)/180*pi)./((u(Wp.xline(kk,:),Wp.yline{kk})).*(vv.^2./(u(Wp.xline(kk,:),Wp.yline{kk})).^2 + 1)).*U(kk,:)+cos(atan(vv./(u(Wp.xline(kk,:),Wp.yline{kk})))+Phi(kk)/180*pi).*dUdv(kk,:); %% diff Ue wrt v
    dUedv(kk,:)=-uu.*sin(atan(vv./uu) + Phi(kk)/180*pi)./(vv.^2+uu.^2).*U(kk,:)+cos(atan(vv./(u(Wp.xline(kk,:),Wp.yline{kk})))+Phi(kk)/180*pi).*dUdv(kk,:); %% diff Ue wrt v
    dUedu(kk,:)=(vv.*sin(atan(vv./(u(Wp.xline(kk,:),Wp.yline{kk}))) + Phi(kk)/180*pi))./((vv.^2 + (u(Wp.xline(kk,:),Wp.yline{kk})).^2)).*U(kk,:)+cos(atan(vv./(u(Wp.xline(kk,:),Wp.yline{kk})))+Phi(kk)/180*pi).*dUdu(kk,:); %% diff Ue wrt u
    dUedPhi(kk,:)=-1/180*pi*sin(atan(vv./(u(Wp.xline(kk,:),Wp.yline{kk}))) + Phi(kk)/180*pi).*U(kk,:); %% diff wrt to Ue
    
        % http://www.nrel.gov/docs/fy05osti/36834.pdf
    %F         = 1.38; THIS NOW BECOMES A TUNING PARAMETER
    
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
    Fthrust             = 1/2*Rho*CT(kk)*((beta(kk)+1)).^2.*Ue(kk,:).*Ue(kk,:); % similar as 1/2*Rho*CT(kk)*(U/(1-a(kk))).^2; %somethin to try 2*60*ones(1,size(U(kk,:),2));%
    dFthrustdbeta       = 1/2*Rho*CT(kk)*Ue(kk,:).^2*(2*beta(kk)+2)+1/2*Rho*dCTdbeta(kk)*Ue(kk,:).^2*(beta(kk)+1).^2; %diff to beta 1/2*Rho*CT(kk)*U^2*(beta(kk)+1).^2;
    dFxdbeta            = dFthrustdbeta*cos((Phi(kk)+phi(kk))*pi/180);
    dFydbeta            = dFthrustdbeta*sin((Phi(kk)+phi(kk))*pi/180);
    
    dSmxdbetakk(Wp.xline(kk,:)-2,Wp.yline{kk}-1) =  -(dFxdbeta.*Wp.dyy2(1,Wp.yline{kk}))';  %JW changed definition of force and location (need correction for yaw
    scale=2;
    dSmydbetakk(Wp.xline(kk,:)-1,Wp.yline{kk}(2:end)-2) = scale*(dFydbeta(2:end).*Wp.dyy2(1,Wp.yline{kk}(2:end)))'; %% JW made this 0
    
    dSm.dbeta(:,kk)   = [vec(dSmxdbetakk');vec(dSmydbetakk')];
    
    %% Gradient of the thrust with respect to yaw angle [to be checked]
    dFthrustdPhi                = Rho*CT(kk)*Ue(kk,:).*((beta(kk)+1)).^2.*dUedPhi(kk,:);%%1/2*Rho*dCTdPhi(kk)*(Ue(kk,:)*(beta(kk)+1)).^2;
    dFxdPhi                     = dFthrustdPhi*cos((Phi(kk)+phi(kk))*pi/180)-pi/180*Fthrust*sin((Phi(kk)+phi(kk))*pi/180);
    dFydPhi                     = +dFthrustdPhi*sin((Phi(kk)+phi(kk))*pi/180)+pi/180*Fthrust*cos((Phi(kk)+phi(kk))*pi/180);
    
    dSmxdphikk(Wp.xline(kk,:)-2,Wp.yline{kk}-1) = -dFxdPhi'.*Wp.dyy2(1,Wp.yline{kk})';
    dSmydphikk(Wp.xline(kk,:)-1,Wp.yline{kk}(2:end)-2) = scale*dFydPhi(2:end)'.*Wp.dyy2(1,Wp.yline{kk}(2:end))';
    
    dSm.dphi(:,kk)   = [vec(dSmxdphikk');vec(dSmydphikk')];
    %% Gradient of the Thrust with respect to the state
    dFxdU = Rho*CT(kk)*Ue(kk,:)*((beta(kk)+1)).^2*cos((phi(kk)+Phi(kk))*pi/180);
    dFydU = Rho*CT(kk)*Ue(kk,:)*((beta(kk)+1)).^2*sin((phi(kk)+Phi(kk))*pi/180);
    dSm.xdu((Wp.xline(kk,:)-3)*(Wp.Ny-2)+Wp.yline{kk}-1,(Wp.xline(kk,:)-3)*(Wp.Ny-2)+Wp.yline{kk}-1) = diag((-dFxdU.*dUedu(kk,:).*Wp.dyy2(1,Wp.yline{kk})));
    %dSm.xdv((Wp.xline(kk,:)-3)*(Wp.Ny-2)+Wp.ylinev{kk}-1,(Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.ylinev{kk}-2) = diag((-dFxdU.*dUedv(kk,:).*Wp.dyy2(1,Wp.yline{kk})*dvvdv));
    dSm.xdv((Wp.xline(kk,:)-3)*(Wp.Ny-2)+Wp.yline{kk}-1,(Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}-2) = 1/2*diag((-dFxdU.*dUedv(kk,:).*Wp.dyy2(1,Wp.yline{kk})));
    dSm.xdv((Wp.xline(kk,:)-3)*(Wp.Ny-2)+Wp.yline{kk}-1,(Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}-1) = dSm.xdv((Wp.xline(kk,:)-3)*(Wp.Ny-2)+Wp.yline{kk}-1,(Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}-1)+ 1/2*diag((-dFxdU.*dUedv(kk,:).*Wp.dyy2(1,Wp.yline{kk})));
    dSm.ydu((Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}(2:end)-2,(Wp.xline(kk,:)-3)*(Wp.Ny-2)+Wp.yline{kk}(2:end)-1) = diag(scale*(dFydU(2:end).*dUedu(kk,2:end).*Wp.dyy2(1,Wp.yline{kk}(2:end))));
   %dSm.ydv((Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}-2,(Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}-2) = diag(scale*(dFydU(2:end).*dUedv(kk,2:end).*Wp.dyy2(1,Wp.yline{kk}(2:end))*dvvdv(2:end,2:end)));
   dSm.ydv((Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}(2:end)-2,(Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}(2:end)-2) = 1/2*diag(scale*(dFydU(2:end).*dUedv(kk,2:end).*Wp.dyy2(1,Wp.yline{kk}(2:end))));
   dSm.ydv((Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}(2:end)-2,(Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}(2:end)-1) =   dSm.ydv((Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}(2:end)-2,(Wp.xline(kk,:)-2)*(Wp.Ny-3)+Wp.yline{kk}(2:end)-1)+ 1/2*diag(scale*(dFydU(2:end).*dUedv(kk,2:end).*Wp.dyy2(1,Wp.yline{kk}(2:end))));
   
   
    %     dSm.xdu(Wp.xline(kk,:)-2,Wp.yline{kk}-1)    = (-dFxdU.*dUedu.*Wp.dyy2(1,Wp.yline{kk}));
    %     dSm.xdv(Wp.xline(kk,:)-1,Wp.ylinev{kk}-2)    = (-dFxdU.*dUedv.*Wp.dyy2(1,Wp.yline{kk})*dvvdv)';
    %     dSm.ydu(Wp.xline(kk,:)-2,Wp.yline{kk}(2:end)-1)    =scale*(dFydU(2:end).*dUedu(2:end).*Wp.dyy2(1,Wp.yline{kk}(2:end)))';
    %     dSm.ydv(Wp.xline(kk,:)-1,Wp.yline{kk}(1:end)-2)    =scale*(dFydU(2:end).*dUedv(2:end).*Wp.dyy2(1,Wp.yline{kk}(2:end))*dvvdv(2:end,2:end))';
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
    
    Fxmean(kk)    = mean(-Fx'.*Wp.dyy2(1,Wp.yline{kk})');
    Fymean(kk)    = mean(Fy'.*Wp.dyy2(1,Wp.yline{kk})');
    
    %%
    CP(kk)        = 4*beta(kk)/(beta(kk)+1)^3*0.768*cos((Phi(kk)-phi(kk))*pi/180)^(1.88);
    Power2(kk)     =1/F* Wp.powerscale*mean(2*Rho*pi*(0.5*Wp.turbine.Drotor)^2*beta(kk)*Ue(kk,:).^3); %% JW ??
    Power(kk)     = 1/F* Wp.powerscale*2*Rho*pi*(0.5*Wp.turbine.Drotor)^2*beta(kk)*mean(Ue(kk,:)).^3; %% JW ??
    %% still have to change the linearisation for the new power definition
    %% Derivative cost function with respect to state and beta
     JXXu(Wp.xline(kk,:)-2,Wp.yline{kk}-1) =  -1/F* Wp.powerscale*pi*(0.5*Wp.turbine.Drotor)^2*6*Rho*beta(kk).*mean(Ue(kk,:)).^2.*dUedu(kk,:)./(length(Wp.yline{kk}));%-3*mean(urotor).^2*pi*(0.5*Wp.turbine.Drotor)^2*2*beta(kk)/(length(solind)/2);
    % JXXv(Wp.xline(kk,:)-1,Wp.yline{kk}-2) = -1/F*pi*(0.5*Wp.turbine.Drotor)^2*6*Rho*beta(kk).*(Ue(kk,:).^2.*dUedv./((length(Wp.yline{kk}))));%-3*mean(urotor).^2*pi*(0.5*Wp.turbine.Drotor)^2*2*beta(kk)/(length(solind)/2);
    JXXv(Wp.xline(kk,:)-1,Wp.ylinev{kk}-2) = -1/F* Wp.powerscale*pi*(0.5*Wp.turbine.Drotor)^2*6*Rho*beta(kk).*mean(Ue(kk,:)).^2.*dUedv(kk,:)./((length(Wp.yline{kk})))*dvvdv;%-3*mean(urotor).^2*pi*(0.5*Wp.turbine.Drotor)^2*2*beta(kk)/(length(solind)/2);
    
    %% two lines above should be checked!!
    
    dSm.dJdx  = [vec(JXXu');vec(JXXv');vec(JXXp')];
    dSm.dJdbeta(kk)   = -1/F* Wp.powerscale*pi*(0.5*Wp.turbine.Drotor)^2*(2*Rho).*mean(Ue(kk,:)).^3;
    dSm.dJdPhi(kk)  = -1/F* Wp.powerscale*6*Rho*pi*(0.5*Wp.turbine.Drotor)^2*beta(kk)*mean(Ue(kk,:)).^2.*mean(dUedPhi(kk,:));
    
    
    %% Input to Ax=b
    Sm.x(Wp.xline(kk,:)-2,Wp.yline{kk}-1)      = -Fx'.*Wp.dyy2(1,Wp.yline{kk})';
    Sm.dx(Wp.xline(kk,:)-2,Wp.yline{kk}-1)     = -dFx'.*Wp.dyy2(1,Wp.yline{kk})'*dbeta(kk);   %dSm.dbeta * dbeta
    
    Sm.y(Wp.xline(kk,:),Wp.yline{kk}(2:end)-2)      = scale*Fy(2:end)'.*Wp.dyy2(1,Wp.yline{kk}(2:end))';

    Sm.y(Wp.xline(kk,:)-2,Wp.yline{kk}(2:end)-2)      = scale*Fy(2:end)'.*Wp.dyy2(1,Wp.yline{kk}(2:end))';
 Sm.y(Wp.xline(kk,:)-1,Wp.yline{kk}(2:end)-2)      = scale*Fy(2:end)'.*Wp.dyy2(1,Wp.yline{kk}(2:end))';

    Sm.dy(Wp.xline(kk,:)-1,Wp.yline{kk}-2)     = dFy'.*Wp.dyy2(1,Wp.yline{kk})'*dbeta(kk);
    
    
    dFu                                   = 4*beta(kk)*Rho*u(Wp.xline(kk,:),Wp.yline{kk})*cos(Phi(kk)*pi/180);
    dFv                                   = 4*beta(kk)*Rho*v(Wp.xline(kk,:),Wp.yline{kk})*sin(Phi(kk)*pi/180);
    Smdu(Wp.xline(kk,:)-2,Wp.yline{kk}-1) = -dFu'.*Wp.dyy2(1,Wp.yline{kk})';
    Smdv(Wp.xline(kk,:)-1,Wp.yline{kk}-2) = dFv'.*Wp.dyy2(1,Wp.yline{kk})';
    
    dSm.dx = blkdiag(diag(vec(Smdu')'),diag(vec(Smdv')'),sparse((Wp.Ny-2)*(Wp.Nx-2),(Wp.Ny-2)*(Wp.Nx-2)));
    
    %% Input to qLPV
    tempx(Wp.xline(kk,:)-2,Wp.yline{kk}-1,kk) = -Fx'.*Wp.dyy2(1,Wp.yline{kk})';
    tempy(Wp.xline(kk,:)-1,Wp.yline{kk}(2:end)-2,kk) = scale*Fy(2:end)'.*Wp.dyy2(1,Wp.yline{kk}(2:end))';
    Sm.xx(:,kk)                               = vec(tempx(:,:,kk)')/beta(kk);
    Sm.yy(:,kk)                               = vec(tempy(:,:,kk)')/beta(kk);
    
    tempdx(Wp.xline(kk,:)-2,Wp.yline{kk}-1,kk) = -dFx'.*Wp.dyy2(1,Wp.yline{kk})';
    tempdy(Wp.xline(kk,:)-1,Wp.yline{kk}-2,kk) = dFy'.*Wp.dyy2(1,Wp.yline{kk})';
    Sm.dxx(:,kk)                               = vec(tempdx(:,:,kk)');
    Sm.dyy(:,kk)                               = vec(tempdy(:,:,kk)');
    
    %% checks
    Sm.UU=U;
    Sm.UU2=Ue;
    Sm.Uv=v(Wp.xline(kk,:),Wp.ylinev{kk});
    Sm.Uu=u(Wp.xline(kk,:),Wp.yline{kk});
    Sm.Power2=Power2;
    
    %% don't really need them
    Umean(kk)     = mean(u(Wp.xline(kk,:),Wp.yline{kk}));
    Vmean(kk)     = mean(v(Wp.xline(kk,:),Wp.yline{kk}));
    Vmean(kk)     = mean(vv);
    Umatrix       = sqrt(u(Wp.xline(kk,:):Wp.xline(kk,:),Wp.yline{kk}).^2+v(Wp.xline(kk,:):Wp.xline(kk,:),Wp.yline{kk}).^2);
    Ur(kk,:)        = mean(Ue(kk,:));
    Ueffect(kk)   = Ur(kk)/(1-a(kk));
    
end

