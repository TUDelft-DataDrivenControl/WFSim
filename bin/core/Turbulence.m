Turbulencemodel = 'WFSim4';

% Define spatial varying mixing-length parameter
if N==1
    x                = [zeros(1,xline(1)+n) linspace(0,lmu,Nx-xline(1)-n)]';
    y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,Ny-yline{1}(end))] ;
    mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
elseif N==2
    %x                = [zeros(1,xline(1)+n) linspace(0,lmu,xline(1+1)-xline(1)-m)]';
    %x                = [x;zeros(m,1);linspace(0,lmu,Nx-xline(2)-n)'];
    x                = [zeros(1,xline(1)+n) linspace(0,lmu,xline(1+1)-xline(1)-4*n)]';
    x                = [x;zeros(4*n,1);linspace(0,lmu,Nx-xline(2)-n)'];
    y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,Ny-yline{1}(end))] ;
    %y                = [zeros(1,yline{1}(1)-2) ones(1,length(yline{1})+2) zeros(1,Ny-yline{1}(end)-1)] ;
    mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
elseif N==3 || N==6
    xline            = sort(unique(xline));
    x                = [zeros(1,xline(1)+n) linspace(0,lmu,xline(2)-xline(1)-m)]';
    x                = [x;zeros(m,1);linspace(0,lmu,xline(3)-xline(2)-m)'];
    x                = [x;zeros(m,1);linspace(0,lmu,Nx-xline(3)-n)'];
    y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,Ny-yline{1}(end))] ;
    %mixing_length    = repmat(x,1,Ny)*0.5*Drotor;
    mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
elseif N==9
    xline            = sort(unique(xline));
    x                = [zeros(1,xline(1)+n) linspace(0,lmu,xline(2)-xline(1)-m)]';
    x                = [x;zeros(m,1);linspace(0,lmu,xline(3)-xline(2)-m)'];
    x                = [x;zeros(m,1);linspace(0,lmu,Nx-xline(3)-n)'];
    y                = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,yline{4}(1)-yline{1}(end)-1) ...
        ones(1,length(yline{4})) zeros(1,yline{3}(1)-yline{4}(end)-1) ...
        ones(1,length(yline{3})) zeros(1,Ny-yline{3}(end))];
    mixing_length    = (repmat(x,1,Ny).*repmat(y,Nx,1))*0.5*Drotor;
else
    mixing_length    = lmu*0.5*Drotor;
end


switch lower(Turbulencemodel)

    case lower('WFSim4')

        % For u-momentum equation
        ax.Tex              = zeros(Nx,Ny);
        ax.Twx              = zeros(Nx,Ny);
        ax.Tnx              = zeros(Nx,Ny);
        ax.Tsx              = zeros(Nx,Ny);
        
        ax.Tnex             = zeros(Nx,Ny);
        ax.Tnwx             = zeros(Nx,Ny);
        ax.Tsex             = zeros(Nx,Ny);
        ax.Tswx             = zeros(Nx,Ny);
        
        ax.Tex(1:Nx-1,1:Ny-1)= Rho*(mixing_length(1:Nx-1,1:Ny-1).^2).*(dyy2(1:Nx-1,1:Ny-1)./(dyy(1:Nx-1,1:Ny-1).*dxx2(1:Nx-1,1:Ny-1))).*abs(u(2:Nx,1:Ny-1)-u(1:Nx-1,1:Ny-1));
        ax.Twx(2:Nx,1:Ny-1)  = Rho*(mixing_length(2:Nx,1:Ny-1).^2).*(dyy2(1:Nx-1,1:Ny-1)./(dyy(1:Nx-1,2:Ny).*dxx2(1:Nx-1,1:Ny-1))).*abs(u(1:Nx-1,1:Ny-1)-u(2:Nx,1:Ny-1));
  
        ax.aE             = ax.aE + ax.Tex;
        ax.aW             = ax.aW + ax.Twx;
        ax.aP             = ax.aP + ax.Twx + ax.Tex;
               
        % For v-momentum equation
        ay.Tey   = zeros(Nx,Ny);
        ay.Twy   = zeros(Nx,Ny);
        ay.Tny   = zeros(Nx,Ny);
        ay.Tsy   = zeros(Nx,Ny);
        
        ay.Tney  = zeros(Nx,Ny);
        ay.Tnwy  = zeros(Nx,Ny);
        ay.Tsey  = zeros(Nx,Ny);
        ay.Tswy  = zeros(Nx,Ny);
        
        ay.Tny(2:Nx,1:Ny-1) = Rho*(mixing_length(2:Nx,1:Ny-1).^2).*(dxx2(1:Nx-1,1:Ny-1)./(dyy(1:Nx-1,1:Ny-1).*dyy2(1:Nx-1,1:Ny-1))).*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1)); 
        ay.Tsy(1:Nx-1,2:Ny) = Rho*(mixing_length(2:Nx,2:Ny).^2).*(dxx2(1:Nx-1,1:Ny-1)./(dyy(1:Nx-1,1:Ny-1).*dyy2(1:Nx-1,1:Ny-1))).*abs(u(1:Nx-1,2:Ny)-u(1:Nx-1,1:Ny-1));
                 
        ay.aN             = ay.aN + ay.Tny;
        ay.aS             = ay.aS + ay.Tsy;
        ay.aP             = ay.aP + ay.Tsy + ay.Tny;
            
    case lower('WFSim3')
    
        % Equal to WFSim 1
        
        
    case lower('WFSim2')
        
        % For u-momentum equation
        ax.Tnx           = zeros(Nx,Ny);
        ax.Tsx           = zeros(Nx,Ny);
        
        ax.Tnx(1:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2)).*abs(u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1));
        ax.Tsx(1:Nx,2:Ny)   = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx(1:Nx,2:Ny)./(dyy(1:Nx,2:Ny).^2)).*abs(u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny));
        
        ax.aN             = ax.aN + ax.Tnx;
        ax.aS             = ax.aS + ax.Tsx;
        ax.aP             = ax.aP + ax.Tnx + ax.Tsx;
        
        % For v-momentum equation
        ay.Tey = zeros(Nx,Ny);
        ay.Twy = zeros(Nx,Ny);
        
        ay.Tey(1:Nx-1,2:Ny)   = Rho*(mixing_length(1:Nx-1,2:Ny).^2)./dyy(1:Nx-1,2:Ny).*abs(u(2:Nx,2:Ny)-u(2:Nx,1:Ny-1)); 
        ay.Twy(2:Nx,2:Ny)     = Rho*(mixing_length(2:Nx,2:Ny).^2)./dyy(1:Nx-1,2:Ny).*abs(u(1:Nx-1,2:Ny)-u(1:Nx-1,1:Ny-1)); 
        
        % Define here Ayo
        Ayo    = sparse((Ny-3)*(Nx-2),(Ny-2)*(Nx-3));
        By_W   = sparse(Ny-3,Ny-2); By_E   = sparse(Ny-3,Ny-2);
        
        for y= 3:Nx-1;
            swy      = spdiags(ay.Twy(y,3:Ny-1)',0,Ny-3,Ny-3); %u_{i,J-1}
            nwy      = -spdiags(ay.Twy(y,3:Ny-1)',0,Ny-3,Ny-3); %u_{i,J}
            By_W     = [swy sparse(Ny-3,1)] + [sparse(Ny-3,1) nwy];
            
            sey      = -spdiags(ay.Tey(y-1,3:Ny-1)',0,Ny-3,Ny-3);%u_{i+1,J-1}
            ney      = spdiags(ay.Tey(y-1,3:Ny-1)',0,Ny-3,Ny-3);%u_{i+1,J}
            By_E     = [sey sparse(Ny-3,1)] + [sparse(Ny-3,1) ney];
            
            Ayo((y-3)*(Ny-3)+1:(y-1)*(Ny-3) ,(y-3)*(Ny-2)+1:(y-2)*(Ny-2)) = [By_E; By_W];
        end
        
        output.Ayo = Ayo;
                
    case lower('WFSim1')
        
        % For u-momentum equation
        ax.Tnx              = zeros(Nx,Ny);
        ax.Tsx              = zeros(Nx,Ny);

        ax.Tnx(1:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2)).*abs(u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1));
        ax.Tsx(1:Nx,2:Ny)   = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx(1:Nx,2:Ny)./(dyy(1:Nx,2:Ny).^2)).*abs(u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny));
        
        %ax.Tnx(1:Nx,2:Ny-1) = Rho*(mixing_length(1:Nx,2:Ny-1).^2).*(dxx(1:Nx,2:Ny-1)./(dyy(1:Nx,2:Ny-1).^2)).*abs(u(1:Nx,3:Ny)-u(1:Nx,2:Ny-1));
        %ax.Tsx(1:Nx,1:Ny-2) = Rho*(mixing_length(1:Nx,1:Ny-2).^2).*(dxx(1:Nx,2:Ny-1)./(dyy(1:Nx,1:Ny-2).^2)).*abs(u(1:Nx,2:Ny-1)-u(1:Nx,1:Ny-2));
        
        ax.aN             = ax.aN + ax.Tnx;
        ax.aS             = ax.aS + ax.Tsx;
        ax.aP             = ax.aP + ax.Tnx + ax.Tsx;
        
        % For v-momentum equation
        ay.Tey            = zeros(Nx,Ny);
        ay.Twy            = zeros(Nx,Ny);
                
        ay.Tey(1:Nx-1,1:Ny) = Rho*(mixing_length(1:Nx-1,1:Ny).^2).*(dyy(1:Nx-1,1:Ny)./(dxx(1:Nx-1,1:Ny).^2)).*abs(v(2:Nx,1:Ny)-v(1:Nx-1,1:Ny));
        ay.Twy(2:Nx,1:Ny)   = Rho*(mixing_length(2:Nx,1:Ny).^2).*(dyy(2:Nx,1:Ny)./(dxx(2:Nx,1:Ny).^2)).*abs(v(1:Nx-1,1:Ny)-v(2:Nx,1:Ny));
                 
        ay.aE             = ay.aE + ay.Tey;
        ay.aW             = ay.aW + ay.Twy;
        ay.aP             = ay.aP + ay.Tey + ay.Twy;
        
        if Linearversion
                       
            % For u-momentum equation
            %dTsxd1           = zeros(Nx,Ny);
            %dTsxd2           = zeros(Nx,Ny);
            %dTnxd1           = zeros(Nx,Ny);
            %dTnxd2           = zeros(Nx,Ny);
            
            % dTsx/du_(i,J)
            %dTsxd1(1:Nx,2:Ny)   = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx2(1:Nx,2:Ny)./(dyy(1:Nx,2:Ny).^2)).*sign((u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny)));
            % dTsx/du_(i,J-1)
            %dTsxd2(1:Nx,2:Ny)   = Rho*(mixing_length(1:Nx,2:Ny).^2).*(dxx2(1:Nx,2:Ny)./(dyy(1:Nx,2:Ny).^2)).*-sign((u(1:Nx,1:Ny-1)-u(1:Nx,2:Ny)));
            % dTnx/du_(i,J)
            %dTnxd1(1:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx2(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2)).*-sign((u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1)));
            % dTnx/du_(i,J+1)
            %dTnxd2(1:Nx,1:Ny-1) = Rho*(mixing_length(1:Nx,1:Ny-1).^2).*(dxx2(1:Nx,1:Ny-1)./(dyy(1:Nx,2:Ny).^2)).*sign((u(1:Nx,2:Ny)-u(1:Nx,1:Ny-1)));
            
            dax.S(1:Nx,2:Ny)   = dax.S(1:Nx,2:Ny)   + ax.Tsx(1:Nx,2:Ny) ;
            dax.N(1:Nx,1:Ny-1) = dax.N(1:Nx,1:Ny-1) + ax.Tnx(1:Nx,1:Ny-1);
            dax.P(1:Nx,1:Ny-1) = dax.P(1:Nx,1:Ny-1) + ax.Tnx(1:Nx,1:Ny-1) + ax.Tsx(1:Nx,1:Ny-1) ;
            
            %dax.S(1:Nx,2:Ny)   = dax.S(1:Nx,2:Ny)+dTsxd2(1:Nx,2:Ny).*u(1:Nx,2:Ny)-dTsxd2(1:Nx,2:Ny).*u(1:Nx,1:Ny-1);
            %dax.N(1:Nx,1:Ny-1) = dax.N(1:Nx,1:Ny-1)-dTnxd2(1:Nx,1:Ny-1).*u(1:Nx,1:Ny-1)+dTnxd2(1:Nx,1:Ny-1).*u(1:Nx,2:Ny);
            %dax.P(1:Nx,1:Ny-1) = dax.P(1:Nx,1:Ny-1)+dTnxd1(1:Nx,1:Ny-1).*u(1:Nx,1:Ny-1)+dTsxd1(1:Nx,1:Ny-1).*u(1:Nx,1:Ny-1) ...
            %-dTnxd1(1:Nx,1:Ny-1).*u(1:Nx,2:Ny)-dTsxd1(1:Nx,1:Ny-1).*u(1:Nx,2:Ny);
            
            % For v-momentum equation
            day.E(1:Nx,2:Ny)   = day.E(1:Nx,2:Ny)   + ay.Tey(1:Nx,2:Ny) ;
            day.W(1:Nx,1:Ny-1) = day.W(1:Nx,1:Ny-1) + ay.Twy(1:Nx,1:Ny-1);
            day.P(1:Nx,1:Ny-1) = day.P(1:Nx,1:Ny-1) + ay.Tey(1:Nx,1:Ny-1) + ay.Twy(1:Nx,1:Ny-1) ;
        end;
        
end
