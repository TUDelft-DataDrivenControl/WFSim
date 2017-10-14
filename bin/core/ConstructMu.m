function mu = ConstructMu(Wp)

m      = Wp.site.m;
n      = Wp.site.n;
xline  = Wp.mesh.xline;
yline  = Wp.mesh.yline;
Nx     = Wp.mesh.Nx;
Ny     = Wp.mesh.Ny;
lmu    = Wp.site.lmu;
Drotor = Wp.turbine.Drotor;

if mod(m,2); m = m+1; end

if strcmp(Wp.name,'6turb_adm_partial')
    xxline            = sort(unique(xline));
    xm1               = [zeros(1,xxline(1)+n) linspace(0,lmu,xxline(2)-xxline(1)-m)]';
    xm2               = [zeros(m,1);linspace(0,lmu,xxline(3)-xxline(2)-m)'];
    xm3               = [zeros(m,1);linspace(0,lmu,Nx-xxline(3)-n)'];
    
    xm                = [xm1;xm2;xm3];
    
    ym1               = [zeros(1,yline{1}(1)-1) ones(1,length(yline{1})) zeros(1,yline{2}(1)-yline{1}(end)-1)...
        ones(1,length(yline{2})) zeros(1,Ny-yline{2}(end))] ;
    ym2               = [zeros(1,yline{3}(1)-1) ones(1,length(yline{3})) zeros(1,yline{4}(1)-yline{3}(end)-1)...
        ones(1,length(yline{4})) zeros(1,Ny-yline{4}(end))] ;
    ym3               = ym1;
    
    ml1               = (repmat(xm1,1,length(ym1)).*repmat(ym1,length(xm1),1))*0.5*Drotor;
    ml2               = (repmat(xm2,1,length(ym2)).*repmat(ym2,length(xm2),1))*0.5*Drotor;
    ml3               = (repmat(xm3,1,length(ym3)).*repmat(ym3,length(xm3),1))*0.5*Drotor;
    mu                = [ml1;ml2;ml3];
    
else
    
    xxline = sort(unique(xline));
    mux    = [];
    
    for kk=1:length(xxline)-1
        a                    = xxline(kk+1)-xxline(kk);
        temp                 = zeros(1,a);
        temp(m/2:end-m/2-1)  = linspace(0,lmu,a-m);
        mux                  = [mux temp];
    end
    
    mux = [zeros(1,xxline(1)) mux zeros(1,m/2-1) linspace(0,lmu,Nx-xxline(end)-m/2+1)];
    
    for kk=1:length(yline)
        yyline(kk,:) = yline{kk};
    end
    yyline = unique(sort(yyline),'rows');
    muy    = zeros(1,Ny);
    
    for kk=1:size(yyline,1)
        muy(yyline(kk,:)) = 1;
    end
    
    mu  = (repmat(mux(1:Nx)',1,Ny).*repmat(muy,Nx,1))*0.5*Drotor; % Dynamic flow viscosity
    
end