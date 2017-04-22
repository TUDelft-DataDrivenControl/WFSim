function mu = ConstructMu(Wp)

m      = Wp.site.m;
xline  = Wp.mesh.xline;
yline  = Wp.mesh.yline;
Nx     = Wp.mesh.Nx;
Ny     = Wp.mesh.Ny;
lmu    = Wp.site.lmu;
Drotor = Wp.turbine.Drotor;

if mod(m,2); m = m+1; end

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