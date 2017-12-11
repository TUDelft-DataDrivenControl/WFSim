function lmu = ConstructLmu( Wp )

Nx      = Wp.mesh.Nx;
Ny      = Wp.mesh.Ny;
xline   = Wp.mesh.xline;
yline   = Wp.mesh.yline;
d       = Wp.site.m+15; % 0< d < xline(i+1)- xline(i)
dprime  = Wp.site.n;
max_lmu = Wp.site.lmu; % = ls*(d-dprime)

uniqueXline = sort(unique(xline));
for kk=1:length(yline);uniqueYline(kk,:)=yline{kk};end
uniqueYline = unique(sort(uniqueYline),'rows');

% lmu in x-direction
lmux    = [];
for kk=1:length(uniqueXline)-1
    a                        = uniqueXline(kk+1)-uniqueXline(kk);
    temp                     = zeros(1,a);
    temp(dprime+1:d-dprime)  = linspace(0,max_lmu,d-2*dprime);
    lmux                     = [lmux temp];
end

lmux = [zeros(1,uniqueXline(1)) lmux zeros(1,dprime) linspace(0,max_lmu,Nx-uniqueXline(end)-dprime)];

% lmu in y-direction
lmuy    = zeros(1,Ny);
for kk=1:size(uniqueYline,1)
    lmuy(uniqueYline(kk,:)) = 1;
end

% lmu in x- and y-direction
lmu  = (repmat(lmux',1,Ny).*repmat(lmuy,Nx,1))*0.5*Wp.turbine.Drotor; 

end

