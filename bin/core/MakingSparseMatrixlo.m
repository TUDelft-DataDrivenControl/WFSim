function [Ax_off,Ay_off] = MakingSparseMatrixlo(Nx,Ny,dax,day)


Ax_off = sparse((Ny-2)*(Nx-3),(Nx-2)*(Ny-3));
Ay_off = sparse((Ny-3)*(Nx-2),(Ny-2)*(Nx-3));

Bx_W   = sparse(Ny-2,Ny-3); Bx_E   = sparse(Ny-2,Ny-3);
By_W   = sparse(Ny-3,Ny-2); By_E   = sparse(Ny-3,Ny-2);


for x= 3:Nx-1;
    swx      = -spdiags(dax.SW(x,3:Ny-1)',0,Ny-3,Ny-3);
    nwx      = -spdiags(dax.NW(x,2:Ny-2)',0,Ny-3,Ny-3);
    
    Bx_W     = [sparse(1,Ny-3); swx] + [nwx; sparse(1,Ny-3)];
    
    sex      = -spdiags(dax.SE(x,3:Ny-1)',0,Ny-3,Ny-3);
    nex      = -spdiags(dax.NE(x,2:Ny-2)',0,Ny-3,Ny-3);
    Bx_E     = [sparse(1,Ny-3); sex] + [nex; sparse(1,Ny-3)];
    
    Ax_off((x-3)*(Ny-2)+1:(x-2)*(Ny-2) ,(x-3)*(Ny-3)+1:(x-1)*(Ny-3)) = [Bx_W Bx_E];
end

for y= 3:Nx-1;
    swy      = -spdiags(day.SW(y,3:Ny-1)',0,Ny-3,Ny-3);
    nwy      = -spdiags(day.NW(y,3:Ny-1)',0,Ny-3,Ny-3); % changed JW %nwy      = -spdiags(day.NW(y,3:Ny-2)',0,Ny-3,Ny-3);
    By_W     = [swy sparse(Ny-3,1)] + [sparse(Ny-3,1) nwy];
    
    sey      = -spdiags(day.SE(y-1,3:Ny-1)',0,Ny-3,Ny-3);
    ney      = -spdiags(day.NE(y-1,3:Ny-1)',0,Ny-3,Ny-3); % changed JW %ney      = -spdiags(day.NE(y-1,2:Ny-2)',0,Ny-3,Ny-3);
    By_E     = [sey sparse(Ny-3,1)] + [sparse(Ny-3,1) ney];
    
    Ay_off((y-3)*(Ny-3)+1:(y-1)*(Ny-3) ,(y-3)*(Ny-2)+1:(y-2)*(Ny-2)) = [By_E; By_W];
end    