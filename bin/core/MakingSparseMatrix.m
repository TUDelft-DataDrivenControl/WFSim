function Ax = MakingSparseMatrix(Nx,Ny,ax,ix,iy,q)
% ix is the index where i begins
% iy is the index where j begins

% Add central components
Ax      =   spdiags(vec(ax.aP(ix:end-q,iy:end-q)'),0,(Nx-ix-1+q)*(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q));

% Add north components
ann     =   vec([zeros(Nx-ix-1+q,1)   ax.aN(ix:end-q,iy:end-q-1) ]') ;
Ax      =   -spdiags(ann,1,(Nx-ix-1+q)*(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q))+Ax;

% Add south components
ass     =   vec([ax.aS(ix:end-q,iy+1:end-q) zeros(Nx-ix-1+q,1) ]') ;
Ax      =   -spdiags(ass,-1,(Nx-ix-1+q)*(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q))+Ax;

% Add east components
aee     =   vec([ zeros(1,Ny-iy-1+q);ax.aE(ix:end-q-1,iy:end-q);]');
Ax      =   -spdiags(aee,Ny-iy-1+q,(Nx-ix-1+q)*(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q))+Ax;

% Add west components
aww     =   vec([ax.aW(ix+1:end-q,iy:end-q);zeros(1,Ny-iy-1+q)]');
Ax      =   -spdiags(aww(1:end),-(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q))+Ax;

%Aw      = -spdiags(aww(1:end),-(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q)); 






