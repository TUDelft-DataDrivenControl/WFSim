function Ax_on = MakingSparseMatrixl(Nx,Ny,dax,ix,iy,q)
% ix is the index where i begins
% iy is the index where j begins

% Add central components
Ax_on   =   spdiags(vec(dax.P(ix:end-q,iy:end-q)'),0,(Nx-ix-1+q)*(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q));

% Add north components
ann     =   vec([zeros(Nx-ix-1+q,1)   dax.N(ix:end-q,iy:end-q-1) ]') ;
Ax_on   =   -spdiags(ann,1,(Nx-ix-1+q)*(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q))+Ax_on;

% Add south components
ass     =   vec([dax.S(ix:end-q,iy+1:end-q) zeros(Nx-ix-1+q,1) ]') ;
Ax_on   =   -spdiags(ass,-1,(Nx-ix-1+q)*(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q))+Ax_on;

% Add east components
aee     =   vec([ zeros(1,Ny-iy-1+q);dax.E(ix:end-q-1,iy:end-q);]');
Ax_on   =   -spdiags(aee,Ny-iy-1+q,(Nx-ix-1+q)*(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q))+Ax_on;

% Add west components
aww     =   vec([dax.W(ix+1:end-q,iy:end-q);zeros(1,Ny-iy-1+q)]');
Ax_on   =   -spdiags(aww(1:end),-(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q),(Nx-ix-1+q)*(Ny-iy-1+q))+Ax_on;


