function mpc = matsys(Wp,mpc)


Am  = mpc.A;
Bm  = mpc.B;
Bmt = mpc.Bt;
Cm  = mpc.C;


Nh = mpc.Nh;
nx = size(Am,1);
nu = Wp.turbine.N;

A = Am;
Ai = A;
AA = Ai;
for ii = 2:Nh
    Ai = A*Ai;
    AA = [AA;Ai];
end
mpc.AA  = AA;


AiB = Bm;
BB = kron(eye(Nh),AiB);
for ii = 1:Nh-1
    AiB = A*AiB;
    BB = BB+kron(diag(ones(Nh-ii,1),-ii),AiB);
end
mpc.BB  = BB;

AiBt = Bmt;
q    = num2cell(AiBt,[1,2]);
BBt  = blkdiag(q{:});
for ii = 1:Nh-1
    AA = A;
    for jj = ii:Nh-1
       BBt(jj*nx+1:(jj+1)*nx , (ii-1)*nu+1:ii*nu) = AA*AiBt(:,:,ii);
       AA = AA*A;
    end   
end
mpc.BBt = BBt;

C      = Cm;
Cc     = kron(eye(Nh),C);
mpc.CC  = Cc;


end

