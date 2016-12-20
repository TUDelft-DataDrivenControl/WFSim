function [Qsp, Bsp]=Solution_space(BQ,BQy,Bt)

A = [BQ', BQy'];

[SLeft SRight] = spspaces(A,2);

Qsp1 = SRight{1}(:,SRight{3});
Qsp  = sparse(size(Qsp1,1),size(Qsp1,2));
%% sparsification
eps=1e-6;

for i=1:size(Qsp1,2)
    ind        = find(sign(abs(Qsp1(:,i))-eps)+1);
    Qsp(ind,i) = Qsp1(ind,i);
end

%% code from Steffen
%Parameter.SOWFI.Q_P = SRight{1}(:,SRight{3});
 
 %   Parameter.SOWFI.Q_P = spfun(@(M)SparcifyMatrix(M,MatrixEPS),Parameter.SOWFI.Q_P);
 
%%

clear Qsp1;

Bsp = A(1:end-1,:)\Bt(1:end-1,:);

%% check

% norm(full([BQ' BQy']*Bsp-Bt))
% norm(full([BQ' BQy']*Qsp))
