function [Qsp, Bsp]=Solution_space(B1,B2,bc)

A = [B1' B2'];

[~,SRight] = spspaces(A,2);

Qsp1 = SRight{1}(:,SRight{3});
Qsp  = sparse(size(Qsp1,1),size(Qsp1,2));

% Sparsification
eps=1e-6;

for i=1:size(Qsp1,2)
    ind        = find(sign(abs(Qsp1(:,i))-eps)+1);
    Qsp(ind,i) = Qsp1(ind,i);
end

Bsp = A(1:end-1,:)\bc(1:end-1,:);

% Check
nnz(full(Qsp'*[B1;B2]));     % momentum eq.
nnz(full([B1' B2']*Qsp));    % continuity eq.

clear Qsp1 SRight

%%
[~,SRight] = spspaces(B1',2);
PP         = SRight{1}(:,SRight{3});
P1         = sparse(size(PP,1),size(PP,2));   % B1' is full rank hence P1 is empty

% Sparsification
for i=1:size(PP,2)
    ind       = find(sign(abs(PP(:,i))-eps)+1);
    P1(ind,i) = PP(ind,i);
end
 
clear PP SRight

[~,SRight] = spspaces(B2',2);
PP         = SRight{1}(:,SRight{3});
P2         = sparse(size(PP,1),size(PP,2));   % B2' is not full rank

% Sparsification
for i=1:size(PP,2)
    ind       = find(sign(abs(PP(:,i))-eps)+1);
    P2(ind,i) = PP(ind,i);
end

clear PP SRight

P1 = sparse(size(P1,1),size(P2,2)); % This since B1' does not have a nullspace

P          = [P1;P2];

% Check
nnz(full(P'*[B1;B2]));     % momentum eq.
nnz(full([B1' B2']*P));    % continuity eq.


%Qsp = P;