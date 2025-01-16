


N = 3;
i = 2;

A = reshape(1:N^i,N*ones(1,i));

P = perms(1:i);

B = A;
for k = 1:size(P,2)
    B = max(B,permute(A,P(k,:)));
end

b = B(:)';

u = unique(b);


M = repmat(u',1,N^i) == b;

s = (1:numel(u))*M;