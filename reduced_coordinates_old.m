function [s,b,u] = reduced_coordinates_old(N,i)
% best explained by example:
% [s,b,u] = reduced_coordinates(3)
%
% u: list of indices of unique entries
%   - in case of ambiguity, largest index is given
% b: full list of indices where ambigious entries' indices are replaced by the index listed in u
% s: b, but values are replaced by numbers from 1 to numel(u) while maintaining the ordering of these values 

if(nargin < 2)
    i = 2;
end

if i == 1
    s = 1:N;
    b = s;
    u = s;
else

%% only i=2
% A = reshape(1:N^2,N,N);
% 
% B = max(A,A');
%%

A = reshape(1:N^i,N*ones(1,i));

P = perms(1:i);

%% initial approach
% B = A;
% for k = 1:size(P,1) % can be optimized by taking max on a binary tree
%% little optimization
B = permute(A,P(1,:));
for k = 2:size(P,1) % can be optimized by taking max on a binary tree
    %%
    B = max(B,permute(A,P(k,:)));
end

b = B(:)';

[u,~,s] = unique(b);


% M = repmat(u',1,N^i) == b;
% 
% s = (1:numel(u))*M;

end