function [s,b,u] = reduced_coordinates(N,i)
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
    B = reshape(1:N^i,N*ones(1,i));
    b = B(:)';

    inds = [];
    b = rec_for_loop(b,inds,1,i,N);


end

    function b = rec_for_loop(b,inds,j,i,N)

        if j==i+1
            inds_sorted = sort(inds);
            if inds_sorted == inds
                % inds
                % B(inds) = counter;
                % counter = counter + 1;
            else
                % inds
                b(subs2ind(N*ones(1,i),inds)) = b(subs2ind(N*ones(1,i),inds_sorted));
            end
        else
            for m=N:-1:1
                inds(j) = m;
                b = rec_for_loop(b,inds,j+1,i,N);
            end
        end
    end

    [u,~,s] = unique(b);


end