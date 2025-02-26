function xs = nnz_xs(nnz_,n,l)
% compute all n-dimensional vectors with nnz nonzero entries
% the nonzero entries are natural numbers and sum up to l

x = [];
xs = rec_nnz_xs(nnz_,n,l,x);

    function xs = rec_nnz_xs(nnz_,n,l,x)
% x is vector of variable length of already determined entries
% 
        if length(x)<=n
            if nnz_ == nnz(x)
                xs = sparse(1:length(x),1,x,n,1);
                % x = [x; zeros(n-length(x),1)];
            else
                xs = [];
                for xi = l-sum(x):-1:0
                    x_ = [x xi];
                    xs = [xs rec_nnz_xs(nnz_,n,l,x_)];
                end
            end
        else
            xs = [];
        end

    end

end