function X = rank_suff_basis(n,is)

% n: ROM dimension
% is: arrey of polynomial degrees

one_hot = @(i,n) sparse(i,1,1,n,1);

ni = numel(is);

for k =1:ni
    i = is(k);

    X = rec_basis(i,n);
end



    function Xi = rec_basis(i,n)

        if i == 1
            Xi = eye(n);
        else
            Xi = [];
            for m = 1:n
                one_hot_m = one_hot(m,m);
                Xim = [one_hot_m, rec_basis(i-1,m) + one_hot_m];
                %% add zeros to fit dimension n
                cols = size(Xim,2);
                Xim = [Xim; zeros(n-m,cols)];
                %%
                Xi = [Xi Xim];
            end
        end
    end

end