function Ai = precompute_rom_operator(f,Vn,i)
% f: i-th order FOM operator as function of i inputs x
% Vn: ROM basis
% i: order of operator

n = size(Vn,2);

ms = [];
Ai = rec_precompute(ms,i,n);
% Ai = rec_precompute(ms,i,n,Vn,f);

    function Ai_ = rec_precompute(ms,i,n)
        % ms: array of mode indices
        % i: order of operator
        % j: current recursion depth

        Ai_ = [];
        %% enforce incrementality structure
        % if isempty(ms) 
            m_max = n;
        % else
        %     m_max = ms(end);
        % end
    %%
        for m = 1:m_max 
            ms_ = [ms m];
            if length(ms_) == i
                Ai_ = [Ai_ Vn'*f(Vn(:,ms_))];
                % ms_
            else
                Ai_ = [Ai_ rec_precompute(ms_,i,n)];
            end
        end

    end

end

% note: the ms correspond to the rank suff basis
% -> can be used to debug/verify