function Ai_ = rec_precompute(ms,i,n,Vn,f)
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
                Ai_ = [Ai_ rec_precompute(ms_,i,n,Vn,f)];
            end
        end
            
    end