function [F2,kps,k_sum] = F2_taylor_approx(k0_,mu_,qH,F2_k,xis,mu0)

kps = zeros(size(xis,1),qH);

for p=1:qH
    kp_ = diff(k0_,mu_,p-1)/factorial(p-1);
    kps(:,p) = eval(kp_(xis,mu0));
end

% k_sum = @(theta) sum(theta'.*kps,2);
k_sum = @(theta) kps*theta(:);

%     function val = F2_(x1,x2,theta)
%         val = 0;
%         for p=1:qH
%             val = val + theta(p)*F2_k(x1,x2,kps(:,p));
%         end
%     end
% 
% F2 = @(x1,x2,theta) F2_(x1,x2,theta); % not sure whether this nestedness is necessary

F2 = @(x1,x2,theta) F2_k(x1,x2,k_sum(theta));

end



