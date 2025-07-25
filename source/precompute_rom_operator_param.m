function Ais = precompute_rom_operator_param(f,Vn,i,q)

n = size(Vn,2);
Ais = zeros(n,n_is(n,i),q);
for k=1:q
    theta_k = one_hot(k,q);
    Ais(:,:,k) = precompute_rom_operator(@(X) f(X,theta_k),Vn,i);
end