function Ais = precompute_rom_operator_param(f,Vn,i,d)


for i=1:d
    mu_i = one_hot(i,d);
    precompute_rom_operator(@(X) f(X,mu_i),Vn,i);
end