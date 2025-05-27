function subs = ind2subs(sizes,ind)

dim = numel(sizes);
prods = cumprod(sizes);

subs = zeros(1,dim);
for i = dim:-1:2
    prod = prods(i-1);
    sub = floor((ind-1)/prod);
    subs(i) = sub;
    ind = ind - prod*sub;
end
subs = subs+1;

subs(1) = ind;




