function ind = subs2ind(sizes,subs)

prods = [1 cumprod(sizes)];
subs = subs-1;
subs(1) = subs(1) + 1;
subs = [subs 0];

ind = dot(prods,subs);