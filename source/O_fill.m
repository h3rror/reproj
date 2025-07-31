function O_f = O_fill(O,is,qs,s)
% add zeros to transform aggregated operator with qs affine terms into 
% operator with s affine terms per operator

O_f = O_change_qs(O,is,qs,s*ones(size(qs));