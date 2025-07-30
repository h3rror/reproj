function O = affine_op(Os,theta)
% 
% [dim1,dim2,q] = size(Os);
% 
% O = zeros(dim1,dim2);
% for i=1:q
%     O = O + theta(i)*Os(:,:,i);
% end

O = sum(pagemtimes(permute(theta(:),[3 2 1]),Os),3);