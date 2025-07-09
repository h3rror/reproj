function P = uniquepowers(X,is)
% computes uniquepower(X,i) for all i in is and stacks the results
% vertically

P = [];
for i=is
    P = [P;uniquepower(X,i)];
end
