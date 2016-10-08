function [ x,w ] = resample( x,w )
%Resampling in Partical filter
n = length(w);
u = sort(rand(n,1));
px_cum = cumsum(w);
[~,idx1] = sort([u;px_cum]);
idx2 = find(idx1<=n);
idx = idx2 - (0:n-1)';
w = ones(n,1)/n;
if size(x,2) == 1
    x = x(idx);
else
    x = x(idx,:);
end

end

