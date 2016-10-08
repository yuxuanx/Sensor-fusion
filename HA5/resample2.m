function [ x,w ] = resample2( x,w,N )
u = sort(rand(N,1));
px_cum = cumsum(w);
[~,idx1] = sort([u;px_cum]);
idx2 = find(idx1<=N);
idx = idx2 - (0:N-1)';
w = ones(N,1)/N;

x = x(idx,:);


end

