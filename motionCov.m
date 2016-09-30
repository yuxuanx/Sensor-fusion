function [ summCov ] = motionCov( chi, x, dim )
%Calculate the covariance summation in prediction step
summCov = zeros(dim,dim);
for i = 1:2*dim
    summCov = summCov + (chi(:,i) - x)*(chi(:,i) - x)';
end

