function [ summCov ] = P_kk( A,chi,x_predict,x_update,dim )
summCov = zeros(dim,dim);
for i = 1:2*dim
    summCov = summCov + (chi(:,i) - x_update)*(A*chi(:,i) - x_predict)';
end


end

