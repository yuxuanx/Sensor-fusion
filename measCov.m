function [ meas_cov ] = measCov( chi, dim, s, y )
measurement = H_meas( chi, dim, s );
meas_cov = zeros(2,2,2*dim);
for i = 1:2*dim
    meas_cov(:,:,i) = (measurement(:,i) - y)*(measurement(:,i) - y)';
end

end

