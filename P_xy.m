function [ Pxy ] = P_xy( chi, x, y, s, dim )
%Calculate the Pxy in updating step
measurements = H_meas(chi,dim,s);
Pxy = zeros(dim,2,dim*2);
for i = 1:2*dim
    Pxy(:,:,i) = (chi(:,i) - x)*(measurements(:,i) - y)';
end
end

