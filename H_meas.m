function [ measurement ] = H_meas( chi, dim, s )
%Calculate the measurements
measurement = zeros(2,2*dim);
for i = 1:2*dim
    H_x = sqrt((chi(1,i)-s(1))^2+(chi(2,i)-s(2))^2);
    H_y = atan2(chi(2,i)-s(2),chi(1,i)-s(1));
%     if H_y < 0
%         H_y = H_y + pi;
%     end
    measurement(:,i) = [H_x;H_y];
end

end

