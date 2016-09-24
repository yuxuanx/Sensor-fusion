function H = H_measurement( x,s1,s2 )
H = zeros(2,1);
x_position = x(1);
y_position = x(2);
H(1) = atan2(y_position - s1(2),x_position - s1(1));
H(2) = atan2(y_position - s2(2),x_position - s2(1));


end

