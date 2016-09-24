function h = H_derivative(x,s,flag)
x_position = x(1);
y_position = x(2);
square = (y_position - s(2))^2 + (x_position - s(1))^2;
if flag == 1
    h = (s(2) - y_position)/square;
else
    h = (x_position - s(1))/square;
end
end