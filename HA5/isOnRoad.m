function [u] = isOnRoad(x,y);

n = prod(size(x));
m = 9;
%column vectors
x = reshape(x,n,1); 
y = reshape(y,n,1);
X = x*ones(1,m);
Y = y*ones(1,m);

bounds = ([1+i 1+9*i 11+9*i 11+i]);
house = zeros(m,5);
house(1,:) = ([2+5.2*i 2+8.3*i 4+8.3*i 4+5.2*i 2+5.2*i]);%House 1
house(2,:) = ([2+3.7*i 2+4.4*i 4+4.4*i 4+3.7*i 2+3.7*i]);%House 2
house(3,:) = ([2+2*i 2+3.2*i 4+3.2*i 4+2*i 2+2*i]);%House 3
house(4,:) = ([5+i 5+2.2*i 7+2.2*i 7+i 5+i]);%House 4
house(5,:) = ([5+2.8*i 5+5.5*i 7+5.5*i 7+2.8*i 5+2.8*i]);%House 5
house(6,:) = ([5+6.2*i 5+9*i 7+9*i 7+6.2*i 5+6.2*i]);%House 6
house(7,:) = ([8+4.6*i 8+8.4*i 10+8.4*i 10+4.6*i 8+4.6*i]);%House 7
house(8,:) = ([8+2.4*i 8+4*i 10+4*i 10+2.4*i 8+2.4*i]);%House 8
house(9,:) = ([8+1.7*i 8+1.8*i 10+1.8*i 10+1.7*i 8+1.7*i]);%House 9

X1 = X >= ones(n,1)*real(house(:,1))';
X2 = X <= ones(n,1)*real(house(:,3))';
Y1 = Y >= ones(n,1)*imag(house(:,1))';
Y2 = Y <= ones(n,1)*imag(house(:,2))';

x3 = x > ones(n,1)*real(bounds(1))';
x4 = x < ones(n,1)*real(bounds(3))';
y3 = y > ones(n,1)*imag(bounds(1))';
y4 = y < ones(n,1)*imag(bounds(2))';

XX = X1.*X2;
YY = Y1.*Y2;
UU = XX.*YY;
u1 = 1-min(1,(sum(UU')))'; 

xx = x3.*x4;
yy = y3.*y4;
u2 = xx.*yy;

u = u1.*u2;


