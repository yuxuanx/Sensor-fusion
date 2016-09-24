x0 = [0 0 1 0 0 0]';
T = 0.01;
A_vel = [eye(2) T*eye(2);zeros(2,2) eye(2)];
A_acc = [eye(2) T*eye(2) T^2/2*eye(2);...
    zeros(2,2) eye(2) T*eye(2);...
    zeros(2,2) zeros(2,2) eye(2)];
Q = [zeros(2,4);zeros(2,2) eye(2)]*T*10^2;
Q = [zeros(2,6);zeros(2,6);zeros(2,4) eye(2)]*T*10^2;
K = 20;
x = zeros(6,K);
x(:,1) = x0;
P = zeros(6,6,K);
for i = 1:K-1
    x(:,i+1) = A_acc*x(:,i);
    P(:,:,i+1) = A_acc*P(:,:,i)*A_acc' + Q;
end

T = 0.01;
A_vel = [eye(2) T*eye(2);zeros(2,2) eye(2)];
for i = 1:K
    X_hat(:,i) = A_vel*X(:,i);
end
diff = X_hat(:,1:end-1) - X(:,2:end);
err = sum(diff.^2,2)/(K-1);
