function [] = PostPlot(x_part,w_part,Xfilt,P_filt,N,K,SIR)
% Illustrating the posterior densities
    %Xp is x_part, W is w_part
    %xkk and Pkk: posterior moments (mean and covariance, K times)
    %SIR if 1: we are using SIR, if 0: we are using SIS
%   Notation, PF:   Xp and W are n x K matrices that contain N 
%                   particles in the interval k = 1,2, ..., K. 
%   Notation, KF:   xkk and Pkk are 1 x K vectors that contain 
%                   the posterior moments for k = 1,2, ..., K. 

Xp=x_part;      W=w_part;
xkk=Xfilt;      Pkk=P_filt;
%   We can plot the densities at all times
for k= 1:K

% Let us first determine the x-interval of interest:
xmin    =   min(Xp(:,k)); xmax  =   max(Xp(:,k));
X       =   linspace(xmin-(xmax-xmin)/3,xmax+(xmax-xmin)/3,800);

%       We can now construct a continuous approximation to the posterior
%       density by placing a Gaussian kernel around each particle
pApprox =   zeros(size(X));   % A vector that will contain the pdf values
if SIR==0
    sigma=0.5;               % The std of the kernel. 
                              % NOTE: the std has to be tuned to obtain a
                              % "nice" illustration.
                              % If you resample, you can set the std to
                              % (xmax-xmin)/sqrt(N);
elseif SIR==1                              
    sigma=(xmax-xmin)/sqrt(N);
end
for     i   =   1   :   N
    pApprox     =   pApprox     +   W(i,k)*normpdf(Xp(i,k),X,sigma);
end

%       We are now ready to plot the densities 
figure;
clf
plot(X,pApprox,'LineWidth',2)   % This is the PF approximation
hold on
plot(X,normpdf(xkk(k),X,sqrt(Pkk(k))),'r-.','LineWidth',2) % KF posterior density
legend('Particle filter approximation','Kalman filter')
title(['p(x_k |  y_{1:k}), k=',num2str(k)],'FontSize',20)
end


%%   Illustrating the weight distribution as a function of time
[Kmesh,Nmesh] = meshgrid((1:K)', (1:N)');
 
figure
clf
mesh(Kmesh,Nmesh,W,'LineWidth',2);
title('Weights, w_k^{(i)}','FontSize',20)
view(-36,22)
xlabel('time, k','FontSize',20)
ylabel('particle number, i','FontSize',20)
axis([1 K 1 N 0 max(max(W))])
end