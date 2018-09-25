clear all
close all
clc
MAX_IMPLEMENT = 100;
K = 4; %BS number
W_max = 1; %Total bandwidth 1GHz
W_min = 0.02; %Minimum bandwidth
% R_all = ones(K,1).*[0.7;0.8;0.9;1]; %Radii of all BSs 1km
R_all = ones(K,1);
lambda_all = ones(K,1)*10/(pi).*[1;0.9;0.8;0.7]; %User arriving rate of all BSs
P_all = ones(K,1)*db2pow(30-30); %Maximum power of  all BSs 30dBm
noise_power = db2pow(30-30);%Noise power -60dBm
alpha = 3;%Path loss 
N_all = pi*R_all.^2.*lambda_all;%Expectation of user number

%% Theoretic
L = 50; % Number of points used in Gaussian-Laguerre quadrature
[x_i, w_i] = pyGaussLaguerre(L); % points and weights in Gaussian-Laguerre                                                                                                                                                                                                                                                                                                                                                                                   quadrature
a=2*gammainc(x_i*(R_all.^alpha)',2/alpha).*gamma(2/alpha);
b=alpha*x_i.^(2/alpha)*(R_all.^2)'*noise_power;
%b=alpha*x_i.^(2/alpha)*(R_all.^2)';
c=alpha*x_i.^(2/alpha+1)*(R_all.^2)';
cvx_begin
cvx_quiet true
cvx_precision best
variable W_i(K,1)

for i = 1:K
    tmp = P_all(i).*a(:,i)./b(:,i)-P_all(i)^2.*a(:,i).*c(:,i).*inv_pos(W_i(i).*b(:,i).^2+P_all(i).*b(:,i).*c(:,i));
    rate(i) = sum(tmp);
end
maximize sum(rate)
subject to
sum(W_i) <= W_max
W_i >= W_min*N_all
cvx_end
R_theoretic = sum(rate)/log(2)
%% 

%% Simulation
% 
MAX_IMPLEMENT = 100;
R_simulation = zeros(K,MAX_IMPLEMENT);
N = zeros(K,MAX_IMPLEMENT);
for i = 1:K
    N(i,:) = poissrnd(lambda_all(i)*pi*R_all(i).^2,[1,MAX_IMPLEMENT]);
end
for j = 1:K
    P_ij = P_all./N(:,i);
    W_ij = W_max/K./N;
parfor i = 1:MAX_IMPLEMENT
    for j = 1:K
        if N(j,i) == 0
            R_simulation(j,i) = 0;
        else
        %clear h r g
        h = 1/sqrt(2)*(randn(N(j,i),1)+i*randn(N(j,i),1));
        r = rand(N(j,i),1).*R_all(j);
        g = arrayfun(@(x) norm(x,2),h).^2./(1+r.^alpha);
        R_simulation(j,i) = sum(W_i(j)*log2(1+P_i(j).*g./(W_i(j)*noise_power)));
        end
    end
end
sum(sum(R_simulation))/MAX_IMPLEMENT