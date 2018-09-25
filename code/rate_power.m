
clear all
close all
clc
K = 4; %BS number
W_max = 1; %Total bandwidth 1GHz
W_min = 0.02; %Minimum bandwidth
% R_all = ones(K,1).*[0.7;0.8;0.9;1]; %Radii of all BSs 1km
R_all = ones(K,1)*1e3;
lambda_all = ones(K,1)*10/(pi*1e6);
%lambda_all = ones(K,1)*10/(pi*1e6).*[1;0.9;0.8;0.7]; %User arriving rate of all BSs
P_all = ones(K,1)*db2pow(30-30); %Maximum power of  all BSs 30dBm
noise_power = db2pow(-0-30);%Noise power -60dBm
alpha = 3;%Path loss 
N_all = pi*R_all.^2.*lambda_all;%Expectation of user number

%% Theoretic
L = 10; % Number of points used in Gaussian-Laguerre quadrature
[x_i, w_i] = pyGaussLaguerre(L); % points and weights in Gaussian-Laguerre quadrature
a = zeros(L,K);
b = zeros(L,K);
for i = 1:L
    for j = 1:K
        a(i,j) = w_i(i)*2*P_all(j)/(alpha*R_all(j)^2)*gammainc(x_i(i)*(R_all(j)^alpha),2/alpha)*gamma(2/alpha)*x_i(i)^(-2/alpha)/noise_power;
        b(i,j) = x_i(i)*P_all(j)/noise_power;
    end
end
%a=2*gammainc(x_i*(R_all.^alpha)',2/alpha).*gamma(2/alpha);
%b=alpha*x_i.^(2/alpha)*(R_all.^2)'*noise_power;
%b=alpha*x_i.^(2/alpha)*(R_all.^2)';
%c=alpha*x_i.^(2/alpha+1)*(R_all.^2)';
%d = P_all'.*a./b;
%e = (P_all.^2)'.*a.*c.*b.^(-2);
%f = P_all'./b.*c;
cvx_begin
cvx_quiet true
cvx_precision best
variable W_i(K,1)
for i = 1:L
    for j = 1:K
        tmp(i,j) = a(i,j)-a(i,j)*b(i,j)*inv_pos(W_i(j)+b(i,j));
    end
end
maximize sum(sum(tmp))
subject to
sum(W_i) <= W_max
W_i >= W_min*N_all
cvx_end
R_theoretic = sum(sum(tmp))/log(2)
%% 

%% Simulation
% 
MAX_IMPLEMENT = 1e5;
R_simulation = zeros(K,MAX_IMPLEMENT);
N = zeros(K,MAX_IMPLEMENT);
P_ij = zeros(K,MAX_IMPLEMENT);
for i = 1:K
    N(i,:) = poissrnd(lambda_all(i)*pi*R_all(i).^2,[1,MAX_IMPLEMENT]);
end
W_ij = W_max/K./N;
for i = 1:K
    P_ij(i,:) = P_all(i)./N(i,:);
end
g=cell(K,MAX_IMPLEMENT);
parfor i = 1:MAX_IMPLEMENT
    for j = 1:K
        if N(j,i) == 0
            R_simulation(j,i) = 0;
        else
        %clear h r g
        h = (randn(N(j,i),1)+1j*randn(N(j,i),1))/sqrt(2);
        r = rand(N(j,i),1).*R_all(j);
        g = h.*conj(h)./(1+r.^alpha);
        rate_all = W_ij(j,i).*log2(1+P_ij(j,i).*g./(W_ij(j,i)*noise_power));
        R_simulation(j,i) = sum(rate_all);
        end
    end
end
R_simulation_average = sum(sum(R_simulation))/MAX_IMPLEMENT


