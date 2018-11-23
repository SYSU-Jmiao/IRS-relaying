%% initialize parameters
clearvars
close all
clc
tic
MAX_IMPLEMENT = 1e1;
K_all = 4:14;
K = 6; %BS number
Wmax_all = 1e8*(0.1:0.1:1);
Wmax = 1e8; %Total bandwidth 100MHz
Rmin_all = 1e6*(0.1:0.1:1);
Rmin = 1e6; %10M bps
radius = 100;
Gamma = db2pow(-174+9.8+8-10+15.3-30);
% l1 = rand(K,1)*0.5+0.5;
l1 = 0.5+0.1*(1:K);
% l2 = rand(K,1)*0.5+0.5;
l2 = 1.1-0.1*(1:K);
R_all = radius.*l1;
lambda_all = 1e3./(1e6).*l2;
noise_power = db2pow(0);%thermal noise power density -174dBm/Hz
noise_power = noise_power*Gamma;
alpha = 3.76;%Path loss
options = optimoptions('fmincon','Algorithm','interior-point','display','off');

g=cell(length(Rmin_all),K,MAX_IMPLEMENT); %CSI
w=cell(length(Rmin_all),K,MAX_IMPLEMENT);
p=cell(length(Rmin_all),K,MAX_IMPLEMENT);
w1=cell(length(Rmin_all),K,MAX_IMPLEMENT); %without reservation
p1=cell(length(Rmin_all),K,MAX_IMPLEMENT);
w2=cell(length(Rmin_all),K,MAX_IMPLEMENT); %without DRA
p2=cell(length(Rmin_all),K,MAX_IMPLEMENT);
w3=cell(length(Rmin_all),K,MAX_IMPLEMENT); %without reservation and DRA
p3=cell(length(Rmin_all),K,MAX_IMPLEMENT);
N = zeros(length(Rmin_all),K,MAX_IMPLEMENT);
N_all = pi*R_all.^2.*lambda_all;%Expectation  of user number
%% Stage 1
for rate_i = 1:length(Rmin_all)
    Rmin = Rmin_all(rate_i);
    W0 = 1/K;
    N0 = W0*Wmax*noise_power;
    P0 = zeros(K,1);
    W0 = W0*ones(K,1);
    for i = 1:K
        f = @(P) rateExpectation(W0(i),P,R_all(i), alpha, noise_power, Wmax)-Rmin*N_all(i);
        P0(i) = fzero(f,1);
    end
    x0 = [W0;P0];
    % x0 = zeros(2*K,1);
    W = W0;
    P = P0;
    [W,P] = MNVO_Allocation([W;P],R_all,N_all,alpha, noise_power, Rmin, Wmax);
    %     tic
    %     [W2,P2] = ADMM_MNVO_Allocation(x0,R_all,N_all,alpha, noise_power, Rmin, Wmax);
    %     toc;
    %%
    %% Stage 2
    for i = 1:K
        N(rate_i,i,:) = poissrnd(lambda_all(i)*pi*R_all(i).^2,[1,MAX_IMPLEMENT]);
    end
    
    %%
    parfor imp_i = 1:MAX_IMPLEMENT
        for j = 1:K
            r_j = R_all(j);
            wmax = W(j);
            wmax0 = W0(j);
            
            n = N(rate_i,j,imp_i);
            if n > 0
                h = (randn(n,1)+1j*randn(n,1))/sqrt(2);
                r = sqrt(rand(n,1).*r_j.^2);
                g{rate_i,j,imp_i} = h.*conj(h)./(1+r.^alpha);
                w0 = wmax/n*ones(n,1);
                p0 = noise_power*Wmax*(2.^(Rmin./(Wmax*w0))-1).*w0./g{rate_i,j,imp_i};
                w2{rate_i,j,imp_i} = w0;
                p2{rate_i,j,imp_i} = p0;
                w1_0 = wmax0/n*ones(n,1);
                p1_0 = noise_power*Wmax*(2.^(Rmin./(Wmax*w1_0))-1).*w1_0./g{rate_i,j,imp_i};
                w3{rate_i,j,imp_i} = w1_0;
                p3{rate_i,j,imp_i} = p1_0;
                if n==1
                    w{rate_i,j,imp_i} = w0;
                    p{rate_i,j,imp_i} = p0;
                else
                    %clear h r g
                    %             tic;
                    %             [w{j,i},p{j,i}]=userAllocation2(W(j),g{j,i},[w0;p0]);
                    %             [w1{j,i},p1{j,i}]=userAllocation2(W0(j),g{j,i},[w1_0;p1_0]);
                    %             [w{j,i},p{j,i}]=userAllocation(W(j),g{j,i});
                    %             [w1{j,i},p1{j,i}]=userAllocation(W0(j),g{j,i});
                    %             tic;
                    [w{rate_i,j,imp_i},p{rate_i,j,imp_i}]=userAllocation3(wmax,g{rate_i,j,imp_i},Wmax, noise_power, Rmin);
                    %             toc;
                    [w1{rate_i,j,imp_i},p1{rate_i,j,imp_i}]=userAllocation3(wmax0,g{rate_i,j,imp_i},Wmax, noise_power, Rmin);
                    %             toc;
                end
            end
        end
    end
end
p_result = (cellfun(@sum,p));
p_result1 = (cellfun(@sum,p1));
p_result2 = (cellfun(@sum,p2));
p_result3 = (cellfun(@sum,p3));

p_sum_result = squeeze(sum(p_result,2));
p_sum_result1 = squeeze(sum(p_result1,2));
p_sum_result2 = squeeze(sum(p_result2,2));
p_sum_result3 = squeeze(sum(p_result3,2));

f=@geo_mean;
p_sum_result_mean = squeeze(f(p_sum_result,2));
p_sum_result_mean1 = squeeze(f(p_sum_result1,2));
p_sum_result_mean2 = squeeze(f(p_sum_result2,2));
p_sum_result_mean3 = squeeze(f(p_sum_result3,2));
%%
figure;
plot(pow2db(p_sum_result_mean),'-')
hold on
plot(pow2db(p_sum_result_mean1),'-*')
plot(pow2db(p_sum_result_mean2),'-+')
plot(pow2db(p_sum_result_mean3),'-x')
legend('Proposed','DRA','Reservation', 'Benchmark');