%% initialize parameters
clearvars
close all
clc
MAX_IMPLEMENT = 100;
hasSimulation = 1;
K_all = 4:14;
K = 1; %BS number
Wmax_all = 1e8*(1:1:10);
Wmax = 1e8; %Total bandwidth 100MHz
Rmin_all = 1e6*(0.2:0.2:2);
l0=1;
% l0 = [2,0.5,0.5,1.0,1.0,2.0];
Rmin = 1e6*l0;
 %1M bps
radius = 100;
Gamma = db2pow(-174+9.8+8-10+15.3-30);
% l1 = rand(K,1)*0.5+0.5;
% l1 = ones(1,K);
% l1 = 1-K/2*0.1+0.1*(1:K);
% l1 = [0.8,0.8,1.0,1.0,1.2,1.2];
l1=1;
l2=1;
% l2 = [1.2,0.8,0.8,1.0,1.0,1.2];
R_all = radius.*l1;
% l2 = 1+K/2*0.1-0.1*(1:K);
% l2 = l1;

% l2 = ones(1,K);
lambda_all = 2e3./(1e6).*l2;
% lambda_all = 2e3./(l1.^2*1e6);
N_all = pi*R_all.^2.*lambda_all;%Expectation  of user numbernoise_power = db2pow(0);%thermal noise power density -174dBm/Hz
f_N = @(x) (ei(x)-log(x)-eulergamma*ones(size(x))).*exp(-x);
N_all_1 = double(1./f_N(vpa(N_all)));

noise_power = db2pow(0);%thermal noise power density -174dBm/Hz
noise_power = noise_power*Gamma;
alpha = 3.76;%Path loss
alpha_all = 3:0.2:5;
options = optimoptions('fmincon','Algorithm','interior-point','display','off');

W_ben_all = zeros(length(Rmin_all),K);
W_eaq_all = zeros(length(Rmin_all),K);
W_ave_all = zeros(length(Rmin_all),K);
W_min_all = zeros(length(Rmin_all),K);
W_ave_test_all = zeros(length(Rmin_all),K);
W_min_test_all = zeros(length(Rmin_all),K);

P_ben_ave_all = zeros(length(Rmin_all),K);
P_ben_min_all = zeros(length(Rmin_all),K);
P_eaq_ave_all = zeros(length(Rmin_all),K);
P_eaq_min_all = zeros(length(Rmin_all),K);
P_ave_all = zeros(length(Rmin_all),K);
P_ave_test_all = zeros(length(Rmin_all),K);
P_min_all = zeros(length(Rmin_all),K);
P_min_test_all = zeros(length(Rmin_all),K);


P_ben_ave = zeros(K,1);
P_ben_min = zeros(K,1);
P_eaq_ave = zeros(K,1);
P_eaq_min = zeros(K,1);

% N = zeros(length(Rmin_all),K,MAX_IMPLEMENT);
p_result_ben = zeros(length(Rmin_all),K,MAX_IMPLEMENT);
p_result_eaq = zeros(length(Rmin_all),K,MAX_IMPLEMENT);
p_result_ave = zeros(length(Rmin_all),K,MAX_IMPLEMENT);
p_result_min = zeros(length(Rmin_all),K,MAX_IMPLEMENT);

load Ng_6_100.mat N g
%% Stage 1
for rate_i = 1:length(Rmin_all)
    tic
    fprintf("r = %d\n", rate_i);
    Rmin = Rmin_all(rate_i)*l0;
    %     Wmax = Wmax_all(rate_i);
    %     alpha = alpha_all(rate_i);
    %     W0 = W;
    W_ben = 1/K*ones(K,1); %% benchmark with equal allocation
    W_ben_all(rate_i,:) = W_ben';
    
    W_eaq = N_all'.*Rmin/sum(N_all'.*Rmin); %% allocation based on user number and QoS
    W_eaq_all(rate_i,:) = W_eaq';
    
    for i = 1:K
        %% average rate constraint
        f_ben_ave = @(x) rateExpectation(W_ben(i),x,R_all(i), N_all_1(i), alpha, noise_power, Wmax)/Rmin(i) - 1;
        P_ben_ave(i) = fzero(f_ben_ave,1);
        
        f_eaq_ave = @(x) rateExpectation(W_eaq(i),x,R_all(i), N_all_1(i), alpha, noise_power, Wmax)/Rmin(i) - 1;
        P_eaq_ave(i) = fzero(f_eaq_ave,1);
        
        %% minimum rate constraint
        f_ben_min = @(x) rateExpectation2(W_ben(i),x,R_all(i), N_all(i), alpha, noise_power, Wmax)/Rmin(i)-1;
        P_ben_min(i) = fzero(f_ben_min,1e-3);
        
        f_eaq_min = @(x) rateExpectation2(W_eaq(i),x,R_all(i), N_all(i), alpha, noise_power, Wmax)/Rmin(i)-1;
        P_eaq_min(i) = fzero(f_eaq_min,1e-3);
    end
    P_ben_ave_all(rate_i,:) = P_ben_ave';
    P_ben_min_all(rate_i,:) = P_ben_min';
    P_eaq_ave_all(rate_i,:) = P_eaq_ave';
    P_eaq_min_all(rate_i,:) = P_eaq_min';
    %     W = W0;
    %     P = P0;
    [W_ave_test,P_ave_test] = MNVO_Allocation([W_ben;P_ben_ave],R_all,N_all_1,alpha, noise_power, Rmin, Wmax);
    [W_ave,P_ave,iter_hist(rate_i)] = ADMM_MNVO_Allocation([W_ben;P_ben_ave],R_all,N_all_1,alpha, noise_power, Rmin, Wmax);
    norm(W_ave-W_ave_test)
%     P_ave_test = P_ave_test.*N_all;
%     P_ave = P_ave.*N_all;
    [W_min,P_min] = ADMM_MNVO_Allocation2([W_ben;P_ben_min],R_all,N_all,alpha, noise_power, Rmin, Wmax);
    [W_min_test,P_min_test] = MNVO_Allocation2([W_ben;P_ben_min],R_all,N_all,alpha, noise_power, Rmin, Wmax);
    norm(W_min-W_min_test)
    
    W_ave_all(rate_i,:) = W_ave';
    W_ave_test_all(rate_i,:) = W_ave_test';
    P_ave_all(rate_i,:) = P_ave';
    P_ave_test_all(rate_i,:) = P_ave_test';
    W_min_all(rate_i,:) = W_min';
    W_min_test_all(rate_i,:) = W_min_test';
    P_min_all(rate_i,:) = P_min';
    P_min_test_all(rate_i,:)=P_min_test';
    %%
    %% Stage 2
%     for i = 1:K
%         N(rate_i,i,:) = poissrnd(N_all(i),[1,MAX_IMPLEMENT]);
%     end
    
    %%
    if hasSimulation
        %%
        parfor imp_i = 1:MAX_IMPLEMENT
            for j = 1:K
                r_j = R_all(j);
                
                n = N(j,imp_i);
                if n > 0
%                     h = (randn(n,1)+1j*randn(n,1))/sqrt(2);
%                     r = sqrt(rand(n,1).*r_j.^2);
%                     g= h.*conj(h)./(1+r.^alpha);
                    if n==1
                        
                        p_result_ben(rate_i,j,imp_i) = noise_power*Wmax*(2.^(Rmin(j)./(Wmax*W_ben(j)))-1).*W_ben(j)./g{j,imp_i};
                        p_result_eaq(rate_i,j,imp_i) = noise_power*Wmax*(2.^(Rmin(j)./(Wmax*W_eaq(j)))-1).*W_eaq(j)./g{j,imp_i};
                        p_result_ave(rate_i,j,imp_i) = noise_power*Wmax*(2.^(Rmin(j)./(Wmax*W_ave(j)))-1).*W_ave(j)./g{j,imp_i};
                        p_result_min(rate_i,j,imp_i) = noise_power*Wmax*(2.^(Rmin(j)./(Wmax*W_min(j)))-1).*W_min(j)./g{j,imp_i};
                    else
                       
                        [w_ben,p_ben]=userAllocation3(W_ben(j),g{j,imp_i},Wmax, noise_power, Rmin(j));
                        p_result_ben(rate_i,j,imp_i) = sum(p_ben);
                        %             toc;
                        [w_eaq,p_eaq]=userAllocation3(W_eaq(j),g{j,imp_i},Wmax, noise_power, Rmin(j));
                        p_result_eaq(rate_i,j,imp_i) = sum(p_eaq);
                        %             toc;
                        [w_ave,p_ave]=userAllocation3(W_ave(j),g{j,imp_i},Wmax, noise_power, Rmin(j));
                        p_result_ave(rate_i,j,imp_i) = sum(p_ave);
                        
                        [w_min,p_min]=userAllocation3(W_min(j),g{j,imp_i},Wmax, noise_power, Rmin(j));
                        p_result_min(rate_i,j,imp_i) = sum(p_min);
                    end
                end
            end
        end
        %%
    end
    toc
end
nowTime = datestr(now,'yyyymmddHHMMSS');
nowTime = ['prediction_compare', nowTime];
save(nowTime);
if hasSimulation
    f=@mean;
    g=@sum;
    p1 = squeeze(f((squeeze(g(p_result_ben,2))),2));
    p2 = squeeze(f((squeeze(g(p_result_eaq,2))),2));
    p3 = squeeze(f((squeeze(g(p_result_ave,2))),2));
    p4 = squeeze(f(squeeze(g(p_result_min,2)),2));
    
    %%
    figure;
    %     xt = alpha_all(1:10);
    %     xlabel('\alpha')
    xt=Rmin_all/1e6;
    xlabel('R_{min}(Mbps)')
    xi = xt;
    y1 = 30+pow2db(P_ave_test_all);
    y2 = 30+pow2db(p3);
    save Prediction_compare.mat xi y1 y2
    %         xt=Wmax_all/1e6;
    %         xlabel('W_{sum}(MHz)')
    hold on
    plot(xt,y1,'-k*')
%     plot(xt,30+pow2db(P_min_all),'-ks')
%     plot(xt,30+pow2db(p1),'-kv')
%     plot(xt,30+pow2db(p2),'-ks')
    plot(xt,y2,'-kv')
%     plot(xt,30+pow2db(p4),'-ko')
    ylabel('P_{sum}(dBm)')
    legend('Prediction','Simulation');
end