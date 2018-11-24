%% initialize parameters
clearvars
close all
clc
MAX_IMPLEMENT = 1e4;
hasSimulation = 1;
K_all = 4:14;
K = 6; %BS number
Wmax_all = 1e8*(1:1:10);
Wmax = 1e8; %Total bandwidth 100MHz
Rmin_all = 1e6*(0.1:0.1:1);
Rmin = 1e6; %1M bps
radius = 100;
Gamma = db2pow(-174+9.8+8-10+15.3-30);
% l1 = rand(K,1)*0.5+0.5;
% l1 = ones(1,K);
% l1 = 1-K/2*0.1+0.1*(1:K);
l1 = [0.8,0.8,1.0,1.0,1.2,1.2];
l2 = [1.2,0.8,0.8,1.0,1.0,1.2];
R_all = radius.*l1;
% l2 = 1+K/2*0.1-0.1*(1:K);
% l2 = l1;

% l2 = ones(1,K);
lambda_all = 2e3./(1e6).*l2;
% lambda_all = 2e3./(l1.^2*1e6);
N_all = pi*R_all.^2.*lambda_all;%Expectation  of user numbernoise_power = db2pow(0);%thermal noise power density -174dBm/Hz
noise_power = db2pow(0);%thermal noise power density -174dBm/Hz
noise_power = noise_power*Gamma;
alpha = 3.76;%Path loss
alpha_all = 3:0.2:5;
options = optimoptions('fmincon','Algorithm','interior-point','display','off');

% g=cell(length(Rmin_all),K,MAX_IMPLEMENT); %CSI
% w=cell(length(Rmin_all),K,MAX_IMPLEMENT);
% p=cell(length(Rmin_all),K,MAX_IMPLEMENT);
% w1=cell(length(Rmin_all),K,MAX_IMPLEMENT); %without reservation
% p1=cell(length(Rmin_all),K,MAX_IMPLEMENT);
% w2=cell(length(Rmin_all),K,MAX_IMPLEMENT); %without DRA
% p2=cell(length(Rmin_all),K,MAX_IMPLEMENT);
% w3=cell(length(Rmin_all),K,MAX_IMPLEMENT); %without reservation and DRA
% p3=cell(length(Rmin_all),K,MAX_IMPLEMENT);
N = zeros(length(Rmin_all),K,MAX_IMPLEMENT);
p_result = zeros(length(Rmin_all),K,MAX_IMPLEMENT);
p_result1 = zeros(length(Rmin_all),K,MAX_IMPLEMENT);
p_result2 = zeros(length(Rmin_all),K,MAX_IMPLEMENT);
p_result3 = zeros(length(Rmin_all),K,MAX_IMPLEMENT);

W_all = zeros(length(Rmin_all),K);
W2_all = zeros(length(Rmin_all),K);
W = 1/K*ones(K,1);
%% Stage 1
for rate_i = 1:length(Rmin_all)
    tic
    fprintf("r = %d\n", rate_i);
%     Rmin = Rmin_all(rate_i);
%         Wmax = Wmax_all(rate_i);
        alpha = alpha_all(rate_i);
    %     W0 = W;
    W0 = 1/K*ones(K,1);
    N0 = W0*Wmax*noise_power;
    P0 = zeros(K,1);
    
    for i = 1:K
        f = @(x) rateExpectation(W0(i),x,R_all(i), alpha, noise_power, Wmax)/(Rmin*N_all(i))-1;
        P0(i) = fzero(f,1);
    end
    %     if sum(isnan(W0)+isnan(P0)) == 0
    %         x0 = [W0;P0];
    %     else
    %         x0 = zeros(2*K,1);
    %     end
    W = W0;
    P = P0;
    [W,P] = MNVO_Allocation([W;P],R_all,N_all,alpha, noise_power, Rmin, Wmax);
    W_all(rate_i,:) = W';
    tic
    [W2,P2,hist(rate_i)] = ADMM_MNVO_Allocation([W0;P0],R_all,N_all,alpha, noise_power, Rmin, Wmax);
    W2_all(rate_i,:) = W2';
    toc;
    %%
    %% Stage 2
    for i = 1:K
        N(rate_i,i,:) = poissrnd(lambda_all(i)*pi*R_all(i).^2,[1,MAX_IMPLEMENT]);
    end
    
    %%
    if hasSimulation
        parfor imp_i = 1:MAX_IMPLEMENT
            for j = 1:K
                r_j = R_all(j);
                wmax = W(j);
                wmax0 = W0(j);
                
                n = N(rate_i,j,imp_i);
                if n > 0
                    h = (randn(n,1)+1j*randn(n,1))/sqrt(2);
                    r = sqrt(rand(n,1).*r_j.^2);
                    g= h.*conj(h)./(1+r.^alpha);
                    w0 = wmax/n*ones(n,1);
                    p0 = noise_power*Wmax*(2.^(Rmin./(Wmax*w0))-1).*w0./g;
                    
                    w2 = w0;
                    p2 = p0;
                    p_result2(rate_i,i,imp_i) = sum(p2);
                    
                    w1_0 = wmax0/n*ones(n,1);
                    p1_0 = noise_power*Wmax*(2.^(Rmin./(Wmax*w1_0))-1).*w1_0./g;
                    
                    w3 = w1_0;
                    p3 = p1_0;
                    p_result3(rate_i,i,imp_i) = sum(p3);
                    if n==1
                        w = w0;
                        p = p0;
                        p_result(rate_i,i,imp_i) = sum(p0);
                    else
                        %clear h r g
                        %             tic;
                        %             [w{j,i},p{j,i}]=userAllocation2(W(j),g{j,i},[w0;p0]);
                        %             [w1{j,i},p1{j,i}]=userAllocation2(W0(j),g{j,i},[w1_0;p1_0]);
                        %             [w{j,i},p{j,i}]=userAllocation(W(j),g{j,i});
                        %             [w1{j,i},p1{j,i}]=userAllocation(W0(j),g{j,i});
                        %             tic;
                        [w,p]=userAllocation3(wmax,g,Wmax, noise_power, Rmin);
                        p_result(rate_i,i,imp_i) = sum(p);
                        %             toc;
                        [w1,p1]=userAllocation3(wmax0,g,Wmax, noise_power, Rmin);
                        p_result1(rate_i,i,imp_i) = sum(p1);
                        %             toc;
                    end
                end
            end
        end
    end
    toc
end
nowTime = datestr(now,'yyyymmddHHMMSS');
nowTime = ['rate_rmin', nowTime];
save(nowTime);
if hasSimulation
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
    xt = alpha_all(1:10);
    xlabel('\alpha')
%     xt=Rmin_all/1e6;
%     xlabel('R(Mbps)')
%         xt=Wmax_all/1e6;
%         xlabel('W_{max}(MHz)')
    hold on
    plot(xt,30+pow2db(p_sum_result_mean),'-k*')
    plot(xt,30+pow2db(p_sum_result_mean1),'-ks')
    plot(xt,30+pow2db(p_sum_result_mean2),'-kv')
    plot(xt,30+pow2db(p_sum_result_mean3),'-ko')
    ylabel('P_{sum}(dBm)')
    legend('Proposed','DRA','Reservation', 'Benchmark');
end
%%
markers = ['+','x','*','^','s','o'];
l_c = strings(6,1);
for i = 1:K
    l_c(i) = '('+string(l1(i))+','+string(l2(i))+')';
end

figure;
hold on
ylabel('Spectrum ratio');
% xlabel('R_{min}(Mbps)');
%         xlabel('W_{max}(MHz)')
xlabel('\alpha')
for i = 1:K
    plot(xt,W_all(:,i),['-k',markers(i)]);
end
legend(l_c)

figure
hold on
ylabel('Spectrum ratio');
xlabel('Iteratons');
iter_i = 5;
for i = 1:K
    plot(hist(iter_i).W(i,:),['-k',markers(i)]);
end
legend(l_c)