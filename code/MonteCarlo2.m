%% MNVO resource allocation with initial point x0
function [W_final,P_all_final] = MonteCarlo2(R_all,N_all,alpha, noise_power, Rmin, Wmax,T)
K = length(R_all);
N = zeros(K,T);
MAX_ITER = 2e3;
ABSTOL   = 1e-5;
RELTOL   = 1e-3;
rho = 1;
tau = 2;
t = 100;
W_all = zeros(K,T);
W_mean = 1/K*ones(K,1);
u = zeros(K,T);
g = cell(K,T);
options = optimoptions('fmincon','MaxFunctionEvaluations',1e5,'OptimalityTolerance', 1e-4,'Algorithm','active-set','display','iter');
% for i = 1:K
%     N(i,:) = poissrnd(N_all(i),[1,T]);
%     r_i = R_all(i);
%     for j = 1:T
%         n = N(i,j);
%         h = (randn(n,1)+1j*randn(n,1))/sqrt(2);
%         r = sqrt(rand(n,1).*r_i.^2);
%         g{i,j} = h.*conj(h)./(1+r.^alpha);
%     end
% end
% load Ng_6_100.mat N g
N_S = sum(N);
N_S_S = sum(sum(N));
W_new = zeros(K,T);
for iter = 1:MAX_ITER
    for i = 1: T
        %% x=[p,w,W]
        W0 = W_mean;
        %% cvx based solution
        cvx_begin
        variable w0(N_S(i),1) nonnegative
        variable p0(N_S(i),1) nonnegative
        variable W(K,1) nonnegative
        minimize sum(p0)+rho/2*pow_pos(norm(W-W_mean+u(:,i)),2)
        sum(W) <= 1
         head = 1;
        for k = 1:K
            tail = head + N(k,i)-1;
            sum(w0(head:tail,1)) <= W(k);
            -rel_entr(w0(head:tail),w0(head:tail)+p0(head:tail).*g{k,i}/(Wmax*noise_power))/log(2) >= Rmin(k)/Wmax
            head = tail + 1;
        end
        cvx_end
        W_new(:,i) = W;
        %% Matlab based solution
%         w0 = zeros(N_S(i),1);
%         p0 = zeros(N_S(i),1);
%         A = zeros(K+1,2*N_S(i)+K);
%         b = [zeros(K,1);1];
%         A(end,2*N_S(i)+1:end) = 1;
%         head = 1;
%         for k = 1:K
%             tail = head + N(k,i)-1;
%             A(k,N_S(i)+(head:tail)) = 1;
%             A(k,2*N_S(i)+k) = -1;
%             w0(head:tail) = W0(k)/N(k,i);
%             p0(head:tail) = noise_power*Wmax*(2.^(Rmin(k)./(Wmax*w0(head:tail)))-1).*w0(head:tail)./g{k,i};
%             head = tail + 1;
%         end
%         x0 = [2*p0;w0;W0];
%         sol = fmincon(@(x)sum(x(1:N_S(i)))+rho/2*norm(x(2*N_S(i)+1:end)-W_mean+u(:,i))^2,x0,...
%             A,b,[],[],zeros(2*N_S(i)+K,1),[],@(x) rate_cons(x,N(:,i),Wmax,Rmin,noise_power,g(:,i)),options);
%         W_new(:,i) = sol(2*N_S(i)+1:end);
    end
    
    W_mean_new = mean(W_new,2);
    
    u_new = u+rho*(W_new - W_mean_new);
    
    r = sum(vecnorm(W_new-W_mean_new).^2);
    tol_pri = sqrt(T)*ABSTOL + RELTOL*(max(max(vecnorm(W_new)),norm(W_mean_new)));
    s = T*rho^2*norm(W_mean_new-W_mean).^2;
    tol_dual =  sqrt(T)*ABSTOL+RELTOL*max(vecnorm(u_new));
    history.r(iter) = r;
    history.s(iter) = s;
    history.eps_pri(iter)  = tol_pri;
    history.eps_dual(iter) = tol_dual;
    if r <= tol_pri && s <= tol_dual
        W_final =  W_mean_new;
        break;
    end
    if r>t*s
        rho = rho*(1+tau);
    elseif s>t*r
        rho = rho/(1+tau);
    end
    
    W_mean = W_mean_new;
    u = u_new;
end
end
% x,N(:,i),Wmax,Rmin,noise_power,g{:,i}
function [c, ceq] = rate_cons(x,N,Wmax,Rmin,noise_power,g)
N_S = sum(N);
ceq = [];
% c = [];
K = length(Rmin);
c = zeros(K,1);
head = 1;
for i = 1:K
    tail = head + N(i)-1;
    pi = x(head:tail);
    wi = x(N_S+(head:tail));
    c(i) = -min(Wmax*wi.*log2(1+pi.*g{i}./(Wmax*noise_power*wi)))/Rmin(i)+1;
    head = tail+1;
end
end
% prob = optimproblem('ObjectiveSense','min');
% W = optimvar('W', K,1, 'LowerBound',0);
% P_all = optimvar('P_all', K,T, 'LowerBound',0);
% cell_spectrum = sum(W) <= 1;
% user_spectrum = optimconstr(K,T);
% user_rate = optimconstr(K,T);
% user_power = optimconstr(K,T);
% for j = 1:T
%     for i = 1:K
%         r_i = R_all(i);
%         n = N(i,j);
%         w = optimvar(['w_',int2str(i),int2str(j)], n,1, 'LowerBound',0);
%         p = optimvar(['p_',int2str(i),int2str(j)], n,1, 'LowerBound',0);
%         h = (randn(n,1)+1j*randn(n,1))/sqrt(2);
%         r = sqrt(rand(n,1).*r_i.^2);
%         g = h.*conj(h)./(1+r.^alpha);
%         user_power(i,j) = P_all(i,j) == sum(p);
%         user_spectrum(i,j) = sum(w) <= W(i);
%         user_rate(i,j) = Wmax*w.*log(1+p.*g./(Wmax*g*noise_power))/log(2)>= R(i);
%     end
% end
% prob.Constraints.cell_spectrum = cell_spectrum;
% prob.Constraints.user_spectrum = user_spectrum;
% prob.Constraints.user_power = user_power;
% prob.Constraints.user_rate = user_rate;
% sol = solve(prob);
% P=sum(P_all,2);