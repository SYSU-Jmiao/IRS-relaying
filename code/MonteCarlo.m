%% MNVO resource allocation with initial point x0
function [W,P] = MonteCarlo(R_all,N_all,alpha, noise_power, Rmin, Wmax,N,g,T)
K = length(R_all);
flag = 0;
if isempty(N)
    flag = 1;
end
if flag
    N = zeros(K,T);
    for i = 1:K
        N(i,:) = poissrnd(N_all(i),[1,T]);
    end
end
N_S = sum(N);
N_S_S = sum(sum(N));
cvx_begin
variable W(K,1) nonnegative
variable P_all(K,T) nonnegative
variable w(N_S_S) nonnegative
variable p(N_S_S) nonnegative
head = 1;
minimise sum(sum(P_all));
for i = 1:K
    for j = 1:T
        tail = head + N(i,j)-1;
        wij = w(head:tail);
        pij = p(head:tail);
        P_all(i,j) == sum(pij);
        head = tail + 1;
        if flag
            n = N(i,j);
            h = (randn(n,1)+1j*randn(n,1))/sqrt(2);
            r = sqrt(rand(n,1).*R_all(i).^2);
            g{i,j} = h.*conj(h)./(1+r.^alpha);
        end
        sum(wij) <= W(i)
        -rel_entr(wij,wij+pij.*g{i,j}/(Wmax*noise_power))/log(2) >= Rmin(i)/Wmax
    end
end
sum(W) <= 1;
cvx_end
P = sum(P_all./N,2);
end