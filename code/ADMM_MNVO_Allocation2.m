function varargout = ADMM_MNVO_Allocation2(x0,R_all,N_all,alpha, noise_power, Rmin, Wmax)
%% initialization
K=0.5*length(x0);
MAX_ITER = 2e3;
ABSTOL   = 1e-6;
RELTOL   = 1e-4;
rho = 1;
tau = 2;
t = 100;
W = x0(1:K);
W_new = zeros(K,1);
P = x0(K+1:end);
P_new = zeros(K,1);
Z = zeros(K,1);
options = optimoptions('fmincon','Algorithm','interior-point','display','off');
miu = zeros(K,1);
history.W(:,1) = W;
%% iteration
for k = 1:MAX_ITER
    for i = 1: K
        WP=fmincon(@(x) x(2)/rho+1/2*norm(x(1)-Z(i)+miu(i))^2,[W(i);P(i)],[],[]...
            ,[],[],[0;0],[1;+Inf],...
            @(x) myfun2(x(1),x(2),R_all(i),N_all(i),alpha, noise_power, Rmin(i), Wmax),options);
        W_new(i) = WP(1);
        P_new(i) = WP(2);        
    end
    
%    %% check
%     W2 = W_new;
%     P2 = P_new;
%     r_=R_all.^alpha;
%     N0 = Wmax*noise_power;
%     for i = 1:K
%         -quadgk(@(x)exp(-x).*((2*r_(i)/(2+alpha))*kummerMs(1+2/alpha,x*r_(i))+kummerMs(2/alpha,x*r_(i))).*(W2(i)*log(1+P2(i)*x/(N0*W2(i)))),0,Inf)/log(2)/(Rmin*N_all(i)/Wmax)
%     end
    
    Z_new = projsplx(W_new+miu);    
    miu_new = miu+(W_new-Z_new);
    
%     r = norm(W_new-Z_new);
%     s = norm(rho*(Z_new-Z));
%     tol_pri = sqrt(K)*ABSTOL+RELTOL*max(norm(W_new), norm(Z_new));
%     tol_dual = sqrt(K)*ABSTOL+RELTOL*norm(miu_new*rho);
%     history.r(k) = r;
%     history.s(k) = s;
%     history.eps_pri(k)  = tol_pri;
%     history.eps_dual(k) = tol_dual;
    
    r = norm((W_new-Z_new)/max(norm(W_new),norm(Z_new)));
    s = norm(rho*(Z_new-Z)/norm(miu_new));
    tol_pri = sqrt(K)*ABSTOL/max(norm(W_new), norm(Z_new))+RELTOL;
    tol_dual = sqrt(K)*ABSTOL/norm(miu_new)+RELTOL;
    history.obj(k) = sum(P_new);
    history.r(k) = r;
    history.s(k) = s;
    history.eps_pri(k)  = tol_pri;
    history.eps_dual(k) = tol_dual;
    history.W(:,k+1) = Z_new;
    if r <= tol_pri && s <= tol_dual
         break;
    end
    if r>t*s
        rho = rho*(1+tau);
    elseif s>t*r
        rho = rho/(1+tau);
    end
    W = W_new;
    P = P_new;
    Z = Z_new;
    miu = miu_new;
end
W = W_new;
P = P_new;
if nargout >= 2
    varargout{1} = W;
    varargout{2} = P;
end

if nargout >= 3
    varargout{3} = history;
end

end
function [c, ceq] = myfun2(W,P,R,N,alpha, noise_power, Rmin, Wmax)
r_ = R^alpha;
ceq = [];
N0 = Wmax*W*noise_power;
% g = @(x) ei(x) - log(x) - eulergamma*ones(size(x));
% f = @(x) 1/log(2).*W*P*exp(-N)*g(N*kummerMs(2/alpha,x*r_).*exp(-x))./(N0+P*x);    
f_h =@(x) (2*r_/(2+alpha))*kummerMs(1+2/alpha,x*r_)+kummerMs(2/alpha,x*r_);
h = @(x) kummerMs(2/alpha,x*r_).*exp(-x);
g = @(x) exp(-N)*(exp(N*h(x))-1);
r = @(x) (W*Wmax*log(1+P*x/N0)/log(2));
f = @(x) f_h(x).*g(x).*r(x);

c = -quadgk(f,0,+Inf)/Rmin+1;
% c = -quadgk(@(x)exp(-x).*((2*r_/(2+alpha))*kummerMs(1+2/alpha,x*r_)+kummerMs(2/alpha,x*r_)).*(W*log(1+P*x/N0)),0,Inf)/log(2)/(Rmin*N/Wmax)+1;
end
