%% MNVO resource allocation with initial point x0
function [W,P] = MNVO_Allocation2(x0,R_all,N_all,alpha, noise_power, Rmin, Wmax)
K = 0.5*length(x0);
options = optimoptions('fmincon','MaxFunctionEvaluations',1e5,'OptimalityTolerance', 1e-5,'Algorithm','interior-point','display','iter','ConstraintTolerance',1e-10);
WP=fmincon(@(x) sum(x(K+1:2*K)),x0,[],[],...
    [ones(1,K),zeros(1,K)],1,zeros(2*K,1),[],...
    @(WP) cellConstraint(WP,R_all,N_all,alpha, noise_power, Rmin, Wmax),options);
W = WP(1:K);
P = WP((K+1):end);
end

function [c, ceq] =cellConstraint(WP,R_all,N_all,alpha, noise_power, Rmin, Wmax)
K = 0.5*length(WP);
c = zeros(K,1);
ceq = zeros(K,1);
W_all = WP(1:K);
P_all = WP(K+1:end);
for i = 1 :K
    N = N_all(i);
    W = W_all(i);
    P = P_all(i);
    r_ = R_all(i)^alpha;
    N0 = Wmax*W*noise_power;
    %     g = @(x) ei(x) - log(x) - double(eulergamma)*ones(size(x));
    %     f = @(x) W(i)*P(i)*exp(-N_all(i))*g(N_all(i)*kummerMs(2/alpha,x*r_).*exp(-x))./(N0+P(i)*x);

    
    f_h = @(x) (2*r_/(2+alpha))*kummerMs(1+2/alpha,x*r_)+kummerMs(2/alpha,x*r_);
    h = @(x) kummerMs(2/alpha,x*r_).*exp(-x);
    g = @(x) exp(-N)*(exp(N*h(x))-1);
    r = @(x) (W*Wmax*log(1+P*x/N0)/log(2));
    f = @(x) f_h(x).*g(x).*r(x);
    
    c(i) = -quadgk(f,0,Inf)/Rmin(i)+1;
    %     c(i) = -quadgk(@(x)exp(-x).*((2*r_/(2+alpha))*kummerMs(1+2/alpha,x*r_)+kummerMs(2/alpha,x*r_)).*(W(i)*log(1+P(i)*x/N0)),0,Inf)/log(2)/(Rmin*N_all(i)/Wmax)+1;
end
end

