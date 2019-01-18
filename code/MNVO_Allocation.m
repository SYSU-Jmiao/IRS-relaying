%% MNVO resource allocation with initial point x0
function [W,P] = MNVO_Allocation(x0,R_all,N_all,alpha, noise_power, Rmin, Wmax)
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
W = WP(1:K);
P = WP(K+1:end);
for i = 1 :K
    r_ = R_all(i)^alpha;
    N0 = Wmax*W(i)*noise_power;
    c(i) = -quadgk(@(x)exp(-x).*((2*r_/(2+alpha))*kummerMs(1+2/alpha,x*r_)+kummerMs(2/alpha,x*r_)).*(W(i)*log(1+P(i)*x/N0)),0,Inf)/log(2)/(Rmin(i)*N_all(i)/Wmax)+1;
end
end

