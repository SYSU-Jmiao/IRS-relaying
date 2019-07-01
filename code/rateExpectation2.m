function y=rateExpectation2(W, P,R,N,alpha, noise_power, Wmax)
r_ = R^alpha;
N0 = W*Wmax*noise_power;
if P <=0 || W<=0
    y=0;
else
%% 1
%     y = quadgk(@(x)exp(-x).*((2*r_/(2+alpha))*kummerMs(1+2/alpha,x*r_)+kummerMs(2/alpha,x*r_)).*(W*Wmax*log(1+P*x/N0)),0,Inf)/log(2);
%% 2
% s_1 = @(x) ei(x) - log(x) - double(eulergamma)*ones(size(x));
% h = @(x) kummerMs(2/alpha,x*r_).*exp(-x);
% f = @(x) Wmax*W*P*exp(-N)*s_1(N*h(x))./(N0+P*x);
%% 3
f_h =@(x) (2*r_/(2+alpha))*kummerMs(1+2/alpha,x*r_)+kummerMs(2/alpha,x*r_);
h = @(x) kummerMs(2/alpha,x*r_).*exp(-x);
g = @(x) exp(-N)*(exp(N*h(x))-1);
r = @(x) (W*Wmax*log(1+P*x/N0)/log(2));
f = @(x) f_h(x).*g(x).*r(x);
% h = @(x) (1+2*r_)/(2+alpha)*kummerMs(1+2/alpha,x*r_)./kummerMs(2/alpha,x*r_);
% g = @(x) (exp(-N*kummerMs(1+2/alpha,x*r_))-1);
% f = @(x) exp(-N)*g(x).*h(x).*(W*Wmax*log(1+P*x/N0)/log(2));
y = quadgk(f,0,Inf);
end
end
