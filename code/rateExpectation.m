function y=rateExpectation(W, P,R, N,alpha, noise_power, Wmax)
r_ = R^alpha;
N0 = W*Wmax*noise_power;
if P <=0 || W<=0
    y=0;
else
    h = @(x) exp(-x).*((2*r_/(2+alpha))*kummerMs(1+2/alpha,x*r_)+kummerMs(2/alpha,x*r_));
    r = @(x) (W*Wmax*log(1+P*x/N0)/N/log(2));
    f = @(x) h(x).*r(x);
    y = quadgk(f,0,Inf);
end
end
