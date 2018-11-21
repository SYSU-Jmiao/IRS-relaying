function y=rateExpectation(W, P,R,alpha, noise_power, Wmax)
r_ = R^alpha;
N0 = W*Wmax*noise_power;
if P <=0 || W<=0
    y=0;
else
    y = quadgk(@(x)exp(-x).*((2*r_/(2+alpha))*kummerMs(1+2/alpha,x*r_)+kummerMs(2/alpha,x*r_)).*(W*Wmax*log(1+P*x/N0)),0,Inf)/log(2);
end
end
