%% user resource allocation3 closed form 
function [w, p] = userAllocation3(wmax,h,Wmax, noise_power, Rmin)
R = Rmin/Wmax;
N0 = Wmax *noise_power;
n = length(h);
if n == 0
    w=0;
    p=0;
    return
end
e = exp(1);
cf_P = @(W) N0*(2.^(R./(W))-1).*W./h;
cf_W = @(x) log(2)*R*ones(n,1)./(ones(n,1)+lambertw(x/(e*N0)*h-1/e*ones(n,1)));
f = @(x) sum((cf_W(x)))-wmax;
hmin = min(h);
hmax=max(h);

tmp = N0*e*(n*log(2)*R/wmax-1)*exp(n*log(2)*R/wmax-1)+N0;
xmin = tmp/hmax;
xmax = tmp/hmin;

if f(xmin)*f(xmax) <0 
    x = fzero(f,[xmin,xmax]);
else
    f(xmin)
    f(xmax)
    x = fzero(f,(xmin+xmax)/2);
end
w = cf_W(x);
p = cf_P(w);
end
