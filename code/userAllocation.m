%% user resource allocation £¨cvx£©
function [w, p] = userAllocation(W,h)
global Wmax noise_power R_min
n = length(h);
cvx_begin

cvx_precision best
cvx_quiet true
variable w(n,1)
variable p(n,1)
minimize sum(p)
subject to
sum(w) <= W
for k = 1:n
    -rel_entr(w(k),w(k)+p(k)*h(k)/(Wmax*noise_power))/log(2)>=R_min/Wmax
end
w>=0
p>=0
cvx_end
end
