MAX_IMP = 1000*MAX_IMPLEMENT;
r_all = zeros(K,MAX_IMP);
r_the = zeros(K,1);
f = @(x) (ei(x)-log(x)-eulergamma)*exp(-x);
for i = 1:K
    
    r_i = R_all(i);
    r_ = r_i^alpha;
    f_i = @(x)1/log(2)*(exp(-x).*((2*r_/(2+alpha))*kummerMs(1+2/alpha,x*r_)+kummerMs(2/alpha,x*r_)).*(Wmax*W(i)*log(1+P(i)*x/(Wmax*W(i)*noise_power))));
    r_the(i) = integral(f_i,0,+inf);
    r_the(i) = double(f(vpa(N_all(i))))*r_the(i);
    parfor imp_i = 1:MAX_IMP
        n = poissrnd(lambda_all(i)*pi*R_all(i).^2);
        w = Wmax*W(i)/n;
        p = P(i)/n;
        d = sqrt(rand(n,1)*r_i^2);
        h = (randn(n,1)+1j*randn(n,1))/sqrt(2);
        g = h.*conj(h)./(1+d.^alpha);
        r_all(i,imp_i) = mean(w*log2(1+p.*g/(w*noise_power)));
    end
end
r_mean = mean(r_all,2);

