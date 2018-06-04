W = 1;P=1;r=1;alpha=3;
f = @(x) exp(-x).*gammainc(x.*r^alpha,2/alpha).*gamma(2/alpha)./(W*x.^(2/alpha)+P*x.^(2/alpha+1));
result1 = integral(f,0,Inf);
result_all = zeros(96,1);
for i = 5:100
    result_all(i-4)=nume(i)-result1;
end
figure
plot(5:100,result_all)
function result2 = nume(L)
    W = 1;P=1;r=1;alpha=3;
% L = 100; % Number of points used in Gaussian-Laguerre quadrature
    [x_i, w_i] = pyGaussLaguerre(L); % points and weights in Gaussian-Laguerre quadrature
    a=gammainc(x_i*r^alpha,2/alpha).*gamma(2/alpha);
    %b=alpha*x_i.^(2/alpha)*(R_all.^2)'*noise_power;
    b=W*x_i.^(2/alpha);
    c=P*x_i.^(2/alpha+1);
    result2 = sum(w_i.*a./(b+c));
end