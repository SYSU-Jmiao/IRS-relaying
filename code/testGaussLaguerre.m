
global W P r alpha N noise_power
W = 1e9;P=1;r=1e3;alpha=3;N=100;noise_power = db2pow(-60);
f = @(x) exp(-x).*x.^(-2/alpha).*gammainc(x.*r^alpha,2/alpha,'lower').*gamma(2/alpha)./(W*noise_power+P.*x)*2*W*P/(alpha*r^2);
[result1,e] = quadgk(f,0,Inf);
e
result_all = zeros(N-4,1);
for i = 5:N
    result_all(i-4)=nume(i);
end
figure
plot(5:N,result1*ones(N-4,1),':')
hold on
plot(5:N,result_all,'-')
legend('real','simulation')
function result2 = nume(L)
     global W P r alpha noise_power
     % L = 100; % Number of points used in Gaussian-Laguerre quadrature
    [x_i, w_i] = pyGaussLaguerre(L); % points and weights in Gaussian-Laguerre quadrature
    a=vpa(gammainc((x_i.*r^alpha),2/alpha,'lower')).*gamma(2/alpha).*2.*W.*P./(alpha.*r.^2);
    %b=alpha*x_i.^(2/alpha)*(R_all.^2)'*noise_power;
    b=W*x_i.^(2/alpha).*noise_power;
    c=P*x_i.^(2/alpha+1);
    result2 = sum(w_i.*a./(b+c));
end