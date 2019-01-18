close all
clear all
global W P r alpha N noise_power point_all weight_all
W = 1e8;
P=db2pow(-10-30);
r=100;
alpha=3.76;
N=500;
noise_power = db2pow(-174-30);
Gamma = db2pow(9.8+8+15.3-10);
noise_power = noise_power*Gamma;
lambda_all = 2e3./(1e6);
n = pi*r.^2.*lambda_all;
W = W/n;
P = P/n;
point_all = dlmread('500_point.txt');
weight_all = dlmread('500_weight.txt');
%%
MAX_IMPLEMENT = 1e8; 
R_simulation = zeros(1,MAX_IMPLEMENT);
h = (randn(MAX_IMPLEMENT,1)+1j*randn(MAX_IMPLEMENT,1))/sqrt(2);
% figure
% histogram(h.*conj(h),'Normalization','cdf')
% hold on
% fh=@(t) 1-exp(-t);
% 
% fplot(fh,[0,max(h.*conj(h))],'-ro')

d = rand(MAX_IMPLEMENT,1)*r^2;
g = h.*conj(h)./(sqrt(d).^alpha);

rate_all = W*log2(1+P.*g/(W*noise_power));
R_simulation_average = sum(rate_all)/MAX_IMPLEMENT;

max_g=max(g);
figure
histogram(g,'Normalization','cdf')
hold on
ft=@(t) 1-2/(alpha).*(gamma(2/alpha)*(gammainc(t*r^alpha,2/alpha))).*(t*r^alpha).^(-2/alpha);
fplot(ft,[0,max_g],'-ro')

%%
global f
f = @(t) (t.^(-2/alpha)).*(gamma(2/alpha)*gammainc(t,2/alpha)).*(2*W*P./(alpha*r^alpha*noise_power*W+alpha*P*t));
% figure
% fplot(f,[0,2*r^alpha])
%f = @(x)x.^(-2/alpha).*(gamma(2/alpha)*gammainc(x*r^alpha,2/alpha)).*(2*W*P./(alpha*r^2*(noise_power*W+P*x)));
[result1] = integral(f,0,Inf)/log(2);
N_points = 450:10:N;
result_HRGH = zeros(length(N_points),1);
% result_GL = zeros(length(N_points),1);
result_logGL = zeros(length(N_points),1);
result_nHRGH = zeros(length(N_points),1);
result_logHRGH = zeros(length(N_points),1);
for i = 1:length(N_points)
    result_HRGH(i) = HRGH(N_points(i));
    result_nHRGH(i) = nHRGH(N_points(i));
    result_logHRGH(i) = logHRGH(N_points(i));
%     result_GL(i) = GL(N_points(i));
    result_logGL(i) = logGL(N_points(i));
end
result_HRGH = result_HRGH/log(2);
result_nHRGH = result_nHRGH/log(2);
result_logHRGH = result_logHRGH/log(2);
% result_GL = result_GL/log(2);
result_logGL = result_logGL/log(2);
figure
plot(N_points,result1*ones(length(N_points),1),'*')
hold on
plot([N_points(1),N],R_simulation_average*ones(2,1),'-o')
plot(N_points,result_HRGH,'-s')
plot(N_points,result_nHRGH,'-+')
plot(N_points,result_logHRGH,'-v')
% plot(N_points,result_GL,'-')
plot(N_points,result_logGL,'-')
legend('GK','simulation','HRGH','nHRGH','logHRGH','logGL')

%%
function result = HRGH(L)
global W P r alpha noise_power point_all weight_all
%% Gaussian-Laguerre quadrature
% [x_i, w_i] = pyGeneralGaussLaguerre(L);
% b = x_i*r^alpha;
%[x_i, w_i] = cgqf( L,5, 0, 0, 0, r^-alpha); % points and weights in Gaussian-Laguerre quadrature
x_i = point_all(L-1,1:L);
x_i = x_i';
w_i = weight_all(L-1,1:L);
w_i = w_i';
nt = L;
mlt = ones(nt,1);
ndx = 1:nt;
[x_i,w_i] = scqf(nt,x_i,mlt,w_i,nt,ndx,6,0,0,0,r^-alpha);

% b = x_i.^2;
% a = 2/alpha.*((b).^(-2/alpha));
% c = (gamma(2/alpha)*gammainc(b,2/alpha));
% d = (gamma(2/alpha+1)*gammainc(b,2/alpha+1));
% e = W*log(1+P*b/(W*noise_power*r^alpha));
% node = (a.*(c+b.^-1.*d*r^alpha)).*e.*r^-alpha.*x_i*2;

%% General Gaussian-Laguerre quadrature
%      ( order, kind, alpha, beta, a, b )
%      5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
%      [x_i, w_i] = cgqf( L,5, -2/alpha, 0, 0, 1+r^-alpha); % points and weights in Gaussian-Laguerre quadrature
%      node = exp(x_i).*(x_i).^(-2/alpha).*(gamma(2/alpha)*gammainc(x_i,2/alpha)).*(2*W*P./(alpha*r^alpha*noise_power*W+alpha*P*x_i));
%node = exp(x_i).*(gamma(2/alpha)*gammainc(x_i,2/alpha)).*(2*W*P./(alpha*r^alpha*noise_power*W+alpha*P*x_i));
%node = r^alpha.^(-2/alpha).*(gamma(2/alpha)*gammainc(x_i*r^alpha,2/alpha)).*(W*P./(noise_power*W+P*x_i));
%.*(x_i*r^alpha).^(-2/alpha)
%% General Gaussian-Hermitian quadrature
%[x_i, w_i] = cgqf( L,6, 0, 0, 0, 1 ); % points and weights in Gaussian-Hermit quadrature
%6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
%     xw = gauss(L,r_hrhermite(L,1e10,1e-12));
%     x_i = xw(:,1);
%     w_i = xw(:,2);
%% General Half-Range-Gaussian-Hermitian quadrature
% x_i = point_all(L-1,1:L);
% x_i = x_i';
% w_i = weight_all(L-1,1:L);
% w_i = w_i';
% r_ = r^alpha;
% snr = W*noise_power/P;
r_ = 1;
snr = W*noise_power*r^alpha/P;
% nt = L;
% mlt = ones(nt,1);
% ndx = 1:nt;
% [x_i,w_i] = scqf(nt,x_i,mlt,w_i,nt,ndx,6,0,0,0,r^-alpha);
node = exp(x_i*r^(-alpha)).*(2/alpha).*((x_i.^2*r_).^(-2/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2*r_,2/alpha)).*(2*W*x_i./(x_i.^2+snr));

%node = (2/alpha).*((x_i.^2).^(-2/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2,2/alpha)).*(2*W*x_i./(x_i.^2+snr));

%node = 2*((x_i).^(1-4/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2*r^alpha,2/alpha)).*(2*W*P/(alpha*r^2))./(P*x_i.^2+noise_power*W);
%node = 2*exp((1-1/r^alpha)*x_i.^2).*(x_i.^(1-4/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2,2/alpha)).*(2*W*P/alpha./(P*x_i.^2+noise_power*W*r^alpha));
%node = 2*(x_i.^(1-4/alpha)).*((r^alpha/(1+r^alpha)).^(-2/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2*(r^alpha/(1+r^alpha)),2/alpha)).*(2*W*P/alpha./(P*x_i.^2+noise_power*W*(1+r^alpha)));

%[x_i, w_i] = pyGaussLaguerre(L);
%node = exp(x_i).*x_i.^(-2/alpha).*(gamma(2/alpha)*gammainc(x_i*r^alpha,2/alpha)).*(2*W*P/(alpha*r^2))./(P*x_i+noise_power*W);
%node = exp(-x_i*(r^-alpha-1)).*(gamma(2/alpha)*gammainc(x_i,2/alpha)).*(2*W*P./(alpha*r^alpha*noise_power*W+alpha*P*x_i));

%b=alpha*x_i.^(2/alpha)*(R_all.^2)'*noise_power;
%[x_i, w_i] = cgqf(L,5,-2/alpha,0,0,r^-alpha);
result = w_i'*node;
end

function result = GL(L)
global W P r alpha noise_power point_all weight_all
%% Gaussian-Laguerre quadrature
[x_i, w_i] = pyGaussLaguerre(L);
% b = x_i*r^alpha;
%[x_i, w_i] = cgqf( L,5, 0, 0, 0, r^-alpha); % points and weights in Gaussian-Laguerre quadrature
% x_i = point_all(L-1,1:L);
% x_i = x_i';
% w_i = weight_all(L-1,1:L);
% w_i = w_i';
% nt = L;
% mlt = ones(nt,1);
% ndx = 1:nt;
% [x_i,w_i] = scqf(nt,x_i,mlt,w_i,nt,ndx,6,0,0,0,r^-alpha);

% b = x_i.^2;
% a = 2/alpha.*((b).^(-2/alpha));
% c = (gamma(2/alpha)*gammainc(b,2/alpha));
% d = (gamma(2/alpha+1)*gammainc(b,2/alpha+1));
% e = W*log(1+P*b/(W*noise_power*r^alpha));
% node = (a.*(c+b.^-1.*d*r^alpha)).*e.*r^-alpha.*x_i*2;

%% General Gaussian-Laguerre quadrature
%      ( order, kind, alpha, beta, a, b )
%      5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
%      [x_i, w_i] = cgqf( L,5, -2/alpha, 0, 0, 1+r^-alpha); % points and weights in Gaussian-Laguerre quadrature
%      node = exp(x_i).*(x_i).^(-2/alpha).*(gamma(2/alpha)*gammainc(x_i,2/alpha)).*(2*W*P./(alpha*r^alpha*noise_power*W+alpha*P*x_i));
%node = exp(x_i).*(gamma(2/alpha)*gammainc(x_i,2/alpha)).*(2*W*P./(alpha*r^alpha*noise_power*W+alpha*P*x_i));
%node = r^alpha.^(-2/alpha).*(gamma(2/alpha)*gammainc(x_i*r^alpha,2/alpha)).*(W*P./(noise_power*W+P*x_i));
%.*(x_i*r^alpha).^(-2/alpha)
%% General Gaussian-Hermitian quadrature
%[x_i, w_i] = cgqf( L,6, 0, 0, 0, 1 ); % points and weights in Gaussian-Hermit quadrature
%6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
%     xw = gauss(L,r_hrhermite(L,1e10,1e-12));
%     x_i = xw(:,1);
%     w_i = xw(:,2);
%% General Half-Range-Gaussian-Hermitian quadrature
% x_i = point_all(L-1,1:L);
% x_i = x_i';
% w_i = weight_all(L-1,1:L);
% w_i = w_i';
% r_ = r^alpha;
% snr = W*noise_power/P;
% r_ = 1;
% snr = W*noise_power*r^alpha/P;
% nt = L;
% mlt = ones(nt,1);
% ndx = 1:nt;
% [x_i,w_i] = scqf(nt,x_i,mlt,w_i,nt,ndx,6,0,0,0,r^-alpha);
% node = (2/alpha).*((x_i.^2*r_).^(-2/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2*r_,2/alpha)).*(2*W*x_i./(x_i.^2+snr));

%node = (2/alpha).*((x_i.^2).^(-2/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2,2/alpha)).*(2*W*x_i./(x_i.^2+snr));

%node = 2*((x_i).^(1-4/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2*r^alpha,2/alpha)).*(2*W*P/(alpha*r^2))./(P*x_i.^2+noise_power*W);
%node = 2*exp((1-1/r^alpha)*x_i.^2).*(x_i.^(1-4/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2,2/alpha)).*(2*W*P/alpha./(P*x_i.^2+noise_power*W*r^alpha));
%node = 2*(x_i.^(1-4/alpha)).*((r^alpha/(1+r^alpha)).^(-2/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2*(r^alpha/(1+r^alpha)),2/alpha)).*(2*W*P/alpha./(P*x_i.^2+noise_power*W*(1+r^alpha)));

%[x_i, w_i] = pyGaussLaguerre(L);
node = x_i.^(-2/alpha).*(gamma(2/alpha)*gammainc(x_i*r^alpha,2/alpha)).*(2*W*P/(alpha*r^2))./(P*x_i+noise_power*W);
%node = exp(-x_i*(r^-alpha-1)).*(gamma(2/alpha)*gammainc(x_i,2/alpha)).*(2*W*P./(alpha*r^alpha*noise_power*W+alpha*P*x_i));

%b=alpha*x_i.^(2/alpha)*(R_all.^2)'*noise_power;
%[x_i, w_i] = cgqf(L,5,-2/alpha,0,0,r^-alpha);
result = w_i'*node;
end

function result = nHRGH(L)
global W P r alpha noise_power point_all weight_all
%% Gaussian-Laguerre quadrature
% [x_i, w_i] = pyGeneralGaussLaguerre(L);
% b = x_i*r^alpha;
%[x_i, w_i] = cgqf( L,5, 0, 0, 0, r^-alpha); % points and weights in Gaussian-Laguerre quadrature
x_i = point_all(L-1,1:L);
x_i = x_i';
w_i = weight_all(L-1,1:L);
w_i = w_i';
% nt = L;
% mlt = ones(nt,1);
% ndx = 1:nt;
% [x_i,w_i] = scqf(nt,x_i,mlt,w_i,nt,ndx,6,0,0,0,r^-alpha);

% b = x_i.^2;
% a = 2/alpha.*((b).^(-2/alpha));
% c = (gamma(2/alpha)*gammainc(b,2/alpha));
% d = (gamma(2/alpha+1)*gammainc(b,2/alpha+1));
% e = W*log(1+P*b/(W*noise_power*r^alpha));
% node = (a.*(c+b.^-1.*d*r^alpha)).*e.*r^-alpha.*x_i*2;

%% General Gaussian-Laguerre quadrature
%      ( order, kind, alpha, beta, a, b )
%      5, Generalized Laguerre, (a,+oo)     (x-a)^alpha*exp(-b*(x-a))
%      [x_i, w_i] = cgqf( L,5, -2/alpha, 0, 0, 1+r^-alpha); % points and weights in Gaussian-Laguerre quadrature
%      node = exp(x_i).*(x_i).^(-2/alpha).*(gamma(2/alpha)*gammainc(x_i,2/alpha)).*(2*W*P./(alpha*r^alpha*noise_power*W+alpha*P*x_i));
%node = exp(x_i).*(gamma(2/alpha)*gammainc(x_i,2/alpha)).*(2*W*P./(alpha*r^alpha*noise_power*W+alpha*P*x_i));
%node = r^alpha.^(-2/alpha).*(gamma(2/alpha)*gammainc(x_i*r^alpha,2/alpha)).*(W*P./(noise_power*W+P*x_i));
%.*(x_i*r^alpha).^(-2/alpha)
%% General Gaussian-Hermitian quadrature
%[x_i, w_i] = cgqf( L,6, 0, 0, 0, 1 ); % points and weights in Gaussian-Hermit quadrature
%6, Generalized Hermite,  (-oo,+oo)   |x-a|^alpha*exp(-b*(x-a)^2)
%     xw = gauss(L,r_hrhermite(L,1e10,1e-12));
%     x_i = xw(:,1);
%     w_i = xw(:,2);
%% General Half-Range-Gaussian-Hermitian quadrature
% x_i = point_all(L-1,1:L);
% x_i = x_i';
% w_i = weight_all(L-1,1:L);
% w_i = w_i';
r_ = r^alpha;
inv_snr = W*noise_power/P;
% r_ = 1;
% snr = W*noise_power*r^alpha/P;
% nt = L;
% mlt = ones(nt,1);
% ndx = 1:nt;
% [x_i,w_i] = scqf(nt,x_i,mlt,w_i,nt,ndx,6,0,0,0,r^-alpha);
node = exp(x_i.^2).*(2/alpha).*((x_i.^2*r_).^(-2/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2*r_,2/alpha)).*(2.*x_i*W./(x_i.^2+inv_snr));

%node = (2/alpha).*((x_i.^2).^(-2/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2,2/alpha)).*(2*W*x_i./(x_i.^2+snr));

%node = 2*((x_i).^(1-4/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2*r^alpha,2/alpha)).*(2*W*P/(alpha*r^2))./(P*x_i.^2+noise_power*W);
%node = 2*exp((1-1/r^alpha)*x_i.^2).*(x_i.^(1-4/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2,2/alpha)).*(2*W*P/alpha./(P*x_i.^2+noise_power*W*r^alpha));
%node = 2*(x_i.^(1-4/alpha)).*((r^alpha/(1+r^alpha)).^(-2/alpha)).*(gamma(2/alpha)*gammainc(x_i.^2*(r^alpha/(1+r^alpha)),2/alpha)).*(2*W*P/alpha./(P*x_i.^2+noise_power*W*(1+r^alpha)));

%[x_i, w_i] = pyGaussLaguerre(L);
%node = exp(x_i).*x_i.^(-2/alpha).*(gamma(2/alpha)*gammainc(x_i*r^alpha,2/alpha)).*(2*W*P/(alpha*r^2))./(P*x_i+noise_power*W);
%node = exp(-x_i*(r^-alpha-1)).*(gamma(2/alpha)*gammainc(x_i,2/alpha)).*(2*W*P./(alpha*r^alpha*noise_power*W+alpha*P*x_i));

%b=alpha*x_i.^(2/alpha)*(R_all.^2)'*noise_power;
%[x_i, w_i] = cgqf(L,5,-2/alpha,0,0,r^-alpha);
result = w_i'*node;
end
function result = logHRGH(L)
global W P r alpha noise_power point_all weight_all
x_i = point_all(L-1,1:L);
x_i = x_i';
w_i = weight_all(L-1,1:L);
w_i = w_i';
% [x_i, w_i] = pyGaussLaguerre(L);
r_ = r^alpha;
snr = P/(W*noise_power);
node = exp(x_i.^2).*((2*r_/(2+alpha))*kummerMs(1+2/alpha,x_i.^2*r_)).*(W*log(1+x_i.^2*snr)).*x_i*2;

% node = ((2*r_/(2+alpha))*kummerMs(1+2/alpha,x_i*r_)+kummerMs(2/alpha,x_i*r_)).*(W*log(1+x_i*snr));
% node = ((2*r_/(2+alpha))*kummerMs(1+2/alpha,x_i.^2*r_)+kummerMs(2/alpha,x_i.^2*r_)).*(W*log(1+x_i.^2*snr)).*x_i*2;
result = w_i'*node;
end

function result = logGL(L)
global W P r alpha noise_power point_all weight_all
% x_i = point_all(L-1,1:L);
% x_i = x_i';
% w_i = weight_all(L-1,1:L);
% w_i = w_i';
[x_i, w_i] = pyGaussLaguerre(L);
% w_i = w_i';
r_ = r^alpha;
snr = P/(W*noise_power);

node = exp(x_i).*((2*r_/(2+alpha))*kummerMs(1+2/alpha,x_i*r_)).*(W*log(1+x_i*snr));

% node = ((2*r_/(2+alpha))*kummerMs(1+2/alpha,x_i*r_)+kummerMs(2/alpha,x_i*r_)).*(W*log(1+x_i*snr));
% node = ((2*r_/(2+alpha))*kummerMs(1+2/alpha,x_i.^2*r_)+kummerMs(2/alpha,x_i.^2*r_)).*(W*log(1+x_i.^2*snr)).*x_i*2;
result = w_i'*node;
end
