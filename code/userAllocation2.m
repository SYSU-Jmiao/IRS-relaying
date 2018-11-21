%% user resource allocation2 £¨fmincon£©
function [w, p] = userAllocation2(wmax,h,x0)
global options
options.Display = 'iter';
options.ConstraintTolerance = 1e-5;
options.OptimalityTolerance = 1e-5;
options.PlotFcn = @optimplotfval;
n = length(h);
wp=fmincon(@(x) sum(x(n+1:2*n)),x0,[],[],...
    [ones(1,n),zeros(1,n)],wmax,zeros(2*n,1),[],...
    @(x)userConstraint(x,h),options);
w = wp(1:n);
p = wp((n+1):end);
end
function [c, ceq] =userConstraint(wp,h)
global  Wmax noise_power Rmin
n = length(h);
c = [];
w = wp(1:n);
p = wp(n+1:end);
ceq = 1-Wmax*w.*log(1+p.*h./( Wmax*w*noise_power))/log(2)/Rmin;

end
