function [points, weights] = GaussLaguerre(n)
    syms t
    w(t) = exp(-t);
    F = laguerreL(0:n-1, t);
    X = sym('X', [1, n+1]);
    L = poly2sym(X, t);
    sys = [int(F.*L.*w(t), t, 0, inf) == 0];
    sys = [sys, int(L^2.*w(t), 0, inf) == 1];
    %S = solve(sys, X);
    %structfun(@display, S);
    sys = [sys, X(1)>0];
    S = solve(sys, X);
    L = subs(L, S);
    laguerreL(n, t);
    x = solve(L);
    points = double(real(vpa(x)));
    %isAlways(in(x, 'real'))
    %xradical = solve(L, 'MaxDegree', n);
    y = sym('y', [n, 1]);
    sys = sym(zeros(n));
    for k=0:n-1
        sys(k+1) = sum(y.*(x.^k)) == int(t^k * w(t), t, 0, inf);
    end
    %sys
    %weights = vpa(vpasolve(sys, y));
    weights = double(struct2array(vpasolve(sys, y))');
    %weights structfun(vpa,weights)
    %[alpha1, alpha2, alpha3, alpha4] = solve(sys, y)
    %S = solve(sys, y)
    %structfun(@double, S)
end