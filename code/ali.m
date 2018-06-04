syms u1 u2 c1 c2 t x c12
%c12 = t*(c1+c2);
%f=x^2/(2*u1*u2);
%f = (2*x-u1)/u2;
%f = (-2*x^2+2*(u1+u2)*x-u1^2-u2^2)/(2*u1*u2);
g=(1-f)*(x-c12);
z=simplify(diff(g,x));
solve(z==0,x)
%subs(g,{t,c1,c2},{1,0,0})