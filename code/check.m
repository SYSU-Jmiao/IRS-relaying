
sum(pow2db(sum((cellfun(@sum,p3)))))/
Wmax*w{2}.*log(1+p{2}.*g{2}./( Wmax*w{2}*noise_power))/log(2)
for i = 1:K
-quadgk(@(x)exp(-x).*((2*r_(i)/(2+alpha))*kummerMs(1+2/alpha,x*r_(i))+kummerMs(2/alpha,x*r_(i))).*(W(i)*log(1+P(i)*x/(N0*W(i)))),0,Inf)/log(2)/(Rmin*N_all(i)/Wmax)
end

for i = 1:K
-quadgk(@(x)exp(-x).*((2*r_(i)/(2+alpha))*kummerMs(1+2/alpha,x*r_(i))+kummerMs(2/alpha,x*r_(i))).*(W2(i)*log(1+P2(i)*x/(N0*W2(i)))),0,Inf)/log(2)/(Rmin*N_all(i)/Wmax)
end