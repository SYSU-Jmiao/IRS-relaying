head = 1;

for k = 1:K
    tail = head + N(k,i)-1;
    p0 = sol(head:tail);
    w0 = sol(N_S(i)+(head:tail));
    W_f(k) = sol(2*N_S(i)+k);
    sum(w0)-W_f(k)
    Rmin(k)-min(Wmax*w0.*log2(1+p0.*g{k,i}./(Wmax*noise_power*w0)))
    head = tail + 1;
end