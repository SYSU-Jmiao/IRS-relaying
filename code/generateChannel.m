max_imp = 100000;
max_user_number = 200;
channel_all = zeros(max_user_number,K,max_imp);
for i = 1:max_imp
    for j = 1:K
        r_j = R_all(j);
        h = (randn(max_user_number,1)+1j*randn(max_user_number,1))/sqrt(2);
        r = sqrt(rand(max_user_number,1).*r_j.^2);
        channel_all(:,j,i) = h.*conj(h)./(1+r.^alpha);
    end
end
save rate_channel_all.mat channel_all
