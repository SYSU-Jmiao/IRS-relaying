Rmin(6)= 2e6;
[W_ave,P_ave,test_iter] = ADMM_MNVO_Allocation([W_ben;P_ben_ave_all(4,:)'],R_all,N_all,alpha, noise_power, Rmin, Wmax);
markers = ['*','d','+','s','x','o'];
l_c = ['The 1st MVNO'; 'The 2nd MVNO'; 'The 3rd MVNO'; 'The 4th MVNO'; 'The 5th MVNO';'The 6th MVNO'];
figure
hold on
ylabel('Bandwidth(MHz)');
xlabel('Iterations');
% iter_i = length(Rmin_all);
[~,max_iter] = size(test_iter.W(i,:));
for i = 1:K
    plot(0:(max_iter-1),test_iter.W(i,:) * 100, ['-k',markers(i)]);
end
legend(l_c)
set(gca,'box','on')