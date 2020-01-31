% load newData.mat
figure;
p1_all = pow2db(squeeze(mean(p_result_ben(1,:,:),3)))+30;
p2_all = pow2db(squeeze(mean(p_result_eaq(1,:,:),3)))+30;
p3_all = pow2db(squeeze(mean(p_result_ave(1,:,:),3)))+30;
p5_all = pow2db(squeeze(mean(p_result_MC(1,:,:),3)))+30;
xt=1:6;
hold on
ylabel('Power(dBm)');
xlabel('MVNO index');
xticks(1:6)
bar([p1_all',p2_all',p3_all',p5_all'])
legend('Benchmark','User Number','Proposed','MC');
set(gca,'box','on')

figure
hold on
ylabel('Spectrum Ratio');
xlabel('MVNO index');
xticks(1:6)
bar([W_ben,W_eaq,W_ave,W_MC])
legend('Benchmark','User Number','Proposed','MC');
set(gca,'box','on')

p1 = squeeze(mean((squeeze(sum(p_result_ben,2))),2));
p2 = squeeze(mean((squeeze(sum(p_result_eaq,2))),2));
p3 = squeeze(mean((squeeze(sum(p_result_ave,2))),2));
p4 = squeeze(mean(squeeze(sum(p_result_min,2)),2));
p5 = squeeze(mean(squeeze(sum(p_result_MC,2)),2));

figure;
xt=Rmin_all/1e6;
xlabel('R_{min}(Mbps)')
save P_method.mat xt p1 p2 p3 p5
hold on
plot(xt,30+pow2db(p1),'-k*')
plot(xt,30+pow2db(p2),'-ks')
plot(xt,30+pow2db(p3),'-kv')
plot(xt,30+pow2db(p5),'-kx')
ylabel('P_{sum}(dBm)')
legend('Benchmark','User Number','Proposed','MC');
