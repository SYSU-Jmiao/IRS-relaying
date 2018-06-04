y=rand(1,300000)*3 + rand(1,300000)*3;
ymin=min(y);
ymax=max(y);
x=linspace(ymin,ymax,20);%将最大最小区间分成20个等分点(19等分),然后分别计算各个区间的个数
yy=hist(y,x);%计算各个区间的个数
yy=yy/length(y);%计算各个区间的个数
bar(x,yy)%画出概率密度分布图

s=0;        
for i=2:length(x)
s=[s,trapz(x([1:i]),yy([1:i]))];    % please remove the " ; "
end
figure;
plot(x,s,x,s,'*')