y=rand(1,300000)*3 + rand(1,300000)*3;
ymin=min(y);
ymax=max(y);
x=linspace(ymin,ymax,20);%�������С����ֳ�20���ȷֵ�(19�ȷ�),Ȼ��ֱ�����������ĸ���
yy=hist(y,x);%�����������ĸ���
yy=yy/length(y);%�����������ĸ���
bar(x,yy)%���������ܶȷֲ�ͼ

s=0;        
for i=2:length(x)
s=[s,trapz(x([1:i]),yy([1:i]))];    % please remove the " ; "
end
figure;
plot(x,s,x,s,'*')