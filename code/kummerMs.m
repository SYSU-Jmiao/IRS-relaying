function [ y ] = kummerMs( s,x )
%KUMMA �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
if x<=0
    y=0;
else
    y = s*x.^(-s).*(gamma(s)*gammainc(x,s));
end
end

