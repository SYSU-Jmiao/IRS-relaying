function [ y ] = kummerMs( s,x )
%KUMMA 此处显示有关此函数的摘要
%   此处显示详细说明
if x<=0
    y=0;
else
    y = s*x.^(-s).*(gamma(s)*gammainc(x,s));
end
end

