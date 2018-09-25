function [points, weights] = pyGaussLaguerre(n)
%PYGAUSSLAGUERRE �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
tmp = py.numpy.polynomial.laguerre.laggauss(int64(n));
points = cellfun(@double,cell(tmp{1}.tolist))';
weights = cellfun(@double,cell(tmp{2}.tolist))';
end

