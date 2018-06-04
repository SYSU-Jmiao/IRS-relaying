function [points, weights] = pyGaussLaguerre(n)
%PYGAUSSLAGUERRE 此处显示有关此函数的摘要
%   此处显示详细说明
tmp = py.numpy.polynomial.laguerre.laggauss(int64(n));
points = cellfun(@double,cell(tmp{1}.tolist))';
weights = cellfun(@double,cell(tmp{2}.tolist))';
end

