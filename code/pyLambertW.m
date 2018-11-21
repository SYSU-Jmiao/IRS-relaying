function y = pyLambertW(x)
x0 = x;
if iscolumn(x)
    x = x';
end
if length(x)>1
    y=real(cellfun(@double,cell(py.scipy.special.lambertw(x).tolist)));
else
    y=real(py.scipy.special.lambertw(x));
end
if iscolumn(x0)
    y = y';
end
end

