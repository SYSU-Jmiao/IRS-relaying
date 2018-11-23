for i = 1:10
    for j = 1:10
        parfor k =1:10
            i*j*k
        end
    end
end
