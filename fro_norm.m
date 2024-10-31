function norm = fro_norm(T)
    norm = 0; 
    [I1, I2, I3] = size(T);
    for i = 1 : I1
        for j = 1 : I2
            for k = 1 : I3
                norm = norm + T(i, j, k) ^ 2;
            end
        end    
    end    
    norm = sqrt(norm);
end