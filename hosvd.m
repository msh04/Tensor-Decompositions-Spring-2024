function [U1,U2,U3] = hosvd(T)
    [I1, I2, I3] = size(T);

    T_unfolding1 = zeros(I1, I2*I3);
    k = 1;
    for i = 1 : I3
        for j = 1 : I2
            T_unfolding1(:, k) = T(:, j, i);
            k = k + 1;
        end
    end
    [U, ~, ~] = svd(T_unfolding1);
    U1 = U(:, 1:3);

    T_unfolding2 = zeros(I2, I1*I3);
    k = 1;
    for i = 1 : I3
        for j = 1 : I1
            T_unfolding2(:, k) = T(j, :, i);
            k = k + 1;
        end
    end
    [U, ~, ~] = svd(T_unfolding2);
    U2 = U(:, 1:3);
    
    T_unfolding3 = zeros(I3, I1*I2);
    k = 1;
    for i = 1 : I2
        for j = 1 : I1
            T_unfolding3(:, k) = T(j, i, :);
            k = k + 1;
        end
    end
    [U, ~, ~] = svd(T_unfolding3);
    U3 = U(:, 1:3);

end

