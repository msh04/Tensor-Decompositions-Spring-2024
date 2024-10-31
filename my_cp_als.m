function [U1, U2, U3] = my_cp_als(T, U1_0, U2_0, U3_0)
    U1 = U1_0;
    U2 = U2_0;
    U3 = U3_0;
    [I1, I2, I3] = size(T);

    T_unfolding1 = zeros(I1, I2*I3);
    k = 1;
    for i = 1 : I3
        for j = 1 : I2
            T_unfolding1(:, k) = T(:, j, i);
            k = k + 1;
        end
    end
    
    T_unfolding2 = zeros(I2, I1*I3);
    k = 1;
    for i = 1 : I3
        for j = 1 : I1
            T_unfolding2(:, k) = T(j, :, i);
            k = k + 1;
        end
    end
    
    T_unfolding3 = zeros(I3, I1*I2);
    k = 1;
    for i = 1 : I2
        for j = 1 : I1
            T_unfolding3(:, k) = T(j, i, :);
            k = k + 1;
        end
    end
    
    max_iteration = 30;
    for i = 1 : max_iteration
        U1 = T_unfolding1 * kr(U3,U2) * pinv((U2.'*U2) .* (U3.'*U3));
        U2 = T_unfolding2 * kr(U3,U1) * pinv((U1.'*U1) .* (U3.'*U3));
        U3 = T_unfolding3 * kr(U2,U1) * pinv((U1.'*U1) .* (U2.'*U2));
    end
end

