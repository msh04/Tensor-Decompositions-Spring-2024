function [G, U1, U2, U3] = my_hooi(X, R)
    % initializing using hosvd
    unfolding1 = tenmat(X, 1);
    [leftSingularVector1, ~, ~] = svd(unfolding1.data);
    U1 = leftSingularVector1(:, 1 : R(1));

    unfolding2 = tenmat(X, 2);
    [leftSingularVector2, ~, ~] = svd(unfolding2.data);
    U2 = leftSingularVector2(:, 1 : R(2));

    unfolding3 = tenmat(X, 3);
    [leftSingularVector3, ~, ~] = svd(unfolding3.data);
    U3 = leftSingularVector3(:, 1 : R(3));
    
    max_iteration = 100;
    for i = 1 : max_iteration
        ten1 = ttm(X,{U2', U3'}, [2, 3]);
        unfolding1 = tenmat(ten1, 1);
        [leftSingularVector1, ~, ~] = svd(unfolding1.data);
        U1 = leftSingularVector1(:, 1 : R(1));

        ten2 = ttm(X, {U1', U3'}, [1, 3]);
        unfolding2 = tenmat(ten2, 2);
        [leftSingularVector2, ~, ~] = svd(unfolding2.data);
        U2 = leftSingularVector2(:, 1 : R(2));

        ten3 = ttm(X, {U1', U2'}, [1, 2]);
        unfolding3 = tenmat(ten3, 3);
        [leftSingularVector3, ~, ~] = svd(unfolding3.data);
        U3 = leftSingularVector3(:, 1 : R(3));
    end
    G = ttm(X, {U1', U2', U3'}, [1, 2, 3]);
end

