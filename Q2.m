X = randomMat(20, 30);
[jU2, jS2, jV2] = jacobi_svd_2sided(X);
[U, S, V] = svd(X);


function [U, S, V] = jacobi_svd_2sided(A)
    [m1, n1] = size(A);
    m = 0;
    n = 0;
    flag = 0;
    S = A;
    if m1 >= n1
        m = m1;
        n = n1;
        S = A;
        flag = 1;
    else
        m = n1;
        n = m1;
        S = A.';
    end    
    U = eye(m);
    V = eye(n);
    tolerance = 1e-3;
    while off_diagonal(S) > tolerance
        for p = 1 : n
            for q = p + 1 : n
                [c, s] = symmetrization(S(p, p), S(p, q), S(q, p), S(q, q));
                J1 = jMatrix(p, q, c, s, m);
                S = J1 * S;
                U = U * J1.';
                [c, s] = rotation(S, p, q);
                J2 = jMatrix(p, q, c, s, n);
                J4 = jMatrix(p, q, c, s, m);
                S = J4.' * S * J2;
                U = U * J4;
                V = V * J2;
            end    
            for q = n + 1 : m
                [c, s] = non_squared_zero(S(p, p), S(q, p));
                J3 = jMatrix(p, q, c, s, m);
                S = J3 * S;
                U = U * J3.';
            end 
        end
    end
    for i = 1 : n
        if S(i, i) < 0
            S(i, i) = - S(i, i);
            U(:, i) = - U(:, i);
        end    
    end
    if flag == 0
        S = S.';
        temp = U;
        U = V;
        V = temp;
    end        
end


function [c, s] = symmetrization(s_pp, s_pq, s_qp, s_qq)
    t = 0;
    if s_pq ~= s_qp 
        t = (s_qp - s_pq) / (s_pp + s_qq);
    end    
    c = 1 / sqrt(1 + t ^ 2);
    s = c * t;
end


function jacobi = jMatrix(p, q, c, s, n)
    jacobi = eye(n);
    jacobi(p, p) = c;
    jacobi(q, q) = c;
    jacobi(p, q) = s;
    jacobi(q, p) = -s;
end


function off = off_diagonal(A)
    off = 0;
    [m, n] = size(A);
    for i = 1:1:m 
        for j = 1:1:n
            if i ~= j
                off = off + A(i, j) ^ 2;
            end    
        end    
    end
    off = sqrt(off);
end


function [c, s] = rotation(A, p, q)
    if abs(A(p, q)) < 1e-6
        c = 1;
        s = 0;
    else 
        tao = (A(q, q) - A(p, p)) / (2 * A(p, q));
        if tao >= 0
            t = 1 / (tao + sqrt(1 + tao ^ 2));
        else 
            t = 1 / (tao - sqrt(1 + tao ^ 2));
        end
        c = 1 / sqrt(1 + t ^ 2);
        s = t * c;
    end    
end


function [c, s] = non_squared_zero(a_pp, a_pq)
    t = 0;
    if a_pp ~= 0 
        t = a_pq / a_pp;
    end    
    c = 1 / sqrt(1 + t ^ 2);
    s = c * t;
end


function A = randomMat(m, n) 
    A = zeros(m, n);
    for i = 1:1:m
        for j = 1:1:n
            A(i, j) = randsample(1:10, 1);   
        end    
    end    
end