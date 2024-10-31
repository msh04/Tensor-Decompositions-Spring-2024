X = randomMat(14, 9);
[jU1, jS1, jV1] = one_sided(X);
[U, S, V] = svd(X);

function jacobi = jMatrix(p, q, c, s, n)
    jacobi = eye(n);
    jacobi(p, p) = c;
    jacobi(q, q) = c;
    jacobi(p, q) = s;
    jacobi(q, p) = -s;
end


function [U, S, V] = one_sided(A) 
    [m1, n1] = size(A);
    m = 0;
    n = 0;
    flag = 0;
    US = A;
    if m1 >= n1
        m = m1;
        n = n1;
        flag = 1;
    else
        m = n1;
        n = m1;
        US = A.';
    end    
    U = zeros(m, m);
    S = zeros(m, n);
    V = eye(n);
    while ~check_orthogonality(US)
        for p = 1 : n - 1
            for q = p + 1 : n
                c = 0;
                s = 0;
                x = US(:, p);
                y = US(:, q);
                a = y.' * x;
                b = y.' * y - x.' * x;
                t = 0;
                if a == 0 
                    c = 1;
                    s = 0;
                else 
                    t = min((- b + sqrt(b ^ 2 + 4 * a ^ 2)) / a, (-b - sqrt(b ^ 2 + 4 * a ^ 2)) / a); 
                end   
                c = 1 / sqrt(1 + t ^ 2);
                s = c * t;
                J = jMatrix(p, q, c, s, n);
                US = US * J;    
                V = V * J;
            end
        end     
    end 
    v = singular_value(US).';
    US = US(:, v);
    V = V(:, v);
    for i = 1 : n
        s_i = norm(US(:, i), 2);
        if s_i ~= 0
            S(i, i) = s_i;
            U(:, i) = US(:, i) / s_i;
        end  
    end
    U = orthogonal(U, rank(U));
    if flag == 0
        S = S.';
        temp = U;
        U = V;
        V = temp;
    end 
end


function bool = check_orthogonality(A)
    bool = false;
    [m, n] = size(A);
    for i = 1 : n
        A(:, i) = A(:, i) / norm(A(:, i), 2);
    end    
    tolerance = 1e-6;
    if norm((A.' * A - eye(n)), "fro") < tolerance
        bool = true;
    end    
end


function v = singular_value(A)
    [m, n] = size(A);
    S = zeros(n, 1);
    v = zeros(n, 1);
    for i = 1 : n
        S(i) = norm(A(:, i), 2);
        v(i) = i;
    end   
    for i = 1 : n - 1
        for j = i + 1 : n
            if S(i) < S(j)
                temp = S(i);
                S(i) = S(j);
                S(j) = temp;
                temp = v(i);
                v(i) = v(j);
                v(j) = temp;
            end    
        end    
    end 
end


function Q = orthogonal(A, r)
    [m, n] = size(A);
    Q = A;
    h = 1;
    for i = r + 1 : n
        q = zeros(m, 1);
        q(h) = 1;
        flag = 0;
        while flag == 0
            for j = 1 : i - 1
                if abs(abs(q.' * Q(:, j)) - norm(q, 2) .* norm(Q(:, j), 2)) < 1e-10
                    flag = 1;
                    break;
                end 
            end
            if flag == 1
                q(h) = 0;
                q(h + 1) = 1;
                h = h + 1;
                flag = 0;
            else 
                u = q;
                for j = 1 : i - 1
                    u = u - (u.' * Q(:, j)) * Q(:, j);
                end    
                Q(:, i) = u / norm(u, 2);
                q(h) = 0;
                q(h + 1) = 1;
                h = h + 1;
                break;
            end    
        end
    end 
end


function A = randomMat(m, n) 
    A = zeros(m, n);
    for i = 1:1:m
        for j = 1:1:n
            A(i, j) = randsample(1:10, 1);   
        end    
    end    
end