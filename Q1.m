X = randomSym(10);
[V, D] = eig(X);
[jV, jD] = Jacobi_eig(X);
[jV_prime, jD_prime] = Jacobi_eig_cyclic(X);

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


function [p, q] = maxIndexFinder(A)
    [m, n] = size(A);
    max = 0;
    p = 0;
    q = 0;
    for i = 1:1:m 
        for j = 1:1:n
            if i ~= j
                if abs(A(i, j)) > max
                    max = abs(A(i, j));
                    p = i;
                    q = j;
                end
            end
        end
    end 
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


function froNorm = fNorm(A)
    [m, n] = size(A);
    froNorm = 0;
    for i = 1:1:n
        for j = 1:1:m
            froNorm = froNorm + A(i, j) ^ 2;
        end    
    end 
    froNorm = sqrt(froNorm);
end


function [jV, jD] = Jacobi_eig(A)
    tolerance = 1e-6;
    [m, n] = size(A);
    jV = eye(n);
    jD = A;
    while off_diagonal(jD) > tolerance * fNorm(jD)
        [p, q] = maxIndexFinder(jD);
        if (p ~= 0)
            [c, s] = rotation(jD, p, q);
            jacobi = jMatrix(p, q, c, s, n);
            jD = jacobi.' * jD * jacobi;
            jV = jV * jacobi;
        else
            break
        end
    end    
end


function [jV, jD] = Jacobi_eig_cyclic(A)
    tolerance = 1e-6;
    [m, n] = size(A);
    jV = eye(n);
    jD = A;
    while off_diagonal(jD) > tolerance * fNorm(jD)
        for p = 1 : n - 1
            for q = p + 1 : n
                [c, s] = rotation(jD, p, q);
                jacobi = jMatrix(p, q, c, s, n);
                jD = jacobi.' * jD * jacobi;
                jV = jV * jacobi;
            end    
        end    
    end 
end

function A = randomSym(n) 
    A = zeros(n, n);
    for i = 1:1:n
        for j = 1:1:n
            if j >= i
                A(i, j) = randsample(1:10, 1);
                A(j, i) = A(i, j);
            end    
        end    
    end    
end