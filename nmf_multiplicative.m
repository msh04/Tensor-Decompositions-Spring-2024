function [B, C] = nmf_multiplicative(A, j, B0, C0) 
    max_iteration = 100;
    B = B0;
    C = C0;
    for i = 1 : max_iteration
        C = C .* (B.' * A) ./ (B.' * B * C + ones(size(C)) * eps);
        B = B .* (A * C.') ./ (B * C * C.' + ones(size(B)) * eps);
    end
end

