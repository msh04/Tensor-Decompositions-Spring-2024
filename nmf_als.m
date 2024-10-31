function [B, C] = nmf_als(A, j, B0, C0)
    max_iteration = 100;
    B = B0;
    C = C0;
    for i = 1 : max_iteration
         B = max(eps, (A * C.') * pinv(C * C.'));
         C = max(eps, pinv(B.' * B) * (B.' * A));
    end
end

