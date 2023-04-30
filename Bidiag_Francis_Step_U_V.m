function [U_next, Bi_next, V_next] = Bidiag_Francis_Step_U_V(U, Bi, V)
    [m, n] = size(Bi);
    U_next = U;
    Bi_next = Bi;
    V_next = V;

    for k = 1:n-1
        % Introduce the bulge and chase it out
        % Compute Givens rotation for the left side (U)
        [c, s] = givens(Bi_next(k, k), Bi_next(k+1, k));
        G = [c, s; -s, c];

        % Apply Givens rotation to the left side (U)
        U_next(:, k:k+1) = U_next(:, k:k+1) * G';
        Bi_next(k:k+1, k:min(k+2, n)) = G * Bi_next(k:k+1, k:min(k+2, n));

        % Compute Givens rotation for the right side (V)
        [c, s] = givens(Bi_next(k, k), Bi_next(k, k+1));
        G = [c, -s; s, c];

        % Apply Givens rotation to the right side (V)
        V_next(:, k:k+1) = V_next(:, k:k+1) * G';
        Bi_next(k:min(k+2, m), k:k+1) = Bi_next(k:min(k+2, m), k:k+1) * G;
    end
end
