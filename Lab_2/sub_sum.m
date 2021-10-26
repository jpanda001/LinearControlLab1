function C = sub_sum(A, B)
    % Determine the basis of the subspace sum
    C = orth([A, B]);
end