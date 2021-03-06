function C = sub_intersect(A, B)
    % Find the nullspace of [A, -B]
    nullspace_basis = null([A, -B]);
    % Extract the part of the nullspace basis corresponding to A these form
    % a set of coefficients for columns of A, -> take as many rows as there
    % are columns of A
    nullspace_basis_a = nullspace_basis(1:size(A, 2), :);
    % Multiply A by the extracted part to get a basis for the intersection
    C = mtimes(A, nullspace_basis_a);
end