% =============================== Output 1 ===============================

% Basis of subspace V
v_cols = [
     1, 1, 3;
    -1, 1, 1;
     0, 0, 0;
     1, 0, 1];
v_basis = orth(v_cols)

% Basis of subspace W
w_cols = [
    1,  1;
    0,  0;
    2, -2;
    1,  0];
w_basis = orth(w_cols)

% Verify that the given vectors are lin. indep.
v_basis_dot = dot(v_basis(:, 1), v_basis(:, 2))
w_basis_dot = dot(w_basis(:, 1), w_basis(:, 2))

% =============================== Output 2 ===============================

% Matrix spanning subspace sum of subspaces V and W
sub_sum_v_w = sub_sum(v_basis, w_basis)

% Matrix spanning intersection of subspaces V and W
sub_intersect_v_w = sub_intersect(v_basis, w_basis)

% =============================== Output 3 ===============================

% Finding the coordinates of x in a given basis for R2
r2_basis = [
    1, 2;
    1, 1];
x = [
    2; 
    1];
z = r2_basis\x

% Verifying that the coordinates are correct
difference_exists = norm(x - mtimes(r2_basis, z)) > 1e-10

% =============================== Output 4 ===============================

% Compute the matrix describing the LT between the given bases
a_std = [
    1,  2,  0, -1;
    0,  1, -1, -2;
    1,  0,  3,  4;
    0, -1,  2,  3;
    0,  0,  2,  2];
r4_basis = [
    1, 1, 1, 1;
    0, 1, 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1];
r5_basis = [
    1,  1, 0,  0, 1;
    1, -1, 0,  0, 0;
    0,  0, 1,  1, 0;
    0,  0, 1, -1, 0;
    0,  0, 0,  0, 1];
a_hat = inv(r5_basis) * a_std * r4_basis                                   % TODO: unsure about expression

% =============================== Output 5 ===============================

% Rank of A
rank_a = rank(a_std)

% dim(ker(A)) = n - rank(A)
dim_ker_a = size(a_std, 2) - rank_a

% Determine injectivity and surjectivity
a_is_injective = dim_ker_a == 0
a_is_surjective = size(orth(transpose(a_std)), 2) == size(a_std, 1)

% =============================== Output 6 ===============================

% First vector b has no solutions
b_1 = [
    1;
    0;
    0;
    0;
    0]
b_1_has_solutions = rank([a_std, b_1]) == rank_a

% Second vector b has infinite solutions
b_2 = [
     1;
     1;
    -2;
    -2;
    -2]
b_2_has_solutions = rank([a_std, b_2]) == rank_a
b_2_solution_1 = a_std\b_2
k = 1;  % k is arbitrary
b_2_solution_2 = k * null(a_std) + b_2_solution_1

% Verify that found solutions work
difference_exists = norm(a_std * b_2_solution_1 - b_2) > 1e-10
difference_exists = norm(a_std * b_2_solution_2 - b_2) > 1e-10


