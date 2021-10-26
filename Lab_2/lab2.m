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
difference = x - mtimes(r2_basis, z)

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


