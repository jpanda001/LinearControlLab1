clear
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

% Verifying that the coordinates are correct -> difference_exists = 0 if
% there is no difference
difference_exists = norm(x - mtimes(r2_basis, z)) > 1e-10

% =============================== Output 4 ===============================

% Compute the matrix describing the LT between the given bases
a_std = [
    1,  2,  0, -1;
    0,  1, -1, -2;
    1,  0,  3,  4;
    0, -1,  2,  3;
    0,  0,  2,  2];

% P 
r4_basis = [
    1, 1, 1, 1;
    0, 1, 0, 0;
    0, 0, 1, 0;
    0, 0, 0, 1];
% Q
r5_basis = [
    1,  1, 0,  0, 1;
    1, -1, 0,  0, 0;
    0,  0, 1,  1, 0;
    0,  0, 1, -1, 0;
    0,  0, 0,  0, 1];
a_hat = inv(r5_basis) * a_std * r4_basis                                   % TODO: unsure about expression - Ben confirmed to be correct

% confirmation:
x = [1; 2; 3; 4];
y_normal = a_std * x;
y_ahat =  r5_basis * a_hat * (r4_basis\x);
difference_exists = norm(y_normal - y_ahat) > 1e-10



% =============================== Output 5 ===============================

% Rank of A
rank_a = rank(a_std)

% dim(ker(A)) = n - rank(A), where n = dimension of domain space = num.cols
dim_ker_a = size(a_std, 2) - rank_a

% Determine injectivity and surjectivity
a_is_injective = dim_ker_a == 0

% surjective if dim(row_space) = number of rows -> rows are lin.indep.
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
b_2_has_infinite_solutions = size(null(a_std), 2) > 0
b_2_solution_1 = a_std\b_2
k = 1;  % k is arbitrary
b_2_solution_2 = k * null(a_std) + b_2_solution_1

% Verify that found solutions work
difference_exists = norm(a_std * b_2_solution_1 - b_2) > 1e-10
difference_exists = norm(a_std * b_2_solution_2 - b_2) > 1e-10

% =============================== Output 7 ===============================
A = [1 2 2;
    1 -3 -5;
    -1 2 4];
V = [0 1;
    1 -2;
    -1 1];
% V is A-invariant if AV is also in the range space of V
V_is_A_invariant = rank([A*V, V]) ==  rank(V)
V_basis = orth(V);
% null(V') is complement basis of V because it contains all elements whose
% dot product with element in V -> complete the basis that V lack
V_compl_basis = null(V');

% P has A-invariant basis on the left, and the other basis on the right
P = [V_basis, V_compl_basis]
A_hat = inv(P) * A * P

n = size(A_hat,1);
k = rank(V);
% matrix is upper triangular if bottom left (n-k) x (k) matrix is zero
is_upper_triangular = all(A_hat(k+1:end, 1:k) < 1e-10 * ones(n-k, k))

% =============================== Output 8 ===============================

A = [5 -11 5;
    0 -6 0;
    -5 5 -1];
B = [1 -2; 
    0 0; 
    1 2];

Qc = ctrb(A,B);

% by contrllability theorem if rank(Qc) = n, for A in R(nxn) ->controllable
AB_is_controllable = rank(Qc) == size(Qc,1);
V = orth(Qc);
W = null(Qc');
P = [V, W];
A_hat = inv(P)*A*P
B_hat = inv(P)*B

n = size(A_hat,1);
k = rank(V);
m = size(B, 2);
% same as textbook if A_hat is upper triangular, and B_hat has zero bottom
same_as_text=all(A_hat(k+1:end,1:k)<1e-10*ones(n-k, k))...
    && all(B_hat(k+1:end,:)<1e-10*ones(n-k,m))

controllable_subsystem_A11 = A_hat(1:k, 1:k)
controllable_subsystem_A12 = A_hat(1:k, k+1:end)
controllable_subsystem_B1 = B_hat(1:k, m)

uncontrollable_subsystem_A22 = A_hat(k+1:end,k+1:end)

