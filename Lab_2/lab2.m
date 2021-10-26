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
