
% RA, 2017-10-26

%%

% INPUT
input_filename_mat  = './OUTPUT/UV/column-d-laplacian.mat';

% PARAMETERS
number_of_eigenvectors = 5;

% OUTPUT
output_filename_unnorm_mat = './OUTPUT/UV/column-e-lapeig-unnorm.mat';
output_filename_normal_mat = './OUTPUT/UV/column-e-lapeig-normal.mat';

%%

% 1. Load the Laplacian
L = load(input_filename_mat); 
L = L.L;

%%

% 2.1.a Compute the eigenvalues (unnormalized)
[V, E] = eig(full(L));
E = diag(E);

% 2.1.b Save
save(output_filename_unnorm_mat, 'L', 'V', 'E', '-v7.3');
clear V E

%%

% 2.2.a Compute the eigenvalues (normalized)
d = sparse(diag(sqrt(diag(L))));
D = d * d;
[U, E] = eig(full(L), full(D));
W = d * U;
E = diag(E);

% U are eigenvectors of (D \ L)
% W are eigenvectors of (d \ L / d)

% 2.2.b Save
save(output_filename_normal_mat, 'L', 'U', 'W', 'E', '-v7.3');
clear U W E

