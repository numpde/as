
% RA, 2017-10-26

% INPUT
input_filename_mat  = './OUTPUT/UV/column-d-laplacian.mat';

% PARAMETERS
number_of_eigenvectors = 5;

% OUTPUT
output_filename_mat = './OUTPUT/UV/column-e-lapeig.mat';


% 1. Load the Laplacian
L = load(input_filename_mat); 
L = L.L;

% 2. Compute first few eigenvalues
[V, E] = eig(full(L));

V = V(:, 1:number_of_eigenvectors);
E = diag(E);

% 3. Save
save(output_filename_mat, 'L', 'V', 'E', '-v7.3');
