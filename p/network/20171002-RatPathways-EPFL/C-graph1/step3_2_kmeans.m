
% RA, 2017-10-26

% [1] http://www.kyb.mpg.de/fileadmin/user_upload/files/publications/attachments/Luxburg07_tutorial_4488%5b0%5d.pdf

%%

close all;
clear all;

%%

% INPUT
input_file_unnorm_mat = './OUTPUT/UV/column-e-lapeig-unnorm.mat';
input_file_normal_mat = './OUTPUT/UV/column-e-lapeig-normal.mat';

% OUTPUT
run = 'run05';
output_path             = ['./OUTPUT/step3_2_kmeans/' run '/'];
output_path_UV          = [output_path 'UV/'];
output_path_clusterflow = [output_path 'clusterflow/'];
output_path_histograms  = [output_path 'histograms/'];

output_file_clusters    = [output_path_UV 'column-f-clusters.mat'];
output_file_clusterflow = [output_path_clusterflow 'flow_all.png'];

%%

% 0. -- Create output directories

mkdir(output_path);
mkdir(output_path_UV);
mkdir(output_path_clusterflow);
mkdir(output_path_histograms);

%%

% 1. -- Load Laplacian and eigenvectors (normalized)

data_normal = load(input_file_normal_mat);
N = data_normal.number_of_eigenvectors;

W = data_normal.W;
assert(~issparse(W));
E = data_normal.E;
assert(isvector(E));

assert(size(W, 2) == N);
n = size(W, 1);

% Restrict the eigenvalue matrix
E = sparse(diag(E(1:N)));

disp(['n x N = ' num2str(n) ' x ' num2str(N)]);

%%

% 2. -- Normalized spectral clustering according to Ng, Jordan, Weiss [1, p.7]

%%

% 2.0 -- Some consistency checks

% Expect a sparse Laplacian
assert(issparse(data_normal.L));

% Check that W are eigenvectors of L_sym with eigenvalues E
inv_sqrt_D = sparse(diag(1 ./ sqrt(diag(data_normal.L))));
L_sym = inv_sqrt_D * data_normal.L * inv_sqrt_D;
assert(issparse(L_sym));
assert(norm(L_sym * W - W * E, 'inf') < 1e-10);
clear inv_sqrt_D


%%

% 2.1 -- Compute clusters


% Cluster with k-means, for different k

opt = statset('MaxIter', 1000);
for k = (1 : N)
    disp(['k = ' num2str(k)]);
    
    % Only use the first k columns for clustering
    w = W(:, 1:k);

    % Normalize rows to get the matrix T
    T = diag(1 ./ sqrt(sum(w.^2, 2))) * w;

    % Cluster, get cluster index for each element
    I = kmeans(T, k, 'Options', opt);
    
    % Relabel clusters by size
    size = histcounts(I, (0 : k) + 0.5);
    assert(length(size) == k);
    [size, rank] = sort(size, 'descend');
    rank(rank) = (1 : k);
    I = rank(I);
    
    Clustering{k} = struct('I', I, 'size', size, 'k', k, 'opt', opt, 'n', n, 'N', N);
    
    % These variables shadow built-in functions
    clear size rank
end

%%

% 2.2 -- Save clusters

save(output_file_clusters, 'Clustering', '-v7.3');

%%

% 2.3.1 -- "Cluster flow" images

load(output_file_clusters);

% Will hold all "cluster flow" images
Z = [];

for a = 1 : min(10, N)
    
    % Will hold one row of Z
    z = [];
    
    for b = 1 : min(10, N)
        Ca = Clustering{a};
        Cb = Clustering{b};
        
        % Sanity check: same number of points was used for clustering
        assert(Ca.n == Cb.n);
        
        % Who goes where?
        S = full(sparse(Ca.I, Cb.I, 1));
        
        % Normalization
        S = S / max(max(S));
        
        border = zeros(size(S, 1), 2);
        z = [z, [border, S, border]];
    end
    
    % Append the row  z  to  Z  with padding
    border = zeros(2, size(z, 2));
    Z = [Z; [border; z; border]];
end

% Color-invert, magnify, and save Z as image
pixel_mag = 20;
imwrite(kron(1-Z, ones(pixel_mag)), output_file_clusterflow);

%%

% 2.3.2 -- Histograms

load(output_file_clusters);

for C = [Clustering{:}]
    k = C.k;
    I = C.I;
    
    % Compute the sizes of the clusters
    H = histcounts(I, (0 : k) + 0.5);
    
    bar(H); 
    set(gca, 'xticklabel', []);
    axis tight;
    
    filename = [output_path_histograms 'k' num2str(k)];
    saveas(gcf, filename, 'epsc');
    close;
end




