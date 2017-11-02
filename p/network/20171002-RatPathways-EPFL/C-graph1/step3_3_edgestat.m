
% RA, 2017-11-02

%%

close all;
clear all;

%%

% INPUT

input_file_graph    = '../A-h5-to-txt/OUTPUT/UV/pathways_mc0_Column.h5.mat';

input_path_kmeans   = './OUTPUT/step3_2_kmeans/run*/';
input_file_clusters = 'UV/column-f-clusters.mat';

% Find input_file_clusters for each run
x = dir(strcat(input_path_kmeans, input_file_clusters));
all_input_file_clusters = strcat({x(:).folder}', '/', {x(:).name}')';

% OUTPUT

output_path_edgestat = './OUTPUT/step3_3_edgestat/';
output_file_inoutedg = [output_path_edgestat 'inoutedg.mat'];

%%

% 0. -- Create output directories

mkdir(output_path_edgestat);

%%

% 1.2.1 -- In-cluster edge statistics

M = getfield(load(input_file_graph), 'M');

INCLUSTA = {};
EXPONENT = [];
for file = all_input_file_clusters
    file = file{1};
    
    % Node n belongs to cluster number
    %   Clustering{k}.I(n)
    % out of k clusters
    Clustering = getfield(load(file), 'Clustering');
    
    % Fraction of in-cluster edges for each cluster set
    [A, B] = find(M);
    inclusta = cellfun(@(C)(nnz(C.I(A) == C.I(B)) / length(A)), Clustering);
    
    % How does fraction of in-cluster edges behave wrt # of clusters?
    % Anticipate a relation FOE ~ (# of clusters)^exponent
    exponent = polyder(polyfit(log([1:length(inclusta)]), log(inclusta), 1));
    
    INCLUSTA{end + 1} = inclusta;
    EXPONENT(end + 1) = exponent;
end

save(output_file_inoutedg, '-v7', 'all_input_file_clusters', 'EXPONENT', 'INCLUSTA');

%%

% 1.2.2 -- Plot in-cluster edge statistics

load(output_file_inoutedg);

close all

figure;
markers = 'ox*^sv';

for inclusta = INCLUSTA
    h = loglog(inclusta{1}, [markers(1) '-']);
    hold on;
    
    % Remove the used marker
    markers = markers(2:end);
end

edgestat_exponent.n    = length(EXPONENT);
edgestat_exponent.mean = mean(EXPONENT);
edgestat_exponent.std  = std(EXPONENT);

disp('Slope:');
disp(edgestat_exponent);
