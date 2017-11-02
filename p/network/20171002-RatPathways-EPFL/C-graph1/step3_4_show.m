
% RA, 2017-11-02

%%

close all;
clear all;

%%

run = 'run05';

% INPUT
input_path_kmeans   = ['./OUTPUT/step3_2_kmeans/' run '/'];
input_file_clusters = [input_path_kmeans 'UV/column-f-clusters.mat'];

input_file_nodesxyz = ['../A-h5-to-txt/OUTPUT/UV/pathways_mc0_Column.h5.mat'];

% OUTPUT

%%

% 0. -- Load data

XYZ        = getfield(load(input_file_nodesxyz), 'XYZ');
M          = getfield(load(input_file_nodesxyz), 'M');
Clustering = getfield(load(input_file_clusters), 'Clustering');

%%

% 1.1 -- Plot clusters

C = Clustering{3};
k = C.k;
I = C.I;

close all

for c = (1 : k)
    figure;
    
    h = plot3(XYZ(I == c, 1), XYZ(I == c, 2), XYZ(I == c, 3), '.');
    h.MarkerSize = 1;
    view(1, -90);
    axis equal
    axis off
end

%%



%%

% 1.3 -- Plot edges

C = Clustering{5};
k = C.k;
I = C.I;

close all

[A, B] = find(M);
nedges = length(A);

AB = [A, B];
AB = AB(randperm(nedges), :);

for e = (1 : nedges)
    edge = AB(e, :);
    if (I(edge(1)) == I(edge(2)))
        % In the same cluster
        c = 'r';
        alpha = 0.2;
    else
        % Between different clusters
        c = 'b';
        alpha = 0.2;
    end
        
    h = plot3(XYZ(edge, 1), XYZ(edge, 2), XYZ(edge, 3), [c '-']);
    h.LineWidth = 1;
    h.Color = [h.Color(1:3), alpha];
    
    view(1, -90);
    axis equal
    axis off
    hold on
    pause(0.01);
    
    disp(['Got: ' num2str(round(e / nedges * 100, 1)) '%']);
end
