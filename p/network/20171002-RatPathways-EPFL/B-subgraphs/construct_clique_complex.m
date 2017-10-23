% Construct the clique complex from an undirected graph.
% Apply only to small graphs.

function [nsimplices, simplices] = construct_clique_complex(G)
 
    nnodes = length(G(1, :));
 
    % Exhaust all cliques.
    nsimplices = 0; simplices = {};
 
    % Nodes.
    for n = 1:nnodes
        simplices{n} = [n];
    end
    nsimplices = nnodes; maxdim = 1;
 
    % Edges.
    [a, b] = find(G > 0);
    subset = find(a < b);
    edges = [a(subset)';b(subset)'];
    nedges = length(edges(1, :));
    if (nedges > 0)
        for n = 1:nedges
            simplices{nsimplices + n} = transpose(edges(:, n));
        end
        nsimplices = nsimplices + nedges;
        maxdim = 2;
    end
 
    % Higher-order cliques.
    % Stop when there are no cliques of a given order.
    maxdim = maxdim + 1; flag = 1;
    while (flag == 1)
     
        % Build cliques on top of lower-order cliques.
        % Only allow the permutations of an increasing order.
        ncliques = 0; cliques = [];
        for n = 1:nsimplices
            if (length(simplices{n}) == (maxdim - 1))
                lowerclique = simplices{n};
                vec = ones(1, nnodes);
                for i = 1:length(lowerclique)
                    vec = vec .* G(lowerclique(i), :);
                end
                subset = find(vec > 0); subset = setdiff(subset, 1:max(lowerclique));
                for i = 1:length(subset)
                    ncliques = ncliques + 1; cliques(ncliques, :) = [lowerclique subset(i)];
                end
             
                % for i=1:nnodes
                % if ((ismember(i,lowerclique)==0)&(sum(G(i,lowerclique)>0)==(maxdim-1))&(i>max(lowerclique)))
                % subset=[lowerclique i];
                % ncliques=ncliques+1; cliques(ncliques,:)=subset;
                % end
                % end
             
            end
        end
     
        % Write cliques to the simplices.
        for n = 1:ncliques
            simplices{nsimplices + n} = cliques(n, :);
        end
        nsimplices = nsimplices + ncliques;
     
        % Debug
        % fprintf('maxdim=%d, ncliques=%d, nsimplices=%d\n',maxdim,ncliques,nsimplices);
     
        % Stop if ncliques=0.
        if (ncliques <= 0)
            flag = 0; % maxdim=maxdim-1;
            % Other increment maxdim by one and proceed.
        else
            maxdim = maxdim + 1;
        end
     
    end
 
 
 
