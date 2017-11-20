function b = betti(file)
    G = getfield(load(file), 'G');
    [nsimplices, simplices] = construct_clique_complex(G);
    b = evaluate_complex_homology(nsimplices,simplices);
    save(file, 'G', 'b');
end


%%% CODE BY C-H YEANG %%%


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
end


% Evaluate the Betti numbers of a simplicial complex.

function bs = evaluate_complex_homology(nsimplices, simplices)

    nnodes=0;
    for n=1:nsimplices
        nnodes=max(nnodes,max(simplices{n}));
    end

    % Construct the boundary maps.
    % rho_{k}.
    maxdim=length(simplices{nsimplices})-1;
    clear bds rs;

    % Debug
    %for k=1:min(3,maxdim)

    for k=1:maxdim
        nlows=0; lows=[]; nhighs=0; highs=[];
        for n=1:nsimplices
            if (length(simplices{n})==k)
                nlows=nlows+1; lows(nlows,:)=simplices{n};
            elseif (length(simplices{n})==(k+1))
                nhighs=nhighs+1; highs(nhighs,:)=simplices{n};
            end
        end

        %bd=zeros(nhighs,nlows);
        %for n=1:nhighs
        %for m=1:nlows
        %if (sum(ismember(lows(m,:),highs(n,:))==1)==k)
        %i=setdiff(highs(n,:),lows(m,:)); i=find(highs(n,:)==i);
        %c=(-1).^(i-1);
        %bd(n,m)=c;
        %end
        %end
        %end
        %bds{k}=bd; rs(k)=nlows;

        bd=zeros(nhighs,nlows);
        tmp1=zeros(nhighs,nnodes);
        tmp2=zeros(nlows,nnodes);
        for n=1:nhighs
            tmp1(n,highs(n,:))=1;
        end
        for n=1:nlows
            tmp2(n,lows(n,:))=1;
        end
        bd=tmp1*transpose(tmp2);
        bd(find(bd<k))=0; bd(find(bd==k))=1;
        [a,b]=find(bd==1);
        for j=1:length(a)
            n=a(j); m=b(j);
            i=setdiff(highs(n,:),lows(m,:)); i=find(highs(n,:)==i);
            c=(-1).^(i-1);
            bd(n,m)=c;
        end
        bds{k}=bd; rs(k)=nlows;


        % Debug
        %fprintf('%d %d %d\n',k,length(bds{k}(:,1)),length(bds{k}(1,:)));

    end

    % Evaluate Betti numbers.
    bs=zeros(1,maxdim);

    % Debug
    %for k=0:min(2,maxdim-1)

    for k=0:(maxdim-1)

        % Debug
        % Do not calculate the k-dimensional Betti number if the prior two consecutive Betti numbers are 0.
        if ((k>=2)&(bs(k)==0)&(bs(k-1)==0))
            noaction=1;
        else
            noaction=0;
        end

        if (noaction==0)

            if (k==0)
                d=rs(k+1);
            else
                %d=rs(k)-rank(transpose(bds{k}));
                tmp=null(transpose(bds{k}));
                d=length(tmp(1,:));
            end
            bs(k+1)=d-rank(transpose(bds{k+1}));

            % Debug
            %fprintf('bs(%d)=%d\n',k,bs(k+1));

        end

    end

    %for k=1:(maxdim-1)
    %bs(k)=rs(k)-rank(bds{k+1})-rank(bds{k});
    %end

end
