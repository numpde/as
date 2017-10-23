% Test the homology calculation codes on large networks.
% Subsample small networks from the large network.
% Control both node and edge sizes of the sampled subnetworks.
% Swap networks and compute the homology again.

outfilename = './ratcolumn_report.txt';
fout = fopen(outfilename, 'w');

maxnsubnodes = 100; maxnsubedges = 500; nsubsamples = 2;
filenames{1} = '../A-h5-to-txt/ratcolumn_data.txt';

densities = [3.0, 3.5, 4.0, 4.5];

for densityind = 1:length(densities)
 
    density = densities(densityind); 
    maxnsubedges = ceil(maxnsubnodes * density);
 
    for fileind = 1:length(filenames)
     
        % Load the edge list file.
        filename = filenames{fileind};
        edges = load(filename); nedges = length(edges(:, 1));
        edges = edges + 1; nnodes = max(max(edges));
        edges((nedges + 1):(2 * nedges), :) = transpose([edges(:, 2)'; edges(:,1)']);
        nedges = nedges * 2;
        nnodes = max(max(edges));
        edges((nedges + 1):(nedges + nnodes), :) = transpose([1:nnodes; 1:nnodes]);
        nedges = nedges + nnodes;
     
        % Construct the graph.
        G = sparse(edges(:, 1), edges(:, 2), ones(nedges, 1), nnodes, nnodes);
     
        [a, b] = find(G > 0); subset = find(a < b); nedges = length(subset);
     
        % Subsample small networks.
        fprintf('file=%s, maxnsubnodes=%d, maxnsubedges=%d\n', filename, maxnsubnodes, maxnsubedges);
        fprintf(fout, 'file=%s, maxnsubnodes=%d, maxnsubedges=%d\n', filename, maxnsubnodes, maxnsubedges);
     
        for sampleind = 1:nsubsamples
         
            valid = 0;
         
            while (valid == 0)
             
                subG = sample_network3(G, nnodes, maxnsubnodes, maxnsubedges);
                nsubnodes = length(subG(1, :));
                [a, b] = find(subG > 0); 
                subset = find(a < b); 
                nsubedges = length(subset);
             
                % Construct the simplicial complexes.
                [nsimplices, simplices] = construct_clique_complex(subG);
             
                dim = length(simplices{nsimplices});
                vec = zeros(1, dim);
                for i = 1:nsimplices
                    l = length(simplices{i});
                    vec(l) = vec(l) + 1;
                end
                k = 0;
                for i = 2:dim
                    k = max(vec(i) * vec(i - 1), k);
                end
                if (k < 1e6)
                    valid = 1;
                end
             
            end
         
            % Evaluate Betti numbers.
            maxdim = length(simplices{nsimplices}) - 1;
            bs = evaluate_complex_homology(nsimplices, simplices);
         
            % Report the Betti numbers.
             fprintf('ind %d, nnodes=%d, nedges=%d, maxcliquesize=%d, bs=( ', sampleind, nsubnodes, nsubedges, maxdim + 1);
             fprintf(fout, 'ind %d, nnodes=%d, nedges=%d, maxcliquesize=%d, bs=( ', sampleind, nsubnodes, nsubedges, maxdim + 1);
            for i = 1:length(bs)
                 fprintf('%d ', bs(i));
                 fprintf(fout, '%d ', bs(i));
            end
            fprintf('), ');
            fprintf(fout, '), ');
         
            % Randomly swap edges in the network.  Report the Betti numbers of the swapped networks.
            templateG = subG; randbs = zeros(1, length(bs));
            swapG = templateG;
            % for n=1:(nsubedges*10)
            for n = 1:(nsubedges * 5)
                flag = 0;
                [a, b] = find(swapG > 0);
                tmp = find(a < b);
                swapedges = transpose([a(tmp)';b(tmp)']); cnt = 0;
                while (flag == 0)
                    vec = randperm(length(swapedges(:, 1))); pair = vec(1:2); cnt = cnt + 1;
                    if (length(intersect(swapedges(pair(1), :), swapedges(pair(2), :))) == 0)
                        v1 = swapedges(pair(1), 1); v2 = swapedges(pair(1), 2);
                        v3 = swapedges(pair(2), 1); v4 = swapedges(pair(2), 2);
                        if ((swapG(v1, v4) == 0) & (swapG(v2, v3) == 0))
                            oldswapG = swapG;
                            swapG(v1, v2) = 0; swapG(v2, v1) = 0;
                            swapG(v3, v4) = 0; swapG(v4, v3) = 0;
                            swapG(v1, v4) = 1; swapG(v4, v1) = 1;
                            swapG(v2, v3) = 1; swapG(v3, v2) = 1;
                            flag = 1;
                        end
                    end
                    if (cnt >= 100)
                        flag = 1;
                    end
                end
            end
            templateG = swapG;
            [nrandsimplices, randsimplices] = construct_clique_complex(swapG);
            randmaxdim = length(randsimplices{nrandsimplices}) - 1;
            tmpbs = evaluate_complex_homology(nrandsimplices, randsimplices);
            randbs(1:length(tmpbs)) = tmpbs;
            fprintf('randmaxcliquesize=%d, randbs=( ', randmaxdim + 1);
            fprintf(fout, 'randmaxcliquesize=%d, randbs=( ', randmaxdim + 1);
            for i = 1:length(bs)
                fprintf('%d ', randbs(i));
                fprintf(fout, '%d ', randbs(i));
            end
            fprintf(')\n');
            fprintf(fout, ')\n');
         
        end
     
    end
 
end

fclose(fout);

% Extract the summary statistics of the sampling results.
clear nnodesmeans nnodesstds nedgesmeans nedgesstds maxcliquesizemeans maxcliquesizestds bsmeans bsstds;
clear randmaxcliquesizemeas randmaxcliquesizestds randbsmeans randbsstds;
fp = fopen(outfilename);
for densityind = 1:length(densities)
    density = densities(densityind); maxnsubedges = ceil(maxnsubnodes * density);
    for fileind = 1:length(filenames)
        s = fgets(fp); nnodesarray = zeros(1, nsubsamples); nedgesarray = zeros(1, nsubsamples);
        maxcliquesizearray = zeros(1, nsubsamples); bsarray = zeros(20, nsubsamples);
        randmaxcliquesizearray = zeros(1, nsubsamples); randbsarray = zeros(20, nsubsamples);
        for sampleind = 1:nsubsamples
            s = fgets(fp);
            str = 'nnodes='; i1 = strfind(s, str); i1 = i1(1); i1 = i1 + length(str);
            str = ','; i2 = strfind(s, str); i2 = i2(find(i2 > i1)); i2 = i2(1); i2 = i2 - 1;
            str = s(i1:i2); val = atof2(str); nnodesarray(sampleind) = val;
            str = 'nedges='; i1 = strfind(s, str); i1 = i1(1); i1 = i1 + length(str);
            str = ','; i2 = strfind(s, str); i2 = i2(find(i2 > i1)); i2 = i2(1); i2 = i2 - 1;
            str = s(i1:i2); val = atof2(str); nedgesarray(sampleind) = val;
            str = 'maxcliquesize='; i1 = strfind(s, str); i1 = i1(1); i1 = i1 + length(str);
            str = ','; i2 = strfind(s, str); i2 = i2(find(i2 > i1)); i2 = i2(1); i2 = i2 - 1;
            str = s(i1:i2); val = atof2(str); maxcliquesizearray(sampleind) = val;
            str = 'bs='; i1 = strfind(s, str); i1 = i1(1); i1 = i1 + length(str);
            str = ','; i2 = strfind(s, str); i2 = i2(find(i2 > i1)); i2 = i2(1); i2 = i2 - 1;
            str = s(i1:i2);
            n = maxcliquesizearray(sampleind) - 1;
            for i = 1:n
                item = getitemval4(str, i, ' '); val = atof2(item);
                bsarray(i, sampleind) = val;
            end
            str = 'randmaxcliquesize='; i1 = strfind(s, str); i1 = i1(1); i1 = i1 + length(str);
            str = ','; i2 = strfind(s, str); i2 = i2(find(i2 > i1)); i2 = i2(1); i2 = i2 - 1;
            str = s(i1:i2); val = atof2(str); randmaxcliquesizearray(sampleind) = val;
            str = 'randbs='; i1 = strfind(s, str); i1 = i1(1); i1 = i1 + length(str);
            % str=','; i2=strfind(s,str); i2=i2(find(i2>i1)); i2=i2(1); i2=i2-1;
            i2 = length(s) - 1;
            str = s(i1:i2);
            n = randmaxcliquesizearray(sampleind) - 1;
            for i = 1:n
                item = getitemval4(str, i, ' '); val = atof2(item);
                randbsarray(i, sampleind) = val;
            end
        end
        nnodesmeans(densityind, fileind) = mean(nnodesarray);
        nnodesstds(densityind, fileind) = std(nnodesarray);
        nedgesmeans(densityind, fileind) = mean(nedgesarray);
        nedgesstds(densityind, fileind) = std(nedgesarray);
        maxcliquesizemeans(densityind, fileind) = mean(maxcliquesizearray);
        maxcliquesizestds(densityind, fileind) = std(maxcliquesizearray);
        for n = 1:20
            bsmeans(densityind, fileind, n) = mean(bsarray(n, :));
            bsstds(densityind, fileind, n) = std(bsarray(n, :));
        end
        randmaxcliquesizemeans(densityind, fileind) = mean(randmaxcliquesizearray);
        randmaxcliquesizestds(densityind, fileind) = std(randmaxcliquesizearray);
        for n = 1:20
            randbsmeans(densityind, fileind, n) = mean(randbsarray(n, :));
            randbsstds(densityind, fileind, n) = std(randbsarray(n, :));
        end
    end
end
fclose(fp);

% Report the summary statistics.
statisticsfilename = './ratcolumn_report_statistics.csv';
fout = fopen(statisticsfilename, 'w');
fprintf(fout, 'dataset%cnnodes%cnedges%cmaxcliquesize%cb0%cb1%cb2%cb3%cb4%crandmaxcliquesize%crandb0%crandb1%crandb2%crandb3%crandb4\n', 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9);
for densityind = 1:length(densities)
    for fileind = 1:length(filenames)
        str = filenames{fileind}; subset = find(str == '/'); k = subset(length(subset)); str = str((k + 1):length(str));
        fprintf(fout, '%s%c%.0f%c%.2f(+-%.2f)%c', str, 9, nnodesmeans(densityind, fileind), 9, nedgesmeans(densityind, fileind), nedgesstds(densityind, fileind), 9);
        fprintf(fout, '%.2f(+-%.2f)%c', maxcliquesizemeans(densityind, fileind), maxcliquesizestds(densityind, fileind), 9);
        for n = 1:5
            fprintf(fout, '%.2f(+-%.2f)%c', bsmeans(densityind, fileind, n), bsstds(densityind, fileind, n), 9);
        end
        fprintf(fout, '%.2f(+-%.2f)%c', randmaxcliquesizemeans(densityind, fileind), randmaxcliquesizestds(densityind, fileind), 9);
        for n = 1:5
            fprintf(fout, '%.2f(+-%.2f)%c', randbsmeans(densityind, fileind, n), randbsstds(densityind, fileind, n), 9);
        end
        fprintf(fout, '\n');
    end
    fprintf(fout, '\n');
end
fclose(fout);


