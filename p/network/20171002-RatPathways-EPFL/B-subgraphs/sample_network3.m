% Subsampling from a large network.
% Pick a seed according to the degree distribution.
% Carry out random walks in the network.
% Stop when the number of subsampled nodes reaches the pre-specified number.
% Difference from sample_network2: attempt to control both node and edge sizes.
% Select the seed according to the node degree distribution.
% Incorporate new nodes by sampling neighbors of the current network according to the distribution of the density of the augmented network.
% f(x;rho)=|x-rho|.  g(x;rho)=max(f(x;rho))-f(x;rho)+delta.  p(x;rho)=g(x;rho)/sum_{x}g(x;rho).

function subG = sample_network3(G, nnodes, nsubnodes, nsubedges)


[a,b]=find(G>0); subset=find(a<b); a=a(subset); b=b(subset);
degs=zeros(1,nnodes);
for i=1:length(a)
degs(a(i))=degs(a(i))+1; degs(b(i))=degs(b(i))+1;
end

found=0;

while (found==0)


%vals=find(degs>=10); 
vals=find(degs>=1);

% For small networks, ensure the candidates are in large enough components.
if (nnodes<=10000)
[ncomps,comps]=find_conn_comps2(nnodes,G);
sel=zeros(1,length(vals));
for i=1:length(vals)
j=vals(i); k=1;
while ((sel(i)==0)&(comps{k}.n>=nsubnodes)&(k<=ncomps))
if (ismember(j,comps{k}.comps)==1)
sel(i)=1;
end
k=k+1;
end
end
vals=vals(find(sel==1));
end


ps=degs(vals); ps=ps/sum(ps);
seed=sample_discrete_rv(vals,ps);


% Perform random walks from the current network.
% Sample nodes according to the deviation from the desired density.
targetdensity=nsubedges/nsubnodes;
subinds=seed; flag=1;
while (flag==1)

[a,b]=find(G(subinds,:)>0); b=unique(b); neighbors=setdiff(b,subinds);
if (length(neighbors)==0)
flag=0;
else
ps=zeros(1,length(neighbors));
for i=1:length(neighbors)
tmpinds=[subinds neighbors(i)];
tmpG=G(tmpinds,tmpinds);
[a,b]=find(tmpG>0); subset=find(a<b);
val=length(subset)/(length(subinds)+1);
ps(i)=abs(val-targetdensity);
end
val=max(ps); ps=val-ps+0.001; ps=ps/sum(ps);
newind=sample_discrete_rv(neighbors,ps);
subinds=union(subinds,newind);
if (length(subinds)>=nsubnodes)
flag=0;
end
end


% Debug
[a,b]=find(G(subinds,subinds)>0); subset=find(a<b);
%fprintf('nsubnodes=%d, nsubedges=%d, density=%.4f, targetdensity=%.4f\n',length(subinds),length(subset),length(subset)/length(subinds),targetdensity);


end


% If the density is much smaller than expected, then resample the subgraph by starting with a seed with higher degrees.
density=length(subset)/length(subinds);
if (density<(0.8*targetdensity))

vals=find(degs>=10); 

% For small networks, ensure the candidates are in large enough components.
if (nnodes<=10000)
[ncomps,comps]=find_conn_comps2(nnodes,G);
sel=zeros(1,length(vals));
for i=1:length(vals)
j=vals(i); k=1;
while ((sel(i)==0)&(comps{k}.n>=nsubnodes)&(k<=ncomps))
if (ismember(j,comps{k}.comps)==1)
sel(i)=1;
end
k=k+1;
end
end
vals=vals(find(sel==1));
end


ps=degs(vals); ps=ps/sum(ps);
seed=sample_discrete_rv(vals,ps);


% Perform random walks from the current network.
% Sample nodes according to the deviation from the desired density.
targetdensity=nsubedges/nsubnodes;
subinds=seed; flag=1;
while (flag==1)

[a,b]=find(G(subinds,:)>0); b=unique(b); neighbors=setdiff(b,subinds);
if (length(neighbors)==0)
flag=0;
else
ps=zeros(1,length(neighbors));
for i=1:length(neighbors)
tmpinds=[subinds neighbors(i)];
tmpG=G(tmpinds,tmpinds);
[a,b]=find(tmpG>0); subset=find(a<b);
val=length(subset)/(length(subinds)+1);
ps(i)=abs(val-targetdensity);
end
val=max(ps); ps=val-ps+0.001; ps=ps/sum(ps);
newind=sample_discrete_rv(neighbors,ps);
subinds=union(subinds,newind);
if (length(subinds)>=nsubnodes)
flag=0;
end
end


% Debug
[a,b]=find(G(subinds,subinds)>0); subset=find(a<b);
%fprintf('nsubnodes=%d, nsubedges=%d, density=%.4f, targetdensity=%.4f\n',length(subinds),length(subset),length(subset)/length(subinds),targetdensity);

end

end


% Extract the subgraph spanned by subinds.
subG=G(subinds,subinds);


% If the number of edges is too large, then continue sampling the network.
[a,b]=find(subG>0); subset=find(a<b); k=length(subset);
if (k<=(nsubedges*2))
found=1;
end

end




