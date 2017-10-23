% Test the homology calculation codes on small networks.
% Swap networks and compute the homology again.

nrandoms=100;
filenames{1}='/home/chyeang/social_network/topodata/toynet.txt';
filenames{2}='/home/chyeang/social_network/topodata/karate.txt';
filenames{3}='/home/chyeang/social_network/data/Newman/dolphins.txt';
filenames{4}='/home/chyeang/social_network/data/Newman/lesmis.txt';
filenames{5}='/home/chyeang/social_network/data/Newman/football.txt';
filenames{6}='/home/chyeang/social_network/data/Newman/polbooks.txt';
filenames{7}='/home/chyeang/social_network/data/Newman/celegansneural.txt';
filenames{8}='/home/chyeang/social_network/topodata/netscience.txt';
filenames{9}='/home/chyeang/social_network/topodata/Bridge_Brook_raw_emat.edgelist2';

for ind=1:length(filenames)

fprintf('ind=%d, filename=%s\n',ind,filenames{ind});

% Import the network edge list.
filename=filenames{ind};
fp=fopen(filename);
fgets(fp); fgets(fp); fgets(fp); fgets(fp); flag=1;
nedges=0; edges=[];
while (flag==1)
s=fgets(fp);
if (s==-1)
flag=0;
else
nedges=nedges+1;
item=getitemval4(s,0,9);
edges(nedges,1)=atof2(item)+1;
item=getitemval4(s,1,9);
edges(nedges,2)=atof2(item)+1;
end
end
fclose(fp);

% Construct the graph.
nnodes=max(max(edges));
G=zeros(nnodes,nnodes);
for n=1:nnodes
G(n,n)=1;
end
for n=1:nedges
i1=edges(n,1); i2=edges(n,2);
G(i1,i2)=1; G(i2,i1)=1;
end

% Apply sequential reduction to the graph.
%newG=sequential_reduction(G);

% Debug
newG=G;

% For netscience, build the network of the largest component.
if (ind==8)
[ncomps,comps]=find_conn_comps(nnodes,G);
subG=G(comps{1}.comps,comps{1}.comps);
newG=subG;
nnodes=length(newG(1,:)); 
nedges=(sum(sum(full(newG)>0))-nnodes)/2;
end

% Construct clique complexes.
[nsimplices,simplices]=construct_clique_complex(newG);

% Evaluate the Betti numbers of the simplicial complex.
bs=evaluate_complex_homology(nsimplices,simplices);


% Randomly swap edges in the network.
templateG=newG; randbs=zeros(nrandoms,length(bs)); 
for randind=1:nrandoms
swapG=templateG; cnt=0;
for n=1:(nedges*10)
flag=0;
[a,b]=find(swapG>0);
tmp=find(a<b);
swapedges=transpose([a(tmp)';b(tmp)']);
while (flag==0)
vec=randperm(length(swapedges(:,1))); pair=vec(1:2);
if (length(intersect(swapedges(pair(1),:),swapedges(pair(2),:)))==0)
v1=swapedges(pair(1),1); v2=swapedges(pair(1),2);
v3=swapedges(pair(2),1); v4=swapedges(pair(2),2);
if ((swapG(v1,v4)==0)&(swapG(v2,v3)==0))
oldswapG=swapG;
swapG(v1,v2)=0; swapG(v2,v1)=0; 
swapG(v3,v4)=0; swapG(v4,v3)=0; 
swapG(v1,v4)=1; swapG(v4,v1)=1;
swapG(v2,v3)=1; swapG(v3,v2)=1;
flag=1; cnt=cnt+1;
end
end
end
end
templateG=swapG;

% Construct clique complexes.
[nrandsimplices,randsimplices]=construct_clique_complex(swapG);

% Evaluate the Betti numbers of the simplicial complex.
tmpbs=evaluate_complex_homology(nrandsimplices,randsimplices);
randbs(randind,1:length(tmpbs))=tmpbs;

% Debug
%fprintf('%d, ',randind);
%for i=1:length(bs)
%fprintf('%d ',randbs(randind,i));
%end
%fprintf('\n');

end

maxdim=length(simplices{nsimplices});

fprintf('filename=%s\n',filename);
fprintf('nnodes=%d, nedges=%d, maxcliquesize=%d, bettis=( ',nnodes,nedges,maxdim);
for i=1:length(bs)
fprintf('%d ',bs(i));
end
fprintf(')\n');
fprintf('rand bettis=( ');
for i=1:length(bs)
fprintf('%.2f(+/-%.2f) ',mean(randbs(:,i)),std(randbs(:,i)));
end
fprintf(')\n');
fprintf('\n');

end

