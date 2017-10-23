% Find all connected components of a graph.

function [nconncomps, conncomps] = find_conn_comps(nnodes, G)

% Recursively label nodes in the graph until all nodes are labeled.
labeled=zeros(1,nnodes); flag=1;

% Label all the singletons.
%k=0;
%for i=1:nnodes
%tmp=find(G(i,:)>0);
%if (length(tmp)==0)
%k=k+1; labeled(i)=k;
%end
%end


% Label all the singletons.
tmpG=G; tmpG=tmpG-diag(diag(tmpG));
%for i=1:nnodes
%tmpG(i,i)=0;
%end
vec=sum(tmpG>0);
subset=find(vec==0); k=0;
for i=1:length(subset)
k=k+1; labeled(subset(i))=k;
end


while (flag==1)

% Check if all nodes are labeled.
n=0;
for i=1:nnodes
if (labeled(i)==0)
n=n+1;
end
end

% If yes, then return.
if (n==0)
flag=0;

% Otherwise pick up an unlabeled node and recursively label its neighbors.
else
i=1; j=0;
while ((i<=nnodes)&(j==0))
if (labeled(i)==0)
j=i;
end
i=i+1;
end
k=k+1;
labeled=recurse_label2(j,nnodes,G,labeled,k);
end

end


% Find the connected components from the labels.
nconncomps=k; clear conncomps;
for i=1:nconncomps
conncomps{i}.n=0;
conncomps{i}.comps=[];
end

for i=1:nnodes
j=labeled(i);
conncomps{j}.n=conncomps{j}.n+1;
conncomps{j}.comps(conncomps{j}.n)=i;
end

% Sort the connected components by size.
tmp=zeros(1,nconncomps);
for i=1:nconncomps
tmp(i)=conncomps{i}.n;
end
[Y,I]=sort(tmp,'descend'); clear tmp2;
for i=1:nconncomps
tmp2{i}=conncomps{I(i)}; 
end
conncomps=tmp2;

% Discard the connected components which do not have type 1 or type 2 nodes.
discard=zeros(1,nconncomps);
for i=1:nconncomps
if (conncomps{i}.n<=0)
discard(i)=1;
end
end
tmp=find(discard==0);
for i=1:length(tmp)
tmp2{i}=conncomps{tmp(i)};
end
nconncomps=length(tmp); conncomps=tmp2;

clear tmp tmp2 Y I discard;


