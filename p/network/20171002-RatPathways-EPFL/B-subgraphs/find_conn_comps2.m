% Find all connected components of a graph.
% Difference from find_conn_comps.m: do not incur recursive functions.

function [nconncomps, conncomps] = find_conn_comps2(nnodes, G)

% Label nodes in the graph until all nodes are labeled.
labeled=zeros(1,nnodes); flag=1;

% Label all the singletons.
tmpG=G; tmpG=tmpG-diag(diag(tmpG));
%for i=1:nnodes
%tmpG(i,i)=0;
%end
degs=sum(tmpG>0);
subset=find(degs==0); k=0;
for i=1:length(subset)
k=k+1; labeled(subset(i))=k;
end


while (flag==1)

% Check if all nodes are labeled.
n=sum(labeled==0);

% Debug
%if (mod(n,1000)==0)
%n
%end

% If yes then return.
if (n==0)
flag=0;

% Otherwise pick up an unlabeled node and make its neighbors all have identical labels.
% Pick up the remaining node with the highest degree.
else
subset=find(labeled==0); 
ind=subset(find(degs(subset)>=max(degs(subset))));
ind=ind(1);

% Ensure all direct and indirect neighbors have identical labels.
subset=union(find(tmpG(ind,:)>0),find(tmpG(:,ind)>0));
neighborlabels=unique(labeled(subset));
if (neighborlabels==0)
k=k+1; labeled(ind)=k; labeled(subset)=k;
else
curlabel=setdiff(neighborlabels,0); curlabel=curlabel(1);
labeled(ind)=curlabel; labeled(subset)=curlabel;
for i=1:length(neighborlabels)
j=neighborlabels(i);
if ((j~=0)&(j~=curlabel))
labeled(find(labeled==j))=curlabel;
end
end
end

end

end


% Extract unique labels.  Each unique label constitutes one component.
uniquelabels=unique(labeled);
nconncomps=length(uniquelabels);
for i=1:nconncomps
j=uniquelabels(i);
subset=find(labeled==j);
conncomps{i}.n=length(subset);
conncomps{i}.comps=subset;
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


