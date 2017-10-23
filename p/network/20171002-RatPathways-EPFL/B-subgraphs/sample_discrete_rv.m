% Sample a value from a discrete distribution.

% RA: discretize(rand(1, n) * sum(p), [0; cumsum(p)])
% https://www.mathworks.com/matlabcentral/fileexchange/21912-sampling-from-a-discrete-distribution

function randval = sample_discrete_rv(vals, ps)

nstates=length(vals);
[Y,I]=sort(ps,'descend');
bds=zeros(1,nstates);
bds(1)=Y(1);
for n=2:nstates
bds(n)=bds(n-1)+Y(n);
end
bds=bds/sum(Y);

val=rand(1); lwd=1; upp=nstates; randval=0;
if (val<=bds(lwd))
randval=I(lwd);
elseif (val>=bds(upp))
randval=I(upp);
end
while (randval==0)
mid=ceil((lwd+upp)/2);
if ((val>bds(mid-1))&(val<=bds(mid)))
randval=I(mid);
elseif ((val>=bds(mid))&(val<bds(mid+1)))
randval=I(mid);
elseif (val>bds(mid))
lwd=mid;
elseif (val<bds(mid))
upp=mid;
end
end

randval=vals(randval);

