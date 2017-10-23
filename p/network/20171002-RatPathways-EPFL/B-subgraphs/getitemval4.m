% Get the selected item in a string.
% Difference from getitemval3: fix the bug that the last entry contains \n.

function item = getitemval4(s, ind, sepch)

tabcnt=0; curind=1;
while ((tabcnt<ind)&(curind<=length(s)))
ch=s(curind);
if (ch==sepch)
tabcnt=tabcnt+1;
elseif ((ch=='\n')&(tabcnt<ind))
tabcnt=ind+1;
end
curind=curind+1;
end

if (tabcnt==ind)
cnt=1; template='';
while ((curind<=length(s))&(s(curind)~=sepch))
%while ((curind<=length(s))&(s(curind)~=sepch)&(s(curind)~='\n'))
template(cnt)=s(curind); curind=curind+1; cnt=cnt+1;
end
if (cnt>1)
item=template;
else
item='';
end
else
item='';
end

if (length(item)>0)
if (isspace(item(length(item)))==1)
item=item(1:(length(item)-1));
end
end
