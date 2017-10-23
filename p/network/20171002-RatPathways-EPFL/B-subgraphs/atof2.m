% Convert a string into a real number.
% Also works for exponential representation.

function val = atof2(s)

% If s=NA then return NaN.
if (strcmp(s,'NA')==1)
val=NaN;

else
val=0; sgn=1; cnt=1; afterpoint=0; afterexp=0; 
expval=0; expsign=1; power=1.0;
while (cnt<=length(s))
ch=s(cnt);
if ((cnt==1)&(ch=='-'))
   sgn=-1;
elseif ((ch=='-')&(afterexp==1))
   expsign=-1;
elseif (ch=='.')
   afterpoint=1; power=0.1;
elseif ((ch=='e')|(ch=='E'))
   afterexp=1; afterpoint=0;
elseif ((ch>='0')&(ch<='9')&(afterpoint==0))
   if (afterexp==0)
     val=val*10+(ch-'0');
   else
     expval=expval*10+(ch-'0');
     end
elseif ((ch>='0')&(ch<='9')&(afterpoint==1))
   if (afterexp==0)
     val=val+(ch-'0')*power;
   else
      expval=expval+(ch-'0')*power;
      end
   power=power*0.1;
   end
cnt=cnt+1;
end

if (sgn==-1)
  val=val*(-1);
  end
if (expsign==-1)
  expval=expval*(-1);
  end

val=val*10^expval;

end
