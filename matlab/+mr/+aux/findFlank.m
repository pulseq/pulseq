function xf=findFlank(x,f,c)
%findFlank: find the x coordinate of the left flank of function f
%   Finds the furst value within x for which abs(f) is greater than 
%   c*max(abs(f)). If xf is not the first element of x, then a linear 
%   interpolation is applied.
m=max(abs(f));
f=abs(f)/m-c;
i=find(f>0,1);
if i>1
    f0=f(i-1);
    f1=f(i);
    xf=(f1*x(i-1)-f0*x(i))/(f1-f0); 
else
    xf=x(1);
end

end

