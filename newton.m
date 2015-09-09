function [y,i] = newton(fhandle,x0,Tol,maxIter)

x = x0;
y = inf;
i = 1;
while and(abs(x-y)>Tol,i<maxIter)
    if i ~= 1
        x = y;
    end
    
    [f,df] = fhandle(x);
    y = x-df*f;
    i = i+1;
end