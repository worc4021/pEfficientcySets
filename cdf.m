function [po,dp] = cdf(SigmaI,z,la,p,sel,offset)
    AbsTol = 1e-32;
    RelTol = 1e-64;
    Sc = sqrt(det(SigmaI))/(2*pi);
    fun = @(x,y)Sc*exp(-.5*(x.*SigmaI(1,1).*x+2*x.*SigmaI(1,2).*y+y.*SigmaI(2,2).*y))+offset;
    po = p-integral2(fun,la(1),z(1),la(2),z(2),'AbsTol',AbsTol,'RelTol',RelTol);
    dp = [integral(@(y)fun(z(1),y),la(2),z(2),'AbsTol',AbsTol,'RelTol',RelTol);
          integral(@(x)fun(x,z(2)),la(1),z(1),'AbsTol',AbsTol,'RelTol',RelTol);]*-1;
    dp = dp(sel);
end