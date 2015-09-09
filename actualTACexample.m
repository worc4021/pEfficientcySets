clear all
close all

n = 100;

A = [1.8,-0.5;
    -1.5,1.2];
B = [2;1];
D = eye(2);

G = [-1.75,1;
    1.2,2];
h = [1.5;1.3];

sigma = 1/25;

wMin = -50*sigma;
wMax = -wMin;
W = [diag([1/wMax,1/wMax]);diag([1/wMin,1/wMin])];

p = .8;

if 1
    
    SigmaI = eye(2)/sigma;
    
    offset = (1+cdf(SigmaI,[1,1]*wMax,[1,1]*wMin,0,[],0))/((wMax-wMin)^2);
    
    
    ll = wMin;
    ul = wMax;
    
    x1 = linspace(ll,ul,n);
    x2 = zeros(1,n);

    x0 = .25;

    opt = optimoptions('fsolve','TolFun',1e-16,'TolX',1e-16,'Display',...
        'off','Jacobian','on','Algorithm','Trust-Region-reflective');

    for i = n:-1:1
        x2(i) = fsolve(@(x)cdf(SigmaI,[x1(i),x],[wMin,wMin],p,1,offset),x0,opt);
        x0 = x2(i);
    end

    X2 = linspace(ll,wMax,n);
    X1 = zeros(1,n);

    x0 = .24;

    for i = n:-1:1
        X1(i) = fsolve(@(x)cdf(SigmaI,[x,X2(i)],[wMin,wMin],p,2,offset),x0,opt);
        x0 = X1(i);
    end


    range = 1:n;
    plot(x1(range),x2(range),'x',X1(range),X2(range),'x');

    
    % P{Ax+Bw<=e}>=p ==> P{oAx+oBw<=oe}>=p ==> P{-X(oAx-oe)>=w}>=p ==>
    % -X(oAx-oe) in Zp
    % w in [-wMin,wMax] ==> Zp = { (w1+wMin)(w2+wMin) / (wMax-wMin)^2 >=p }

    P = Polyhedron('V',[x1(range),X1(range);x2(range),X2(range)]','R',eye(2));
    V = P.A;
    v = P.b;

    [V,v] = inequalityReduction([V;W],[v;ones(4,1)]);
    save('ZpSet.mat','V','v');
else
    load('ZpSet.mat');
end

if 1
    X = sdpvar(2);
    opt = sdpsettings('solver','mosek','verbose',0);

    Cons = X(:)>=0;
    Obj = trace((X*G-eye(2))*(X*G-eye(2))');

    diagnostic = optimize(Cons,Obj,opt);

    M = value(X);

    [Lambda,lambda] = inequalityReduction(-V*M*A,v-V*M*h);

    [~,K,~,gamma2,~] = terminalController(A,B,D,eye(2),1,300);


    Xk = Polyhedron(Lambda,lambda);
    Sk = Polyhedron();
    Ws = Polyhedron(W,ones(4,1));

    iter = 1;
    iterMax = 50;

    while and(~(Xk<=Sk),iter<iterMax)
        if iter~=1
            Xk = Sk.minHRep;
        end
        Sk = (A+B*K)*Xk - D*Ws;
        iter = iter + 1;
    end
end

