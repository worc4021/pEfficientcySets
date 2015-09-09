clear all
close all
clc




A = [1.8,-0.5;
    -1.5,1.2];
B = [2;1];
D = eye(2);

G = [-1.75,1;
    1.2,2];
h = [1.5;1.3];

p = .2;

N = 300;

wMin = -2.68;
wMax = -wMin;
% P{Ax+Bw<=e}>=p ==> P{oAx+oBw<=oe}>=p ==> P{-X(oAx-oe)>=w}>=p ==>
% -X(oAx-oe) in Zp
% w in [-wMin,wMax] ==> Zp = { (w1+wMin)(w2+wMin) / (wMax-wMin)^2 >=p }
W = [diag([1/wMax,1/wMax]);diag([1/wMin,1/wMin])];
w1 = linspace(.9*wMin,wMax,N);
w2 = (wMax-wMin)^2*p./(w1-wMin)+wMin;
E = [w1',w2'];

%     fprintf('0 1 0\n');
%     fprintf('0 1 0\n');
% for i = 1:N
%     fprintf('1');
%     for j = 1:2
%         [n,d] = rat(E(i,j));
%         fprintf(' %d/%d', n,d);
%     end
%     fprintf('\n');
% end

% [ZpA,ZpB] = facetEnumeration([eye(2);E],[0;0;ones(N,1)]);
% 
% [V,v] = inequalityReduction([W;ZpA],[ones(4,1);ZpB]);
% 
Zp = Polyhedron('V',E,'R',eye(2));
V = Zp.A;
v = Zp.b;

[V,v] = inequalityReduction([W;V],[ones(4,1);v]);

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



w = sdpvar(2,1);
Cons = W*w<=ones(4,1);

LambdaNext = [];
lambdaNext = [];

iter = 1;
while and(~isContained(Lambda,lambda,LambdaNext,lambdaNext),iter<iterMax)
    
    if iter~=1
        Lambda = LambdaNext;
        lambda = lambdaNext;
    end
    
    LambdaNext=zeros(size(Lambda));
    lambdaNext=zeros(size(lambda));
    
        for i = 1:length(lambda)
            obj = -Lambda(i,:)*D*w;
            diagnostic = optimize(Cons,obj,opt);
            LambdaNext(i,:) = Lambda(i,:)*(A+B*K);
            lambdaNext(i) = lambda(i)-Lambda(i,:)*D*value(w);
        end
    Dummy = Polyhedron([Lambda;LambdaNext],[lambda;lambdaNext]);
    Dummy = Dummy.minHRep;
    [LambdaNext,lambdaNext] = inequalityReduction([Lambda;LambdaNext],[lambda;lambdaNext]);
    iter = iter+1;
end

Lambda = LambdaNext;
lambda = lambdaNext;
