% TACexample

vV = vertexEnumeration(Lambda,lambda);

la = sdpvar(size(vV,1),1);
Cons = [ones(1,length(la))*la == 1;la(:)>=0;la(:)<=1];
obj = (rand(1,length(la))*la)^2;
optimize(Cons,obj,opt);

n = 250;
x = [(value(la)'*vV)',zeros(2,n)];

for i = 1:n
    wR = rand(2,1)*(wMax-wMin)+wMin;
    x(:,i+1) = (A+B*K)*x(:,i) + D*wR;
end

plot(Polyhedron(Lambda,lambda),'alpha',.2);
hold('on');
plot(x(1,:),x(2,:),'x');
hold('off')