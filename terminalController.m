function [Pt,Kt,KWt,gamma2,info] = terminalController(A,B,D,Q,R,gamma)

X = sdpvar(size(A,1));
Y = sdpvar(size(B,2),size(B,1));
P = sdpvar(size(A,1));
if isempty(gamma)
    gammaSQ = sdpvar(1);
else
    gammaSQ = gamma;
end

BLK11 = blkdiag(X,gammaSQ*eye(size(X)));
BLK12 = [X*A'+Y'*B',(sqrt(Q)*X)',(sqrt(R)*Y)';
         D',zeros(size(X)),zeros(size(Y'))];
BLK22 = blkdiag(X,eye(size(Q)),eye(size(R)));

SCHUR = [P,eye(size(P));eye(size(P)),X];

SDP = [X >= 0, [BLK11,BLK12;BLK12',BLK22] >= 0, SCHUR>=0];

opt = sdpsettings('verbose',0); %'solver','mosek-sdp',

info = optimize(SDP,trace(P),opt);

gamma2 = value(gammaSQ);
Xt = value(X);
Pt = Xt\eye(size(Xt));
Kt = value(Y)*Pt;
KWt = (gamma2*eye(size(Pt))-D'*Pt*D)\D'*Pt*(A+B*Kt);