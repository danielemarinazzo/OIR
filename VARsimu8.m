%% simulation of three scalar processes - VAR(3)

function [Am,Su,Ak,Mv,lambdamax] = VARsimu8

Q=8; % number of series
Mv=[2 2 3 1]'; % number of series in each block

ca=[0.5 0.5];
cb=[0.4 0.4]; %cb=[0 0];
cc=[0.4 0.4 0.4];
cd=[0.7]; %cd=0;
ck=[-0.5 0.5 0.5 0.5];

pmax=2; % maximum lag
r1=0.9; f1=0.1; % oscillation channel 1
r3=0.9; f3=0.25; % oscillation channel 3

a1=ca(1); a2=ca(2);
b1=cb(1); b2=cb(2);
c1=cc(1); c2=cc(2); c3=cc(3);
d1=cd(1);
k1=ck(1); k2=ck(2); k3=ck(3); k4=ck(4);
Su=eye(Q); %Residual covariance ma trix (DIAGONAL)
B0=zeros(Q,Q); %Matrix of instantaneous effects:no instantaneous effects
Bk=NaN*ones(Q,Q,pmax);  
%effects at lag 1
Bk(1,:,1)=[2*r1*cos(2*pi*f1) 0 0 0 0 0 0 0];
Bk(2,:,1)=[k1 0 0 0 0 0 0 0];
Bk(3,:,1)=[0 0 2*r3*cos(2*pi*f3) 0 0 0 0 d1];
Bk(4,:,1)=[0 a2 0 0 0 0 0 0];
Bk(5,:,1)=[0 0 c1 0 0 0 0 0];
Bk(6,:,1)=[0 0 0 0 k2 0 0 0];
Bk(7,:,1)=[0 0 c3 0 0 k4 0 0];
Bk(8,:,1)=[0 0 0 b2 0 0 0 0];
%effects at lag 2
Bk(1,:,2)=[-r1^2 0 0 0 0 0 0 0];
Bk(2,:,2)=[0 0 0 0 0 0 0 0];
Bk(3,:,2)=[0 a1 -r3^2 0 0 0 0 0];
Bk(4,:,2)=[0 0 0 0 0 0 0 0];
Bk(5,:,2)=[0 0 0 0 0 k3 0 0];
Bk(6,:,2)=[0 0 c2 0 0 0 0 0];
Bk(7,:,2)=[0 0 0 0 0 0 0 0];
Bk(8,:,2)=[0 0 b1 0 0 0 0 0];
% concateno in matrice Bm
Bm=[];
for k=1:pmax
    Bm=[Bm Bk(:,:,k)];
end
Am=Bm; Ak=Bk; % no inst eff
E=eye(Q*pmax);AA=[Am;E(1:end-Q,:)];lambda=eig(AA);
lambdamax=max(abs(lambda));

end