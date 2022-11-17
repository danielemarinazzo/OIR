%% simulation of three scalar processes - VAR(3)

function [Am,Su,Ak,Mv,lambdamax] = VARsimu10
% clear; close all; clc

Q=10; % number of series
Mv=[4 1 2 1 2]'; % number of series in each block

pmax=2; % maximum lag
r1=0.9; f1=0.1; % oscillation channel 1
r2=0.8; f2=0.25; % oscillation channel 3

coup=[1 2 1 0.5; 1 4 2 -0.5; 2 3 1 0.5; 3 4 1 0.2;...
      6 7 1 0.3; 7 6 2 0.3;...
      9 10 1 0.4; 10 9 2 -0.2;
      3 8 1 0.3; 2 8 2 0.4;...      
      5 8 1 -0.4;...
      7 8 1 0.3;...
      8 9 1 0.7;...
      10 4 1 0.5;...
      ]; % in each row: "i j k c" to impose coupling from i to j at lag k with coeff c 

Su=eye(Q); %Residual covariance ma trix (DIAGONAL)
B0=zeros(Q,Q); %Matrix of instantaneous effects:no instantaneous effects
Bk=NaN*ones(Q,Q,pmax);  
%effects at lag 1
Bk(1,:,1)=[2*r1*cos(2*pi*f1) 0 0 0 0 0 0 0 0 0];
Bk(2,:,1)=[0 0 0 0 0 0 0 0 0 0];
Bk(3,:,1)=[0 0 0 0 0 0 0 0 0 0];
Bk(4,:,1)=[0 0 0 0 0 0 0 0 0 0];
Bk(5,:,1)=[0 0 0 0 2*r1*cos(2*pi*f1) 0 0 0 0 0];
Bk(6,:,1)=[0 0 0 0 0 0 0 0 0 0];
Bk(7,:,1)=[0 0 0 0 0 0 2*r1*cos(2*pi*f1) 0 0 0];
Bk(8,:,1)=[0 0 0 0 0 0 0 2*r2*cos(2*pi*f2) 0 0];
Bk(9,:,1)=[0 0 0 0 0 0 0 0 0 0];
Bk(10,:,1)=[0 0 0 0 0 0 0 0 0 0];
%effects at lag 2
Bk(1,:,2)=[-r1^2 0 0 0 0 0 0 0 0 0];
Bk(2,:,2)=[0 0 0 0 0 0 0 0 0 0];
Bk(3,:,2)=[0 0 0 0 0 0 0 0 0 0];
Bk(4,:,2)=[0 0 0 0 0 0 0 0 0 0];
Bk(5,:,2)=[0 0 0 0 -r1^2 0 0 0 0 0];
Bk(6,:,2)=[0 0 0 0 0 0 0 0 0 0];
Bk(7,:,2)=[0 0 0 0 0 0 -r1^2 0 0 0];
Bk(8,:,2)=[0 0 0 0 0 0 0 -r2^2 0 0];
Bk(9,:,2)=[0 0 0 0 0 0 0 0 0 0];
Bk(10,:,2)=[0 0 0 0 0 0 0 0 0 0];

%set time-lagged coupling coefficients
for ic=1:size(coup,1)
    Bk(coup(ic,2),coup(ic,1),coup(ic,3))=coup(ic,4);
end


% concateno in matrice Bm
Bm=[];
for k=1:pmax
    Bm=[Bm Bk(:,:,k)];
end
Am=Bm; Ak=Bk; % no inst eff
E=eye(Q*pmax);AA=[Am;E(1:end-Q,:)];lambda=eig(AA);
lambdamax=max(abs(lambda));

end