%% OIR - frequency domain bivariate Granger Causality from reduced ISS model

% A,CR,KR,VR - parameters of the ISS reduced model to analyze
% j1,j2 - vector indexes of the two blocks for which to compute the bivariate measures
% j1 and j2 must cover the whole dimension of the reduced model: length(j1)+length(j2)=Mr, [j1 j2]=1:Mr

function [f12,f1_2,f2_1,f1o2,S,freq] = oir_fdGC(A,CR,KR,VR,j1,j2,Fs,nfft)

narginchk(6,8)
if nargin<8, nfft=512; end
if nargin<7, Fs=1; end

Mr=length(j1)+length(j2);
assert(Mr==size(CR,1));
assert(sum([j1 j2]==1:Mr)==Mr);

bVR{1,1}=VR(j1,j1); bVR{1,2}=VR(j1,j2);
bVR{2,1}=VR(j2,j1); bVR{2,2}=VR(j2,j2);

freq = (0:nfft-1)*(Fs/(2*nfft)); % frequency axis
s = exp(1j*2*pi*freq/Fs); % vector of complex exponentials
H = issMA(A,CR,KR,s); % transfer function from ISS parameters (SSGC toolbox Barnett)

S=zeros(Mr,Mr,nfft); % Spectral Matrix S(w)
bH=cell(2,2,nfft); 
bS=cell(2,2,nfft); 
f12=nan*ones(nfft,1);
f1_2=nan*ones(nfft,1); f2_1=f1_2; % Geweke spectral causality from sources to destinations
for n=1:nfft % at each frequency
    S(:,:,n)  = H(:,:,n)*VR*H(:,:,n)'; %PSD
    
    %%% extraction of blocks indexes
    bH{1,1,n}=H(j1,j1,n); bH{1,2,n}=H(j1,j2,n);
    bH{2,1,n}=H(j2,j1,n); bH{2,2,n}=H(j2,j2,n);  
    bS{1,1,n}=S(j1,j1,n); bS{1,2,n}=S(j1,j2,n);
    bS{2,1,n}=S(j2,j1,n); bS{2,2,n}=S(j2,j2,n);
    
    %instterm=bH{2,2,n}*bVR{2,1}*bH{2,1,n}' + bH{2,1,n}*bVR{1,2}*bH{2,2,n}';
    %f(n) = log( abs(det(bS{2,2,n})) / abs(det(bS{2,2,n}-bH{2,1,n}*bVR{1,1}*bH{2,1,n}'- instterm)) );
    f2_1(n) = log( abs(det(bS{1,1,n})) / abs(det(bH{1,1,n}*bVR{1,1}*bH{1,1,n}')));
    f1_2(n) = log( abs(det(bS{2,2,n})) / abs(det(bH{2,2,n}*bVR{2,2}*bH{2,2,n}')));
    
    f12(n)=log( abs(det(bS{1,1,n}))*abs(det(bS{2,2,n})) / abs(det(S([j1 j2],[j1 j2],n))) );
                   
end
f1o2=f12-f1_2-f2_1;

end


%% Taken from SSGC Toolbox of Lionel Barnett
% L. Barnett and A. K. Seth, Granger causality for state-space models, Phys.Rev. E 91(4) Rapid Communication, 2015.

function H = issMA(A,C,K,z)

% Compute moving-average representation (transfer function) for a state space
% model in innovations form (eq. 4 in the reference article).
%
% A,C,K - innovations form state space parameters
% z     - a vector of points on the unit circle in the complex plane
%
% H     - transfer function

[m,m1]  = size(A); assert(m1 == m);
[n,m1]  = size(C); assert(m1 == m);
[m1,n1] = size(K); assert(n1 == n && m1 == m);
z = z(:);

h = length(z);
In = eye(n);
Im = eye(m);
H = zeros(n,n,h);
for k = 1:h
    H(:,:,k) = In + C*((z(k)*Im-A)\K); % eq. 4
end
end

