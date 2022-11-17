%% OIR - computation of delta of the O-information rate when the block Xj is added to the group X-j

% ij -  complete vector of indexes
% j - index of the target to analyze within ij

function out=oir_deltaO(A,C,K,Su,Mv,ij,j,Fs,nfft,do_td)

narginchk(7,10)
if nargin<10, do_td='n'; end
if nargin<9, nfft=512; end
if nargin<8, Fs=1; end %default non compute time domain measures through submodels

% if ~exist(nargin)

% ij=[1 2 3]; % vector index of the blocks X={X_1...X_N} to consider for OIR computation
% j=1; % index of the target block Xj

%%% computation of dO(X-j;Xj)
assert(ismember(j,ij)) %verify target belongs to group
ii=setdiff(ij,j); % driver indexes
N=length(ij); % order of the OIR to compute

out1=oir_mir(A,C,K,Su,Mv,ii,j,Fs,nfft,do_td); %MIR between Xj and X-j
i_cs=nchoosek(ii,N-2); % n. of combinations to compute the sum of the deltaO 
dO12=0; dO1_2=0; dO2_1=0; dO1o2=0;
dO12f=zeros(nfft,1); dO1_2f=dO12f; dO2_1f=dO12f; dO1o2f=dO12f;
for cnt=1:N-1
    outtmp=oir_mir(A,C,K,Su,Mv,i_cs(cnt,:),j,Fs,nfft,do_td); %MIR between Xj and X-ij
    % time domain measures
    dO12=dO12+outtmp.I12;
    dO1_2=dO1_2+outtmp.T1_2;
    dO2_1=dO2_1+outtmp.T2_1;
    dO1o2=dO1o2+outtmp.I1o2;
    % spectral functions
    dO12f=dO12f+outtmp.f12;
    dO1_2f=dO1_2f+outtmp.f1_2;
    dO2_1f=dO2_1f+outtmp.f2_1;
    dO1o2f=dO1o2f+outtmp.f1o2; 
end
% time domain measures: last term 
dO12=dO12+(2-N)*out1.I12; % dO(X-j;Xj) - Eq.(5)
dO1_2=dO1_2+(2-N)*out1.T1_2; % Eq.(8a)
dO2_1=dO2_1+(2-N)*out1.T2_1; % Eq.(8b)
dO1o2=dO1o2+(2-N)*out1.I1o2; % Eq.(8c)
% spectral functions: last term
dO12f=dO12f+(2-N)*out1.f12; % Eq.(23)
dO1_2f=dO1_2f+(2-N)*out1.f1_2; % Eq.(24)
dO2_1f=dO2_1f+(2-N)*out1.f2_1; % Eq.(24)
dO1o2f=dO1o2f+(2-N)*out1.f1o2; % Eq.(24)

% output structure
out.dO1_2=dO1_2;
out.dO2_1=dO2_1;
out.dO1o2=dO1o2;
out.dO12=dO12;
out.dO1_2f=dO1_2f;
out.dO2_1f=dO2_1f;
out.dO1o2f=dO1o2f;
out.dO12f=dO12f;
out.freq=out1.freq;

end