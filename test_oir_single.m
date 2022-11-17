%% test script for the O-information rate
% computes only one delta OIR: deltaO obtained adding the block Xj to the group X-j
clear; close all; clc;

%% Parameters
nfft=512;   % frequency resolution for spectral causalities
Fs=1;
do_td='y';% if 'y', compute time domain measures through submodels (to verify the equivalence with mean of frequency domain measures)

% SIMULATION
[Am,Su,Ak,Mv] = VARsimu8;
%decommentare le seguenti se voglio testare la non-strict causality
% Su(1,4)=0.7; Su(4,1)=0.7;
% Su(2,8)=0.4; Su(8,2)=0.4;
% Su(6,8)=0.4; Su(8,6)=0.4;
% Su(1,8)=0.4; Su(8,1)=0.4;

% Note: Mv contains the dimension of each block

%% analysis
sources=[1 2 4]; % indexes of the N-1 blocks forming X-j
target=3; % Index of the block Xj
multiplet=[sources target]; %assigned group of block processes X^N

[A,C,K,~] = oir_ar2iss(Am); % Eq. (26) - equivalent ISS model
out=oir_deltaO(A,C,K,Su,Mv,multiplet,target,Fs,nfft,do_td);
dO12=out.dO12; % Eq.(5)
dO1_2=out.dO1_2; % Eq.(8a)
dO2_1=out.dO2_1; % Eq.(8b)
dO1o2=out.dO1o2; % Eq.(8c)
dO12f=out.dO12f; % Eq.(23)
dO1_2f=out.dO1_2f; % Eq.(24)
dO2_1f=out.dO2_1f; % Eq.(24)
dO1o2f=out.dO1o2f; % Eq.(24)

%% plots and disps
lw=1.2;
nfreq=out.freq*2*pi; %omega
plot(nfreq',dO12f,'linewidth',lw);
hold on;
plot(nfreq',dO1_2f,'linewidth',lw);
plot(nfreq',dO2_1f,'linewidth',lw);
plot(nfreq',dO1o2f,'linewidth',lw);
legend('\delta_{X^N_{-j} ; X_j}','\delta_{X^N_{-j} \rightarrow X_j}','\delta_{X_j \rightarrow X^N_{-j}}' ,'\delta_{X_j \cdot X^N_{-j}}')
xlabel('\omega')
xlim([0 pi]);

%%% display time-domain dO values (both as computed and averaged from the spectral measures)
disp('delta OIR Xj ; X-j')
disp([dO12 mean(dO12f)/2])
disp('delta OIR X-j -> Xj')
disp([dO1_2 mean(dO1_2f)/2])
disp('delta OIR Xj -> X-j')
disp([dO2_1 mean(dO2_1f)/2])
disp('delta OIR Xj o X-j')
disp([dO1o2 mean(dO1o2f)/2])
