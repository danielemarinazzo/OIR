%% test script - computation of MIR and its decompositions for two (block) processes taken from a multivariate simulated process
clear; close all; clc;

%% Parameters
i_1=[1]; %index of first (vector) process Z1
i_2=[2]; %index of second (vector) process Z2

nfft=512;   % frequency resolution for spectral causalities
Fs=1;
do_td='y';% if 'y', compute time domain measures through submodels (to verify the equivalence with mean of frequency domain measures)

% SIMULATION
%%% simple bivariate process
[Am,Su,Ak,Mv] = VARsimu2;
%decommentare le seguenti se voglio testare la non-strict causality
Su(1,2)=0.5; Su(2,1)=0.5;

%%% 8-variate process of Faes Biol.Cyb 2013
% [Am,Su,Ak,Mv] = VARsimu8;
%decommentare le seguenti se voglio testare la non-strict causality
% Su(1,4)=0.7; Su(4,1)=0.7;
% Su(2,8)=0.4; Su(8,2)=0.4;
% Su(6,8)=0.4; Su(8,6)=0.4;
% Su(1,8)=0.4; Su(8,1)=0.4;
% decommentare questa se voglio assumere tutti i processi scalari (no blocchi)
% Mv=ones(sum(Mv),1);


% Note: Mv contains the dimension of each block
%% analysis
% function which calls others to perform the complete bivariate block computation
%%% equivalent ISS model
[A,C,K,~] = oir_ar2iss(Am); % Eq. (26)
ret=oir_mir(A,C,K,Su,Mv,i_1,i_2,Fs,nfft,do_td); 
I12=ret.I12; % Eq.(12a)
T1_2=ret.T1_2; T2_1=ret.T2_1; % Eq.(12b)
I1o2=ret.I1o2; % Eq.(12c)
f12=ret.f12; % Eq.(16)
f1_2=ret.f1_2; f2_1=ret.f2_1; % Eq.(18)
f1o2=ret.f1o2; % Eq.(19)


%% plots and disps
%%% display time-domain dO values (both as computed, and as average of the spectral measure)
disp('F1;2')
disp([mean(f12)/2 I12])
disp('F1->2')
disp([mean(f1_2)/2 T1_2])
disp(' ')
disp('F2->1')
disp([mean(f2_1)/2 T2_1])
disp('F1o2')
disp([mean(f1o2)/2 I1o2])

%%% plot spectral measures
ymin=1.05*min([min(f12) min(f1_2) min(f2_1) min(f1o2)]);
ymax=1.05*max([max(f12) max(f1_2) max(f2_1) max(f1o2)]);
nfreq=ret.freq*2*pi; %omega
lw=1.5;
figure(1)
subplot(1,3,1); plot(nfreq,f12,'k','linewidth',lw); legend('f_{Z_1;Z_2}');xlabel('\omega')
xlim([0 pi]); ylim([ymin ymax]); title('total coupling')
subplot(1,3,2); plot(nfreq,f1_2,'linewidth',lw);hold on; plot(nfreq,f2_1,'linewidth',lw);
xlim([0 pi]); ylim([ymin ymax]); title('causal coupling')
legend('f_{Z_1 \rightarrow Z_2}','f_{Z_2 \rightarrow Z_1}');xlabel('\omega')
subplot(1,3,3); plot(nfreq,f1o2,'color',[0.5 0.5 0.5],'linewidth',lw);hold on; legend('f_{Z_1 \cdot Z_2}');xlabel('\omega')
xlim([0 pi]); ylim([ymin ymax]); title('mixing')

