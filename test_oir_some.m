%% test script for the O-information rate - complete and decomposed
% given a sequence of indexes, computes MIR for the first triplet, and OIR and deltaOIR for the multiplets with N>=3
clear; close all; clc;

%% Parameters
nfft=512;   % frequency resolution for spectral causalities
Fs=1; %sampling frequency
do_td='y';% if 'y', compute time domain measures through submodels (to verify the equivalence with mean of frequency domain measures)

% SIMULATION
[Am,Su,Ak,Mv] = VARsimu8;
%decommentare le seguenti se voglio testare la non-strict causality
% Su(1,4)=0.7; Su(4,1)=0.7;
% Su(2,8)=0.4; Su(8,2)=0.4;
% Su(6,8)=0.4; Su(8,6)=0.4;
% Su(1,8)=0.4; Su(8,1)=0.4;

Mv=[1 1 2 2 2]; % override Mv to have more blocks
% Mv=[3 3 2];

M=length(Mv);

%% analysis
% selection of a combination of indexes
rng('shuffle');% to generate different permutation at each run of the script
% iM=randperm(M); %random ordering of variable indexes - show one random path
iM=[5 1 4 3 2];


%%% equivalent ISS model
[A,C,K,~] = oir_ar2iss(Am); % Eq. (26)
%%% MIR for the triplet in iM(1:3)
out3_2=oir_mir(A,C,K,Su,Mv,iM(3),iM(2),Fs,nfft,do_td);
I32=out3_2.I12; T3_2=out3_2.T1_2; T2_3=out3_2.T2_1; I3o2=out3_2.I1o2;
f32=out3_2.f12; f3_2=out3_2.f1_2; f2_3=out3_2.f2_1; f3o2=out3_2.f1o2;
out3_1=oir_mir(A,C,K,Su,Mv,iM(3),iM(1),Fs,nfft,do_td);
I31=out3_1.I12; T3_1=out3_1.T1_2; T1_3=out3_1.T2_1; I3o1=out3_1.I1o2;
f31=out3_1.f12; f3_1=out3_1.f1_2; f1_3=out3_1.f2_1; f3o1=out3_1.f1o2;
out3_12=oir_mir(A,C,K,Su,Mv,iM(3),[iM(1) iM(2)],Fs,nfft,do_td);
I312=out3_12.I12; T3_12=out3_12.T1_2; T12_3=out3_12.T2_1; I3o12=out3_12.I1o2;
f312=out3_12.f12; f3_12=out3_12.f1_2; f12_3=out3_12.f2_1; f3o12=out3_12.f1o2;

% out3=oir_deltaO(Am,Su,Mv,iM(1:3),iM(3),nfft,Fs);
% [out3_1.F12+out3_2.F12-out3_12.F12 out3.dO12]

%%% delta OIR when adding iM(N) to iM(1:N)
dO=nan*ones(M,1); dO1_2=dO; dO2_1=dO; dO1o2=dO;
dOf=cell(M,1); dO1_2f=dOf; dO2_1f=dOf; dO1o2f=dOf;
for N=3:M % vary target inside the multiplet
    out=oir_deltaO(A,C,K,Su,Mv,iM(1:N),iM(N),Fs,nfft,do_td);
    dO(N-1)=out.dO12;
    dO1_2(N-1)=out.dO1_2; % Eq.(8a)
    dO2_1(N-1)=out.dO2_1; % Eq.(8b)
    dO1o2(N-1)=out.dO1o2; % Eq.(8c)
    dOf{N-1}=out.dO12f; % Eq.(23)
    dO1_2f{N-1}=out.dO1_2f; % Eq.(24)
    dO2_1f{N-1}=out.dO2_1f; % Eq.(24)
    dO1o2f{N-1}=out.dO1o2f; % Eq.(24)
end

%%% OIR using recursion
OIR=nan*ones(M,1); %OIR
OIRf=cell(M,1); %spectral OIR
OIR(3)=dO(2);
OIRf{3}=dOf{2};
for N=4:M
    OIR(N)=OIR(N-1)+dO(N-1);
    OIRf{N}=OIRf{N-1}+dOf{N-1};
end



%% plots
lw=1.2;
nfreq=(out.freq*2*pi)'; %omega


figure(1)
subplot(2,2,1);
xtesto=0.1; ymax=max(f312); ymin=min(f312); dy=(ymax-ymin)/12;
plot(nfreq,f312,'k','linewidth',lw);hold on;
text(xtesto,ymax,num2str(mean(f312)/2),'color','k')
plot(nfreq,f31,'b','linewidth',lw);
text(xtesto,ymax-dy,num2str(mean(f31)/2),'color','b')
plot(nfreq,f32,'g','linewidth',lw);hold on;
text(xtesto,ymax-2*dy,num2str(mean(f32)/2),'color','g')
plot(nfreq,f31+f32-f312,'r','linewidth',lw);
text(xtesto,ymax-3*dy,num2str(mean(f31+f32-f312)/2),'color','r')
plot(nfreq,dOf{2},'r--','linewidth',lw);hold on;
text(xtesto,ymax-4*dy,num2str(mean(dOf{2})/2),'color','r') %verification: should see only one red line
xlabel('\omega'); xlim([0 pi]);
legend(['f_{X_{' int2str(iM(3))  '} ; X_{' int2str(iM(1:2)) '}}:' num2str(I312)],...
        ['f_{X_{' int2str(iM(3))  '} ; X_' int2str(iM(1)) '}:' num2str(I31)],...
        ['f_{X_{' int2str(iM(3))  '} ; X_' int2str(iM(2)) '}:' num2str(I32)],...
        ['\delta_{X_{' int2str(iM(1:2))  '} ; X_' int2str(iM(3)) '}:' num2str(dO(2))]); legend('boxoff')
title(['OIR(3) decomposition for X_' int2str(iM(1)) ',X_' int2str(iM(2)) ',X_' int2str(iM(3)) ])

xtesto=0.1; ymax=max(f312); ymin=min(f312); dy=(ymax-ymin)/15;
subplot(2,2,2);
plot(nfreq,f312,'k','linewidth',lw);hold on;
text(xtesto,ymax,num2str(mean(f312)/2),'color','k')
plot(nfreq,f3_12,'m','linewidth',lw);
text(xtesto,ymax-dy,num2str(mean(f3_12)/2),'color','m')
plot(nfreq,f12_3,'c','linewidth',lw);
text(xtesto,ymax-2*dy,num2str(mean(f12_3)/2),'color','c')
plot(nfreq,f3o12,'y','linewidth',lw);
text(xtesto,ymax-3*dy,num2str(mean(f3o12)/2),'color','y')
xlabel('\omega'); xlim([0 pi]);
legend(['f_{X_{' int2str(iM(3))  '};X_{' int2str(iM(1:2)) '}}:' num2str(I312)],...
        ['f_{X_{' int2str(iM(3))  '} \rightarrow X_{' int2str(iM(1:2)) '}}:' num2str(T3_12)],...
        ['f_{X_{' int2str(iM(1:2))  '} \rightarrow X_' int2str(iM(3)) '}:' num2str(T12_3)],...
        ['f_{X_{' int2str(iM(3))  '} \cdot X_{' int2str(iM(1:2)) '}}:' num2str(I3o12)]); legend('boxoff')
title(['causal decomposition of MIR between X_' int2str(iM(3)) ' and X_{' int2str(iM(1:2)) '}'])
    
xtesto=0.1; ymax=max(f31); ymin=min(f31); dy=(ymax-ymin)/15;
subplot(2,2,3);
plot(nfreq,f31,'b','linewidth',lw);hold on;
text(xtesto,ymax,num2str(mean(f31)/2),'color','b')
plot(nfreq,f3_1,'m','linewidth',lw);
text(xtesto,ymax-dy,num2str(mean(f3_1)/2),'color','m')
plot(nfreq,f1_3,'c','linewidth',lw);
text(xtesto,ymax-2*dy,num2str(mean(f1_3)/2),'color','c')
plot(nfreq,f3o1,'y','linewidth',lw);
text(xtesto,ymax-3*dy,num2str(mean(f3o1)/2),'color','y')
xlabel('\omega'); xlim([0 pi]);
legend(['f_{X_{' int2str(iM(3))  '};X_' int2str(iM(1)) '}:' num2str(I31)],...
        ['f_{X_{' int2str(iM(3))  '} \rightarrow X_' int2str(iM(1)) '}:' num2str(T3_1)],...
        ['f_{X_{' int2str(iM(1))  '} \rightarrow X_' int2str(iM(3)) '}:' num2str(T1_3)],...
        ['f_{X_{' int2str(iM(3))  '} \cdot X_' int2str(iM(1)) '}:' num2str(I3o1)]); legend('boxoff')
title(['causal decomposition of MIR between X_' int2str(iM(3)) ' and X_{' int2str(iM(1)) '}'])

xtesto=0.1; ymax=max(f32); ymin=min(f32); dy=(ymax-ymin)/15;
subplot(2,2,4);
plot(nfreq,f32,'g','linewidth',lw);hold on;
text(xtesto,ymax,num2str(mean(f32)/2),'color','g')
plot(nfreq,f3_2,'m','linewidth',lw);
text(xtesto,ymax-dy,num2str(mean(f3_2)/2),'color','m')
plot(nfreq,f2_3,'c','linewidth',lw);
text(xtesto,ymax-2*dy,num2str(mean(f2_3)/2),'color','c')
plot(nfreq,f3o2,'y','linewidth',lw);
text(xtesto,ymax-3*dy,num2str(mean(f3o2)/2),'color','y')
xlabel('\omega'); xlim([0 pi]);
legend(['f_{X_{' int2str(iM(3))  '};X_' int2str(iM(2)) '}:' num2str(I32)],...
        ['f_{X_{' int2str(iM(3))  '} \rightarrow X_' int2str(iM(2)) '}:' num2str(T3_2)],...
        ['f_{X_{' int2str(iM(2))  '} \rightarrow X_' int2str(iM(3)) '}:' num2str(T2_3)],...
        ['f_{X_{' int2str(iM(3))  '} \cdot X_' int2str(iM(2)) '}:' num2str(I3o2)]); legend('boxoff')
title(['causal decomposition of MIR between X_' int2str(iM(3)) ' and X_{' int2str(iM(2)) '}'])  




for N=4:M
    xtesto=0.1; ymax=max(OIRf{N}); ymin=min(OIRf{N}); dy=(ymax-ymin)/15;
    figure(N-2);
    subplot(1,2,1);
    plot(nfreq,OIRf{N-1},'r','linewidth',lw);hold on;
    text(xtesto,ymax,num2str(mean(OIRf{N-1})/2),'color','r')
    plot(nfreq,dOf{N-1},'b','linewidth',lw);
    text(xtesto,ymax-dy,num2str(mean(dOf{N-1})/2),'color','b')
    plot(nfreq,OIRf{N},'k','linewidth',lw);
    text(xtesto,ymax-2*dy,num2str(mean(OIRf{N})/2),'color','k')
    plot(nfreq,OIRf{N-1}+dOf{N-1},'k--','linewidth',lw) %verification: should see only one black line
    xlabel('\omega'); xlim([0 pi]);
    legend(['\Omega_{X_{' int2str(iM(1:N-1))  '}}:' num2str(OIR(N-1))],...
        ['\delta_{X_{' int2str(iM(1:N-1)) '};X_' int2str(iM(N)) '}:' num2str(dO(N-1))],...
        ['\Omega_{X_{' int2str(iM(1:N))  '}}:' num2str(OIR(N))]); legend('boxoff')
    title(['OIR recursion, from X_{' int2str(iM(1:N-1)) '} to X_{' int2str(iM(1:N)) '}'])
    
    xtesto=0.1; ymax=max(dOf{N-1}); ymin=min(dOf{N-1}); dy=(ymax-ymin)/15;
    subplot(1,2,2);
    plot(nfreq,dOf{N-1},'b','linewidth',lw); hold on;
    text(xtesto,ymax,num2str(mean(dOf{N-1})/2),'color','b')
    plot(nfreq,dO1_2f{N-1},'c','linewidth',lw);
    text(xtesto,ymax-dy,num2str(mean(dO1_2f{N-1})/2),'color','c')
    plot(nfreq,dO2_1f{N-1},'m','linewidth',lw);
    text(xtesto,ymax-2*dy,num2str(mean(dO2_1f{N-1})/2),'color','m')
    plot(nfreq,dO1o2f{N-1},'g','linewidth',lw);
    text(xtesto,ymax-3*dy,num2str(mean(dO1o2f{N-1})/2),'color','g')
    plot(nfreq,dO1_2f{N-1}+dO2_1f{N-1}+dO1o2f{N-1},'b--','linewidth',lw); %verification: should see only one blue line
    xlabel('\omega'); xlim([0 pi]);
    legend(['\delta_{X_{' int2str(iM(1:N-1)) '};X_' int2str(iM(N)) '}:' num2str(dO(N-1))],...
    ['\delta_{X_{' int2str(iM(1:N-1)) '} \rightarrow X_' int2str(iM(N)) '}:' num2str(dO1_2(N-1))],...
    ['\delta_{X_' int2str(iM(N)) ' \rightarrow X_{' int2str(iM(1:N-1)) '}}:' num2str(dO2_1(N-1))],...
    ['\delta_{X_{' int2str(iM(1:N-1)) '} \cdot X_' int2str(iM(N)) '}:' num2str(dO1o2(N-1))]); legend('boxoff')
    title(['causal decomposition of \Delta OIR between X_{' int2str(iM(1:N-1)) '} and X_{' int2str(iM(N)) '}'])  
end

disp(['iM = ' int2str(iM)])
