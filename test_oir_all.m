%% test script for the O-information rate - complete and decomposed
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

Mv=[1 1 2 2 2]; % override Mv to have more blocks

M=length(Mv);

%% MIR analysis
%%% equivalent ISS model
[A,C,K,~] = oir_ar2iss(Am); % Eq. (26)

disp('MI analysis...')
Iij=nan*ones(M,M); Ti_j=nan*ones(M,M); Iioj=nan*ones(M,M);
fij=cell(M,M); fi_j=cell(M,M); fioj=cell(M,M);
for i1=1:M
    for i2=i1+1:M
        outtmp=oir_mir(A,C,K,Su,Mv,i1,i2,Fs,nfft,do_td);
        Iij(i1,i2)=outtmp.I12; Iij(i2,i1)=Iij(i1,i2);
        Ti_j(i1,i2)=outtmp.T1_2;  Ti_j(i2,i1)=outtmp.T2_1;
        Iioj(i1,i2)=outtmp.I1o2; Iioj(i2,i1)=Iioj(i1,i2); 
        fij{i1,i2}=outtmp.f12;
        fi_j{i1,i2}=outtmp.f1_2;  fi_j{i2,i1}=outtmp.f2_1;
        fioj{i1,i2}=outtmp.f1o2; 
    end
end

%% OIR analysis
allN=cell(M,1);
dO=cell(M,1); %deltaOIR
dOd=cell(M,1); %deltaOIR decomposed (three terms)
dOf=cell(M,1); %spectral deltaOIR 
dOfd=cell(M,1); %spectral deltaOIR decomposed (three terms)
for N=3:M
    comb=nchoosek(1:M,N); % all multiplets of size N
    allN{N}=comb;
        
    dO12=nan*ones(size(comb)); dO1_2=dO12; dO2_1=dO12; dO1o2=dO12;
    dO12f=cell(size(comb)); dO1_2f=dO12f; dO2_1f=dO12f; dO1o2f=dO12f;
    for n=1:size(comb,1)
        ij=comb(n,:);
        for k=1:size(comb,2) % vary target inside the multiplet
            clc; disp(['multiplet ' int2str(n) ' of size ' int2str(size(comb,2)) ', target ' int2str(k)])
            out=oir_deltaO(A,C,K,Su,Mv,ij,ij(k),Fs,nfft,do_td);
            dO12(n,k)=out.dO12;
            dO1_2(n,k)=out.dO1_2; % Eq.(8a)
            dO2_1(n,k)=out.dO2_1; % Eq.(8b)
            dO1o2(n,k)=out.dO1o2; % Eq.(8c)
            dO12f{n,k}=out.dO12f; % Eq.(23)
            dO1_2f{n,k}=out.dO1_2f; % Eq.(24)
            dO2_1f{n,k}=out.dO2_1f; % Eq.(24)
            dO1o2f{n,k}=out.dO1o2f; % Eq.(24)
        end
    end
    dO{N-1}=dO12;
    dOd{N-1}.dO1_2=dO1_2;
    dOd{N-1}.dO2_1=dO2_1;
    dOd{N-1}.dO1o2=dO1o2;
    dOf{N-1}=dO12f;
    dOfd{N-1}.dO1_2f=dO1_2f;
    dOfd{N-1}.dO2_1f=dO2_1f;
    dOfd{N-1}.dO1o2f=dO1o2f;
end


OIR=cell(M,1); %OIR
OIRf=cell(M,1); %spectral OIR
OIR{3}=dO{2}(:,1); %other columns are the same
OIRf{3}=dOf{2}(:,1);
for N=4:M
    t1=1; %index of Xj (actually can be any number btw 1 and N)
    t2=setdiff(1:N,t1); %index of X-j
    for cnt=1:size(allN{N},1)
        tmp=allN{N}(cnt,:);
        iN1=find(sum(tmp(t2)==allN{N-1},2)==N-1); %position of X-j in the O-info one step back
        OIR{N}(cnt,1)=OIR{N-1}(iN1)+dO{N-1}(cnt,t1);
        OIRf{N}{cnt,1}=OIRf{N-1}{iN1}+dOf{N-1}{cnt,t1};
    end
end

%% plots for one representative sequence of processes
close all;
rng('shuffle');% to generate different permutation at each run of the script
% iM=randperm(M); %random ordering of variable indexes - show one random path
iM=[5 1 4 3 2];

lw=1.2;
nfreq=(outtmp.freq*2*pi)'; %omega



for N=4:M
    iN=find(sum((sort(iM(1:N))==allN{N}),2)==N); % position of iM(1:N) in allN{N}
    iN1= find(sum(sort(iM(1:N-1))==allN{N-1},2)==N-1); % position of iM(1:N-1) in allN{N-1}
    iD=find(allN{N}(iN,:)==iM(N)); % position of iM(N) in allN{N}(iN,:)
    
    ONsel=OIR{N}(iN); % selected OIR^N to print
    ON1sel=OIR{N-1}(iN1); % selected OIR^(N-1) to print
    dOsel=dO{N-1}(iN,iD); % selected delta OIR to print
    dO1_2sel=dOd{N-1}.dO1_2(iN,iD);
    dO2_1sel=dOd{N-1}.dO2_1(iN,iD);
    dO1o2sel=dOd{N-1}.dO1o2(iN,iD);
    
    ONfsel=OIRf{N}{iN}; % selected OIR^N(w) to plot
    ON1fsel=OIRf{N-1}{iN1}; % selected OIR^(N-1)(w) to plot
    dOfsel=dOf{N-1}{iN,iD}; % selected delta OIR(w) to plot
    dOf1_2sel=dOfd{N-1}.dO1_2f{iN,iD};
    dOf2_1sel=dOfd{N-1}.dO2_1f{iN,iD};
    dOf1o2sel=dOfd{N-1}.dO1o2f{iN,iD};
    
    xtesto=0.3;
    figure(N-2);
    subplot(1,2,1);
    plot(nfreq,ON1fsel,'r','linewidth',lw);hold on;
    text(xtesto,1.02*max(ON1fsel),num2str(mean(ON1fsel)/2),'color','r')
    plot(nfreq,dOfsel,'b','linewidth',lw);
    text(2*xtesto,1.02*max(dOfsel),num2str(mean(dOfsel)/2),'color','b')
    plot(nfreq,ONfsel,'k','linewidth',lw);
    text(3*xtesto,1.02*max(ONfsel),num2str(mean(ONfsel)/2),'color','k')
    plot(nfreq,ON1fsel+dOfsel,'k--','linewidth',lw) %verification: should see only one black line
    xlabel('\omega'); xlim([0 pi]);
    legend(['\Omega_{X_{' int2str(iM(1:N-1))  '}}:' num2str(ON1sel)],...
        ['\delta_{X_{' int2str(iM(1:N-1)) '};X_' int2str(iM(N)) '}:' num2str(dOsel)],...
        ['\Omega_{X_{' int2str(iM(1:N))  '}}:' num2str(ONsel)]);
    subplot(1,2,2);
    plot(nfreq,dOfsel,'b','linewidth',lw); hold on;
    text(xtesto,1.02*max(dOfsel),num2str(mean(dOfsel)/2),'color','b')
    plot(nfreq,dOf1_2sel,'c','linewidth',lw);
    text(2*xtesto,1.02*max(dOf1_2sel),num2str(mean(dOf1_2sel)/2),'color','c')
    plot(nfreq,dOf2_1sel,'m','linewidth',lw);
    text(3*xtesto,1.02*max(dOf2_1sel),num2str(mean(dOf2_1sel)/2),'color','m')
    plot(nfreq,dOf1o2sel,'g','linewidth',lw);
    text(4*xtesto,1.02*max(dOf1o2sel),num2str(mean(dOf1o2sel)/2),'color','g')
    plot(nfreq,dOf1_2sel+dOf2_1sel+dOf1o2sel,'b--','linewidth',lw); %verification: should see only one blue line
    xlabel('\omega'); xlim([0 pi]);
    legend(['\delta_{X_{' int2str(iM(1:N-1)) '};X_' int2str(iM(N)) '}:' num2str(dOsel)],...
    ['\delta_{X_{' int2str(iM(1:N-1)) '} \rightarrow X_' int2str(iM(N)) '}:' num2str(dO1_2sel)],...
    ['\delta_{X_' int2str(iM(N)) ' \rightarrow X_{' int2str(iM(1:N-1)) '}}:' num2str(dO2_1sel)],...
    ['\delta_{X_{' int2str(iM(1:N-1)) '} \cdot X_' int2str(iM(N)) '}:' num2str(dO1o2sel)]);
end
