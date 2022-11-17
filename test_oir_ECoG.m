%% OIR analysis on ECoG signals
close all;clear all;clc

%% load data
load('ECoG_data.mat');
%% parameters
Fs=250; % sampling frequency
nfft=1000;
%% OIR analysis
Mv=[2 2 2 2 2]; %number of series in each block
momax=16;% maximum model order
M=length(Mv);
allN=cell(M,1);
dO_r=cell(M,1); %deltaOIR - rest
dOf_r=cell(M,1); %spectral deltaOIR - rest
dO_a=cell(M,1); %deltaOIR -anestesia
dOf_a=cell(M,1); %spectral deltaOIR - anestesia

for N=3:M
    
    comb=nchoosek(1:M,N); % all multiplets of size N
    allN{N}=comb;
    DO_rest=cell(size(comb)); DOF_rest=cell(size(comb));
    DO_anes=cell(size(comb)); DOF_anes=cell(size(comb));
    
    for n=1:size(comb,1)
        ij=comb(n,:);
        for k=1:size(comb,2) % vary target inside the multiplet
            clc; disp(['multiplet ' int2str(n) ' of size ' int2str(size(comb,2)) ', target ' int2str(k)])
            %%% model order estimation
            [~,pottmdl_r] = oir_mosVAR(REST',momax);
            [~,pottmdl_a] = oir_mosVAR(ANES',momax);
            %%% VAR identification procedure
            [Am_r,Su_r]=idMVAR(REST',pottmdl_r,0);
            [Am_a,Su_a]=idMVAR(ANES',pottmdl_a,0);
            %%% equivalent ISS model
            [A_r,C_r,K_r,~] = oir_ar2iss(Am_r);
            [A_a,C_a,K_a,~] = oir_ar2iss(Am_a);
            %%% delta OIR (rest)
            out_r=oir_deltaO(A_r,C_r,K_r,Su_r,Mv,ij,ij(k),Fs,nfft);
            DOF_r=out_r.dO12f;
            DO_r=out_r.dO12;
            %%% delta OIR (anes)
            out_a=oir_deltaO(A_a,C_a,K_a,Su_a,Mv,ij,ij(k),Fs,nfft);
            DOF_a=out_a.dO12f;
            DO_a=out_a.dO12;
            freq=(out_r.freq)';
            
            %%% collecting all the measures
            DOF_rest{n,k} = DOF_r;
            DO_rest{n,k} = DO_r;
            DOF_anes{n,k} = DOF_a;
            DO_anes{n,k} = DO_a; 
        end
    end
    dO_r{N-1}=DO_rest;
    dOf_r{N-1}=DOF_rest;
    dO_a{N-1}=DO_anes;
    dOf_a{N-1}=DOF_anes;
end

OIR_r=cell(M,1); %OIR
OIRf_a=cell(M,1); %spectral OIR
OIR_r{3}=dO_r{2}(:,1); %other columns are the same
OIRf_r{3}=dOf_r{2}(:,1);
OIR_a{3}=dO_a{2}(:,1); %other columns are the same
OIRf_a{3}=dOf_a{2}(:,1);

for N=4:M
    t1=1; %index of Xj (actually can be any number btw 1 and N)
    t2=setdiff(1:N,t1); %index of X-j
    for cnt=1:size(allN{N},1)
        tmp=allN{N}(cnt,:);
        iN1=find(sum(tmp(t2)==allN{N-1},2)==N-1); %position of X-j in the O-info one step back
        OIR_r{N}{cnt,1}=OIR_r{N-1}{iN1}+dO_r{N-1}{cnt,t1};
        OIRf_r{N}{cnt,1}=OIRf_r{N-1}{iN1}+dOf_r{N-1}{cnt,t1};
        OIR_a{N}{cnt,1}=OIR_a{N-1}{iN1}+dO_a{N-1}{cnt,t1};
        OIRf_a{N}{cnt,1}=OIRf_a{N-1}{iN1}+dOf_a{N-1}{cnt,t1};
    end
end

%% plot multiplets of order 3 (Figure 5)
figure
stringa3=cell(size(allN{3},1),1);
IND=[1,8]; %plot only specific combination as in figure 5
for m=1:size(IND,2)
    subplot(1,2,m);
    %%% rest
    plot(freq,mean(OIRf_r{3}{IND(m)},2),'Color','r','LineWidth',1.5);
    hold on
    %%% anesthesia
    plot(freq,mean(OIRf_a{3}{IND(m)},2),'Color','b','LineWidth',1.5);
    xlim([0 70])
    set(gca,'Xtick',[0:10:150])
    xlabel('f');
    stringa3{m}=['\nu_{X_{' int2str(allN{3}(IND(m),:)) '}}'];
    title(stringa3{m})
end
sgtitle('OIR of order 3')
legend({'rest','anes'})
%% plot multiplets of order 4
figure
stringa4=cell(size(allN{4},1),1);
IND=[1,5]; %plot only specific combination as in figure 5
for m=1:size(IND,2)
    subplot(1,2,m);
    plot(freq,mean(OIRf_r{4}{IND(m)},2),'Color','r','LineWidth',1.5);
    hold on
    %%% anesthesia
    plot(freq,mean(OIRf_a{4}{IND(m)},2),'Color','b','LineWidth',1.5);
    xlim([0 70])
    set(gca,'Xtick',[0:10:150])
    xlabel('f');
    stringa4{m}=['\nu_{X_{' int2str(allN{4}(IND(m),:)) '}}'];
    title(stringa4{m})
end
sgtitle('OIR of order 4')
legend({'rest','anes'})
%% plot multiplets of order 5
clear h
figure
m=1;
plot(freq,mean(OIRf_r{5}{m},2),'Color','r','LineWidth',1.5);
%%% anesthesia
hold on
plot(freq,mean(OIRf_a{5}{m},2),'Color','b','LineWidth',1.5);
xlim([0 70])
set(gca,'Xtick',[0:10:150])
xlabel('f');
stringa5{m}=['\nu_{X_{' int2str(allN{5}(m,:)) '}}'];
title(stringa5{m})
legend({'rest','anes'})
