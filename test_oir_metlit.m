%% OIR analysis on ECoG signals
close all;clear all;clc
datadir='C:\Users\daniele\Dropbox\code\MEG_ASSR\sample_data';
sub=1;
ROIs=[1 5 31 41 55 61];
%% load data
load([datadir filesep 'COB' num2str(sub,'%2.2i') '_Lit_adj_Corr_sample_data.mat']);
lit_tot=sample_data(ROIs,:,:);
ntrials_lit=size(lit_tot,3);
load([datadir filesep 'COB' num2str(sub,'%2.2i') '_Met_adj_Corr_sample_data.mat']);
met_tot=sample_data(ROIs,:,:);
ntrials_met=size(met_tot,3);
min_trials=min(ntrials_lit,ntrials_met);
lit_tot=lit_tot(:,:,1:min_trials);
met_tot=met_tot(:,:,1:min_trials);

met=squeeze(met_tot(:,:,1))'; met=zscore(detrend(met),[],1);
lit=squeeze(lit_tot(:,:,1))'; lit=zscore(detrend(lit),[],1);
%% parameters
Fs=600; % sampling frequency
nfft=1000;
%% OIR analysis
Mv=[1 1 1 1 1 1]; %number of series in each block
momax=16;% maximum model order
M=length(Mv);
allN=cell(M,1);
dO_lit=cell(M,1); %deltaOIR - lit
dOf_lit=cell(M,1); %spectral deltaOIR - lit
dO_met=cell(M,1); %deltaOIR -met
dOf_met=cell(M,1); %spectral deltaOIR - met

for N=3:M
    
    comb=nchoosek(1:M,N); % all multiplets of size N
    allN{N}=comb;
    DO_lit=cell(size(comb)); DOF_lit=cell(size(comb));
    DO_met=cell(size(comb)); DOF_met=cell(size(comb));
    
    for n=1:size(comb,1)
        ij=comb(n,:);
        for k=1:size(comb,2) % vary target inside the multiplet
            clc; disp(['multiplet ' int2str(n) ' of size ' int2str(size(comb,2)) ', target ' int2str(k)])
            %%% model order estimation
            [~,pottmdl_lit] = oir_mosVAR(lit',momax);
            [~,pottmdl_met] = oir_mosVAR(met',momax);
            %%% VAR identification procedure
            [Am_r,Su_lit]=idMVAR(lit',pottmdl_lit,0);
            [Am_a,Su_met]=idMVAR(met',pottmdl_met,0);
            %%% equivalent ISS model
            [A_lit,C_lit,K_lit,~] = oir_ar2iss(Am_r);
            [A_met,C_met,K_met,~] = oir_ar2iss(Am_a);
            %%% delta OIR (lit)
            out_lit=oir_deltaO(A_lit,C_lit,K_lit,Su_lit,Mv,ij,ij(k),Fs,nfft);
            DOf_lit=out_lit.dO12f;
            DO_lit=out_lit.dO12;
            %%% delta OIR (met)
            out_met=oir_deltaO(A_met,C_met,K_met,Su_met,Mv,ij,ij(k),Fs,nfft);
            DOf_met=out_met.dO12f;
            DO_met=out_met.dO12;
            freq=(out_lit.freq)';
            
            %%% collecting all the measures
            DOF_lit{n,k} = DOf_lit;
            DO_lit{n,k} = DO_lit;
            DOF_met{n,k} = DOf_met;
            DO_met{n,k} = DO_met; 
        end
    end
    dO_lit{N-1}=DO_lit;
    dOf_lit{N-1}=DOF_lit;
    dO_met{N-1}=DO_met;
    dOf_met{N-1}=DOF_met;
end

OIR_lit=cell(M,1); %OIR
OIRf_met=cell(M,1); %spectral OIR
OIR_lit{3}=dO_lit{2}(:,1); %other columns are the same
OIRf_lit{3}=dOf_lit{2}(:,1);
OIR_met{3}=dO_met{2}(:,1); %other columns are the same
OIRf_met{3}=dOf_met{2}(:,1);

for N=4:M
    t1=1; %index of Xj (actually can be any number btw 1 and N)
    t2=setdiff(1:N,t1); %index of X-j
    for cnt=1:size(allN{N},1)
        tmp=allN{N}(cnt,:);
        iN1=find(sum(tmp(t2)==allN{N-1},2)==N-1); %position of X-j in the O-info one step back
        OIR_lit{N}{cnt,1}=OIR_lit{N-1}{iN1}+dO_lit{N-1}{cnt,t1};
        OIRf_lit{N}{cnt,1}=OIRf_lit{N-1}{iN1}+dOf_lit{N-1}{cnt,t1};
        OIR_met{N}{cnt,1}=OIR_met{N-1}{iN1}+dO_met{N-1}{cnt,t1};
        OIRf_met{N}{cnt,1}=OIRf_met{N-1}{iN1}+dOf_met{N-1}{cnt,t1};
    end
end

%% plot multiplets of order 3 (Figure 5)
figure
stringa3=cell(size(allN{3},1),1);
IND=[1,8]; %plot only specific combination as in figure 5
for m=1:size(IND,2)
    subplot(1,2,m);
    %%% lit
    plot(freq,mean(OIRf_lit{3}{IND(m)},2),'Color','r','LineWidth',1.5);
    hold on
    %%% metthesia
    plot(freq,mean(OIRf_met{3}{IND(m)},2),'Color','b','LineWidth',1.5);
    xlim([0 70])
    set(gca,'Xtick',[0:10:150])
    xlabel('f');
    stringa3{m}=['\nu_{X_{' int2str(allN{3}(IND(m),:)) '}}'];
    title(stringa3{m})
end
sgtitle('OIR of order 3')
legend({'lit','met'})
%% plot multiplets of order 4
figure
stringa4=cell(size(allN{4},1),1);
IND=[1,5]; %plot only specific combination as in figure 5
for m=1:size(IND,2)
    subplot(1,2,m);
    plot(freq,mean(OIRf_lit{4}{IND(m)},2),'Color','r','LineWidth',1.5);
    hold on
    %%% metthesia
    plot(freq,mean(OIRf_met{4}{IND(m)},2),'Color','b','LineWidth',1.5);
    xlim([0 70])
    set(gca,'Xtick',[0:10:150])
    xlabel('f');
    stringa4{m}=['\nu_{X_{' int2str(allN{4}(IND(m),:)) '}}'];
    title(stringa4{m})
end
sgtitle('OIR of order 4')
legend({'lit','met'})
%% plot multiplets of order 5
clear h
figure
m=1;
plot(freq,mean(OIRf_lit{5}{m},2),'Color','r','LineWidth',1.5);
%%% metthesia
hold on
plot(freq,mean(OIRf_met{5}{m},2),'Color','b','LineWidth',1.5);
xlim([0 70])
set(gca,'Xtick',[0:10:150])
xlabel('f');
stringa5{m}=['\nu_{X_{' int2str(allN{5}(m,:)) '}}'];
title(stringa5{m})
legend({'lit','met'})
