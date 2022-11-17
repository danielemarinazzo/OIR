%% OIR decomposition in time and spectral domain for a single subject
clear; close all; clc;

%% Data analysis parameters
nfft=1000; % number of points on frequency axis 
pfilter=0.94; % high-pass filter parameter for detrending
c=4; % number of experimental conditions

% selection of blocks and order for the computation of OIR 
Mv=[1 1 1]; % number of series in each block --> 3 blocks, each containing 1 series
iM=[1 2 3]; % indexes of the series to analyze
j=3; % index of the series target

% spectral analysis
rangeLF = [0.04 0.12]; % fixed
HFwidth=0.04; % HF bandwidth 
fresp=0.1645; % computed a priori


%% Load data
load('data_HPRESPSAP.txt');

series{1}=data_HPRESPSAP(:,(1:3)); %SB
series{2}=data_HPRESPSAP(:,(4:6)); %CB10
series{3}=data_HPRESPSAP(:,(7:9)); %CB15
series{4}=data_HPRESPSAP(:,(10:12)); %CB20

%% Preprocessing and data analysis
% init
dO1_2=cell(1,c); dO2_1=dO1_2; dO1o2=dO1_2; dO12=dO1_2;
dO1_2f=dO1_2; dO2_1f=dO1_2; dO1o2f=dO1_2; dO12f=dO1_2;
dO1_2f_HF=dO1_2; dO2_1f_HF=dO1_2; dO1o2f_HF=dO1_2; dO12f_HF=dO1_2;
dO1_2f_LF=dO1_2; dO2_1f_LF=dO1_2; dO1o2f_LF=dO1_2; dO12f_LF=dO1_2;
for ic=1:c % cycle over 4 conditions

    So=series{ic}; %1st series: HP; 2nd series: RESP; 3rd series: SAP
    [N,Q]=size(So);
    % rearrange matrix to have RESP=X1; SAP=X2; RR=X3:
    So=So(:,[2 3 1]);

    % mean value of Heart Period
    mean_HP=mean(So(:,3)); 
    
    % sampling frequency 
    fs=1/mean_HP;
    
    % remove mean
    S=So-mean(So);
    
    % VAR identification
    p(ic) = oir_mosVAR(S',10); % model order selection (Akaike)
    [Am,Su]=oir_idVAR(S',p(ic)); 
    
    % AR to ISS model
    [A,C,K,~] = oir_ar2iss(Am);
    
    % computation of the HF band around the respiratory peak
    nlHFwidth=round((nfft*2/fs)*HFwidth); % from width [Hz] to freq. points
    nresp=round((nfft*2/fs)*fresp); % retrieve frequency bin of the respiratory peak
    spectBandHF=[nresp-nlHFwidth nresp+nlHFwidth];
    if spectBandHF(2)>nfft, spectBandHF(2)=nfft; end
    
    % fixed LF band:
    spectBandLF=round((nfft*2/fs)*rangeLF) + [1 0];
    
    %% Computation of OIR measures
    out_OIR=oir_deltaO(A,C,K,Su,Mv,iM,j,fs,nfft);
    f=out_OIR.freq;

    % time domain measures:
    dO1_2{ic}=out_OIR.dO1_2; 
    dO2_1{ic}=out_OIR.dO2_1; 
    dO1o2{ic}=out_OIR.dO1o2;
    dO12{ic}=out_OIR.dO12;
    % spectral distributions:
    dO1_2f{ic}=out_OIR.dO1_2f; 
    dO2_1f{ic}=out_OIR.dO2_1f; 
    dO1o2f{ic}=out_OIR.dO1o2f;
    dO12f{ic}=out_OIR.dO12f;
    
    % average values within the HF band
    dO12f_HF{ic}=sum(dO12f{ic}(spectBandHF(1):spectBandHF(2)))/(2*nfft);
    dO1_2f_HF{ic}=sum(dO1_2f{ic}(spectBandHF(1):spectBandHF(2)))/(2*nfft);
    dO2_1f_HF{ic}=sum(dO2_1f{ic}(spectBandHF(1):spectBandHF(2)))/(2*nfft);
    dO1o2f_HF{ic}=sum(dO1o2f{ic}(spectBandHF(1):spectBandHF(2)))/(2*nfft);
    
    % average values within the LF band
    dO12f_LF{ic}=sum(dO12f{ic}(spectBandLF(1):spectBandLF(2)))/(2*nfft);
    dO1_2f_LF{ic}=sum(dO1_2f{ic}(spectBandLF(1):spectBandLF(2)))/(2*nfft);
    dO2_1f_LF{ic}=sum(dO2_1f{ic}(spectBandLF(1):spectBandLF(2)))/(2*nfft);
    dO1o2f_LF{ic}=sum(dO1o2f{ic}(spectBandLF(1):spectBandLF(2)))/(2*nfft);

end

%% Plot spectral distributions
% parameters for figure 
col{1}=[192 0 0]./255; col{2}=[0 80 150]./255;
col{3}=[60 180 200]./255; col{4}=[10 255 250]./255;
tti{1}='SB'; tti{2}='CB10'; tti{3}='CB15'; tti{4}='CB20';
DimensioneFont=15.5;
lw=1.2;

h0=figure('numbertitle','off','WindowState','maximized','Color','w');

figure(h0)
for m=1:c
    s{m}=subplot(3,4,m);
    plot(f,dO12f{m},'color',col{m},'LineWidth',lw);
    xlim([0 fs/2]); ylim([-0.2 4]); title(tti{m}); 
    if m==1
        ylabel(['\delta_{X_{' int2str(iM(3))  '};X_{' int2str(iM(1)) '},X_{' int2str(iM(2)) '}}'])
        yticks([0 1 2 3 4])
    else
        yticks([ ]);
    end
    xticks([])
    set(gca,'FontSize',DimensioneFont,'FontName','Times');

    s{m+c}=subplot(3,4,m+c);
    plot(f,dO1_2f{m}+dO1o2f{m},'color',col{m},'LineWidth',lw);
    xlim([0 fs/2]); ylim([-1 3]); 
    if m==1
        ylabel(['\delta_{X_{' int2str(iM(1))  '} , X_{' int2str(iM(2)) '}\rightarrow\cdotX_{' int2str(iM(3)) '}}'])
        yticks([-1 0 1 2])
    else
        yticks([ ]);
    end
    xticks([])
    set(gca,'FontSize',DimensioneFont,'FontName','Times');

    s{m+2*c}=subplot(3,4,m+2*c);
    plot(f,dO2_1f{m},'color',col{m},'LineWidth',lw);
    if m==1
        ylabel(['\delta_{X_{' int2str(iM(3))  '}\rightarrowX_{' int2str(iM(1)) '},X_{' int2str(iM(2)) '}}'])
        yticks([-1 0 1 2])
    else
        yticks([ ]);
    end
    xlabel('Freq. [Hz]'); xlim([0 fs/2]); ylim([-1 3.1]);
    set(gca,'FontSize',DimensioneFont,'FontName','Times');

end

set(s{1},'Position',[0.1300 0.67 0.2 0.27])
set(s{2},'Position',[0.3361 0.67 0.2 0.27])
set(s{3},'Position',[0.542234042553192	0.67 0.2 0.27])
set(s{4},'Position',[0.748351063829787	0.67 0.2 0.27])
set(s{5},'Position',[0.130000000000000	0.39 0.2 0.27])
set(s{6},'Position',[0.336117021276596	0.39 0.2 0.27])
set(s{7},'Position',[0.542234042553192	0.39 0.2 0.27])
set(s{8},'Position',[0.748351063829787	0.39 0.2 0.27])
set(s{9},'Position',[0.130000000000000	0.11 0.2 0.27])
set(s{10},'Position',[0.336117021276596	0.11 0.2 0.27])
set(s{11},'Position',[0.542234042553192	0.11 0.2 0.27])
set(s{12},'Position',[0.748351063829787	0.11 0.2 0.27])

