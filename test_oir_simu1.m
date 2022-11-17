%% test MIR/OIR decomposition on a modified version of the simulation of Faes et al, PTRSA 2021
clear; close all; clc;

%% SIMULATION setting
M=3; Mv=[1 1 1];
par.poles{1}=([0.7 0.1; 0.9 0.35]); % Self-oscillation 1 Oscillations RR
par.poles{2}=([0.7 0.1]); % Self-oscillation 2
par.poles{3}=([]); % Self-oscillations 3
ncoeff=20;
b1=fir1(ncoeff,2*0.2,'low');
b2=fir1(ncoeff,2*0.2,'high');
par.coup=[2 3 1 1];
c=0.4;
for k=1:ncoeff+1
    par.coup=[par.coup; 1 3 k 1*c*b1(k)];
    par.coup=[par.coup; 1 2 k 1*(1-c)*b2(k)];
end
par.Su=[2 0.5 2]; %variance of innovation processes

[Am,Su,Ak]=VARsimuPTRSA2021(M,par); % % VAR parameters
Am=Am'; Su=Su';

Fs=1;


%% analysis

%%% equivalent ISS model
[A,C,K,~] = oir_ar2iss(Am); % Eq. (26)
%%% MIR for the triplets
fij=cell(M,M); fi_j=cell(M,M); fj_i=cell(M,M); fioj=cell(M,M);
for i=1:M
    for j=i+1:M
        outtmp=oir_mir(A,C,K,Su,Mv,i,j);
        fi_j{i,j}=outtmp.f1_2;
        fi_j{j,i}=outtmp.f2_1;
        fij{i,j}=outtmp.f12;
        fioj{i,j}=outtmp.f1o2;
    end
end

%%% retrieve spectra
outY_X=oir_mir(A,C,K,Su,Mv,3,1);
S(:,3)=squeeze(abs(outY_X.S(1,1,:))); %SY
S(:,1)=squeeze(abs(outY_X.S(2,2,:))); %SX
outY_Z=oir_mir(A,C,K,Su,Mv,3,2);
S(:,2)=squeeze(abs(outY_Z.S(2,2,:))); %SZ

%%% deltaOIR(2) = OIR(3)
dOikjf=cell(M,1);dOik_jf=cell(M,1); dOj_ikf=cell(M,1); dOikojf=cell(M,1);
for j=1:M
    outtmp=oir_deltaO(A,C,K,Su,Mv,1:3,j);
    dOikjf{j}=outtmp.dO12f;
    dOik_jf{j}=outtmp.dO1_2f;
    dOj_ikf{j}=outtmp.dO2_1f;
    dOikojf{j}=outtmp.dO1o2f;
end

% out1=oir_deltaO(A,C,K,Su,Mv,1:3,1);
% dO23_1f=out1.dO1_2f; % Eq.(24)
% dO1_12f=out1.dO2_1f; % Eq.(24)
% dO23o1f=out1.dO1o2f; % Eq.(24)
% 
% out2=oir_deltaO(A,C,K,Su,Mv,1:3,2);
% dO13_2f=out2.dO1_2f; % Eq.(24)
% dO2_13f=out2.dO2_1f; % Eq.(24)
% dO13o2f=out2.dO1o2f; % Eq.(24)
% 
% out3=oir_deltaO(A,C,K,Su,Mv,1:3,3);
% dO12_3f=out3.dO1_2f; % Eq.(24)
% dO3_12f=out3.dO2_1f; % Eq.(24)
% dO13o3f=out3.dO1o2f; % Eq.(24)


%% plots
col1=[0 0 0]/255; % black
col1b=[96 96 96]/255; % dark grey
col1c=[192 192 192]/255; % light grey
col2=[192 64 64]/255; % dark red
col2b=[192 64 160]/255; % viola
col2c=[192 160 64]/255; % miele
col3=[64 128 192]/255; % dark azure
col3b=[64 64 192]/255; % blu
col3c=[64 192 192]/255; % celeste
col4=[64 160 64]/255; % dark green
col5=[224 224 32]/255; % yellow

lw=1.2;
nfreq=(outtmp.freq*2*pi/Fs)'; %omega
freq=(outtmp.freq)'; %omega
nfft=length(freq);
% ymin=min([fyxz;fyz;fyx;dOf]); ymax=max([fyxz;fyz;fyx;dOf]); deltay=(ymax-ymin)/20;


ymax=max([fij{1,2};fij{1,3};fij{2,3}]);
mappa=[nan 2 3; 5 nan 7; 9 10 nan];
figure(1)
for i=1:M
    subplot(M,M+1,5*(i-1)+1)
    plot(freq,S(:,i),'color',col1,'linewidth',lw);
    xlabel('f'); xlim([0 Fs/2]); legend(['S_X{_' int2str(i) '}'])
    for j=1:M
        if i~=j
            subplot(M,M+1,mappa(i,j))
            plot(freq,fi_j{i,j},'color',col2,'linewidth',lw); hold on;
            xlabel('f'); xlim([0 Fs/2]); ylim([0 2]);
            if i==1 && j==2, ylim([0 ymax]); end
            if i<j
                legend(['f_{X_' int2str(i) '\rightarrow X{_' int2str(j) '}}'])
            end
        end
        if i>j
            plot(freq,fioj{j,i},'color',col3,'linewidth',lw,'linestyle','--');
            legend(['f_{X_' int2str(i) '\rightarrow X{_' int2str(j) '}}'],...
                ['f_{X_' int2str(i) '\cdot X{_' int2str(j) '}}']);
            ylim([-0.02 2]);
        end
    end
end
subplot(M,M+1,[4 8 12]);
j=1;
plot(freq,dOikjf{j},'color',col3,'linewidth',lw); hold on;
plot(freq,dOik_jf{j},'color',col2,'linewidth',lw);
plot(freq,dOj_ikf{j},'color',col4,'linewidth',lw,'linestyle','--');
plot(freq,dOikojf{j},'color',col5,'linewidth',lw,'linestyle','--');
xlim([0 Fs/2]);
legend('\nu_{X_{123}}',...
    ['\delta_{X_{' int2str(setdiff(1:3,j)) '}\rightarrow X_' int2str(j) '}'],...
    ['\delta_{X_{' int2str(j) '}\rightarrow X_{' int2str(setdiff(1:3,j)) '}}'],...
    ['\delta_{X_{' int2str(setdiff(1:3,j)) '}\cdot X_' int2str(j) '}']...
    );
set(gcf,'color','w');

% figure(2)
% for j=1:M
%     subplot(1,M,j);
%     plot(freq,dOikjf{j},'color',col3,'linewidth',lw); hold on;
%     plot(freq,dOik_jf{j},'color',col2,'linewidth',lw);
%     plot(freq,dOj_ikf{j},'color',col4,'linewidth',lw,'linestyle','--');
%     plot(freq,dOikojf{j},'color',col5,'linewidth',lw,'linestyle','--');
%     legend('\nu_{X_{123}}',...
%         ['\delta_{X_{' int2str(setdiff(1:3,j)) '}\rightarrow X_' int2str(j) '}'],...
%         ['\delta_{X_{' int2str(j) '}\rightarrow X_{' int2str(setdiff(1:3,j)) '}}'],...
%         ['\delta_{X_{' int2str(setdiff(1:3,j)) '}\cdot X_' int2str(j) '}']...
%         );
% end

[h1,f]=freqz(b1,1,2*nfft,'whole',Fs);
[h2]=freqz(b2,1,2*nfft,'whole',Fs);
H1=abs(h1); H1=c*H1(1:nfft); %consider only up to Nyquist frequency
H2=abs(h2); H2=(1-c)*H2(1:nfft);
figure(3);
subplot(1,2,1);plot(freq,H1,'linewidth',lw)
xlim([0 Fs/2]);ylim([0 1]);legend('H_{31}'); xlabel('frequency')
subplot(1,2,2);plot(freq,H2,'linewidth',lw)
xlim([0 Fs/2]);ylim([0 1]);legend('H_{21}'); xlabel('frequency')
set(gcf,'color','w');

% figure(89); plot(f,angle(h1));

for j=1:M
    outtmp=oir_deltaO(A,C,K,Su,Mv,1:3,j);
    dOikjf{j}=outtmp.dO12f;
    dOik_jf{j}=outtmp.dO1_2f;
    dOj_ikf{j}=outtmp.dO2_1f;
    dOikojf{j}=outtmp.dO1o2f;
end


disp(['I_{X1,X2} = ' num2str(mean(fij{1,2})/2)]);
disp(['I_{X1,X3} = ' num2str(mean(fij{1,3})/2)]);
disp(['I_{X2,X3} = ' num2str(mean(fij{2,3})/2)]);
disp(' ')
disp(['T_{X1->X2} = ' num2str(mean(fi_j{1,2})/2)]);
disp(['T_{X1->X3} = ' num2str(mean(fi_j{1,3})/2)]);
disp(['T_{X2->X3} = ' num2str(mean(fi_j{2,3})/2)]);
disp(' ')
disp(['Omega_{X1,X2,X3} = ' num2str(mean(dOikjf{1})/2)]);
disp(['Delta_{X2,X3}->X1 = ' num2str(mean(dOik_jf{1})/2)]);
disp(['Delta_X1->{X2,X3} = ' num2str(mean(dOj_ikf{1})/2)]);
disp(['Delta_{X2,X3}.X1 = ' num2str(mean(dOikojf{1})/2)]);

%%%% correggere con le bande LF, HF...
dOikjf_LF=mean(dOikjf{1}(freq>0.04 & freq<0.12));
dOikjf_HF=mean(dOikjf{1}(freq>0.31 & freq<0.39));
disp(' ')
disp(['Omega_{X1,X2,X3}(LF) = ' num2str(dOikjf_LF)]);
disp(['Omega_{X1,X2,X3}(HF) = ' num2str(dOikjf_HF)]);




