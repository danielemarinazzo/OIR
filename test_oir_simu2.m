%% test MIR/OIR decomposition on a 10-variate VAR model
clear; close all; clc;

Fs=100;
alpha=[8 12];
beta=[18 30];

%% SIMULATION 
[Am,Su,Ak,Mv,lambdamax] = VARsimu10; % VAR parameters
Q=size(Am,1);
M=length(Mv);


%% MIR analysis (without collecting the causal decompositions)
%%% equivalent ISS model
[A,C,K,~] = oir_ar2iss(Am); % Eq. (26)

disp('MI analysis...')
MIR=nan*ones(M,M); fij=cell(M,M); S=cell(M,M);
for i1=1:M
    for i2=i1+1:M
        outtmp=oir_mir(A,C,K,Su,Mv,i1,i2,Fs);
        S{i1,i2}=outtmp.S;
        MIR(i1,i2)=outtmp.I12; MIR(i2,i1)=MIR(i1,i2);
        fij{i1,i2}=outtmp.f12;
    end
end


%% OIR analysis (without collecting the causal decompositions)
allN=cell(M,1);
dO=cell(M,1); %deltaOIR
dOf=cell(M,1); %spectral deltaOIR 
for N=3:M
    comb=nchoosek(1:M,N); % all multiplets of size N
    allN{N}=comb;
        
    dO12=nan*ones(size(comb)); dO12f=cell(size(comb));
    for n=1:size(comb,1)
        ij=comb(n,:);
        for k=1:size(comb,2) % vary target inside the multiplet
            clc; disp(['multiplet ' int2str(n) ' of size ' int2str(size(comb,2)) ', target ' int2str(k)])
            out=oir_deltaO(A,C,K,Su,Mv,ij,ij(k),Fs);
            dO12(n,k)=out.dO12;
            dO12f{n,k}=out.dO12f; % Eq.(23)
        end
    end
    dO{N-1}=dO12;
    dOf{N-1}=dO12f;
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


%%% integral inside alpha and beta bands
nfreq=(out.freq*2*pi/Fs)'; %omega
freq=(out.freq)'; %omega
nfft=length(freq);

alpha_range=find(freq<=alpha(2) & freq>=alpha(1));
beta_range=find(freq<=beta(2) & freq>=beta(1));
OIRalpha=cell(M,1); OIRbeta=cell(M,1); %OIR in alpha and beta bands
for N=3:M
    for cnt=1:size(OIRf{N},1)
        OIRalpha{N}(cnt,1)=sum(OIRf{N}{cnt}(alpha_range))/(2*nfft);
        OIRbeta{N}(cnt,1)=sum(OIRf{N}{cnt}(beta_range))/(2*nfft);
    end
end

% mean(OIRf{N}{cnt})/2
% sum(OIRf{N}{cnt}(:))/(2*nfft)
% sum(OIRf{N}{cnt}(:))*((Fs/2-0)/nfft)/Fs
% sum(OIRf{N}{cnt}(alpha_range))/(2*nfft)
% (sum(OIRf{N}{cnt}(alpha_range))*(Fs/2-0)/nfft)/Fs
% mean(OIRf{N}{cnt}(alpha_range)) %questo non è corretto (non è nats)




%% plot OIR
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

%%%%%%%%%%%%% plot MIR and spectra
%%%%%% compute spectra for blocks not trivial code)
SZ=cell(M,M);
for i_1=1:M
    for i_2=i_1+1:M
        SZ1=nan*ones(length(freq),1); SZ2=SZ1;
        for n=1:length(freq)
            [~,~,j1,j2]=oir_subindexes(Mv,i_1,i_2);
            SZ1(n)=abs(det(S{i_1,i_2}(j1,j1,n)));
            SZ2(n)=abs(det(S{i_1,i_2}(j2,j2,n)));
        end
        SZ{i_1,i_2}(:,1)=SZ1;
        SZ{i_1,i_2}(:,2)=SZ2;
    end
end
set(gcf,'color','w');

%%% plot spectra and MIR
figure(1)
for i1=1:M
    subplot(M,M,M*(i1-1)+i1)
    if i1==M
        plot(freq,SZ{i1-1,i1}(:,2),'color',col3,'linewidth',lw);
    else
        plot(freq,SZ{i1,i1+1}(:,1),'color',col3,'linewidth',lw);
    end
    set(gca,'box','off') ; set(gca, 'color', 'none'); set(gca,'XTick',[], 'YTick', [])
    xlabel('f'); xlim([0 Fs/2]);
    for i2=i1+1:M
        subplot(M,M,M*(i1-1)+i2)
        plot(freq,fij{i1,i2},'color',col1b,'linewidth',lw);
        if max(fij{i1,i2})<0.0001, ylim([0 1]); end
        set(gca,'box','off') ; set(gca, 'color', 'none'); set(gca,'XTick',[], 'YTick', [])
        xlabel('f'); xlim([0 Fs/2]);
    end
end
set(gcf,'color','w');



%%%%%%%%%%% plot OIR
figure(2)
subplot(1,3,3);
plot(freq,OIRf{5}{1},'color',col1,'linewidth',lw);
xlabel('f'); xlim([0 Fs/2]);
stringa5=['O_{X_{' int2str(allN{5}) '}}:' num2str(OIR{5})];
legend(stringa5); title('OIR of order 5')

subplot(1,3,2);
colori2=[col1b; col2; col3; col4;col5]; stringa4=cell(size(allN{4},1),1);
for m=1:size(allN{4},1)
    plot(freq,OIRf{4}{m},'color',colori2(m,:),'linewidth',lw);hold on;
    stringa4{m}=['O_{X_{' int2str(allN{4}(m,:)) '}}'];
    %stringa4{m}=['O_{X_{' int2str(allN{4}(m,:)) '}}:' num2str(OIR{4}(m))];
end
xlabel('f'); xlim([0 Fs/2]);legend(stringa4); title('OIR of order 4')

subplot(1,3,1);
colori3=[colori2; col1c; col2b; col2c; col3b; col3c]; stringa3=cell(size(allN{3},1),1);
for m=1:size(allN{3},1)
    plot(freq,OIRf{3}{m},'color',colori3(m,:),'linewidth',lw);hold on;
    stringa3{m}=['O_{X_{' int2str(allN{3}(m,:)) '}}'];
    %stringa3{m}=['O_{X_{' int2str(allN{3}(m,:)) '}}:' num2str(OIR{3}(m))];
end
xlabel('f'); xlim([0 Fs/2]);legend(stringa3); title('OIR of order 3')

%%%% plot OIR in average and in bands
figure(3)
subplot(1,3,1)
bar([OIR{3} OIRalpha{3} OIRbeta{3}])
xticklabels(stringa3)
legend('whole band','\alpha band','\beta band')

subplot(1,3,2)
bar([OIR{4} OIRalpha{4} OIRbeta{4}])
xticklabels(stringa4)
legend('whole band','\alpha band','\beta band')

subplot(1,3,3)
bar([OIR{5} OIRalpha{5} OIRbeta{5}]')
xticklabels(stringa5)

set(gcf,'color','w');



