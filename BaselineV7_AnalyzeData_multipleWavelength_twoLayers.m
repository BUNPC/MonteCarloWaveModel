%% read history file
% clc
% clear
% jobid=1;


config=1; % configuration number
nDets=7;
Ntissue=2;
filenm='Tissue1';
fid = fopen([filenm,'.his'], 'rb');
his = fread(fid, inf, 'float');
totalP=length(his);
VariableNumber=3;
his=reshape(his,[1+Ntissue*VariableNumber,totalP/(1+Ntissue*VariableNumber)]);
% calculate g1 and g2
prob_cell=1;
D=10^(-6); %%mm^2/s
Dtissue=[D,D];% diffusion coefficient in the tissue types

k0=2*pi/(0.8*10^(-3));

%figure
%k0=2*pi/(0.8*10^(-3)); %wavelength
tmua=[0.01,0.01];
tall=0:10^(-6):10^(-2); % s

mu_a=0.01; %mm^(-1)
l_coherence=90;
load('FWHM_NearGaussian_320ps.mat'); % the IRF function we are using

t_s0all=0.5:1:2.5;

%%
    for ts0=1:length(t_s0all)

       t_s0=t_s0all(ts0); %measure at 1 ns

    clear mask Nscatteringall pathlengthall
    %tau=10.^([-7:0.1:-3]);%:(0:20)*10; %s
    %tau=10.^([-6:0.1:0.1]);
    %Ntau=length(tau);
   
 for ii=3 % this determines the source-detector separation we are using. If only one detector is used in the Monte carlo simulation, ii=1

     index=find(his(1,:)==ii-1);
     Nphotons=length(index);
     %NphotonsAll(ii)=Nphotons;
     %g1=ones(Nphotons,Ntau);
     %g1_original=ones(Nphotons,Ntau);
      %g1_cellSwelling=ones(Nphotons,Ntau);
      clear momentumTransfer momentumTrasfer_cell MomentumTransferall
      MomentumTransferall=zeros(Ntissue,Nphotons);
      Nscatteringall=zeros(Ntissue,Nphotons);
      pathlengthall=zeros(Ntissue,Nphotons);
      
    for jj=1:Ntissue
        %Delta_r2=6*D*tau;
        %DeltaR_cell=(1+DeltaV).^(1/3)-1;
        mua=tmua(jj);
        pathlength=his(jj+1,index);
        momentumTrasfer=his(jj+Ntissue+1,index);
        Nscattering=his(jj+2*Ntissue+1,index);
       MomentumTransferall(jj,:)=momentumTrasfer;
       Nscatteringall(jj,:)=Nscattering;
       pathlengthall(jj,:)=pathlength;
        
    end
 end
 
 %%
 Nscatteringplot=sum(Nscatteringall);
 pathlengthplot=sum(pathlengthall);
 %hist(pathlengthplot,50)
%  subplot(1,2,1)
%  hist(Nscatteringplot,50); xlabel('N_s'); ylabel('number of photons')
%  subplot(1,2,2)
%  hist(pathlengthplot,50); xlabel('L_n (mm)'); ylabel('number of photons')

 %%
 
 Nphotons=size(pathlengthall,2);
 
 centralWavelength= 800; %nm
  pathlengthall=pathlengthall(:,1:Nphotons);
  Nscatteringall=Nscatteringall(:,1:Nphotons);
 

%%
c=300; %mm/ns
n=1.33;
v=c/n;
 
%for l_coherence=100





%k_c=1/2/pi/l_coherence; % convert everything into mm^{-1}

k_c=pi/2/l_coherence; % Hui's formula

k_0=2*pi/(0.8*10^(-3)); % mm^{-1}

kall= k_0-3*k_c:k_c/10:k_0+3*k_c; % all the k vectors
amplitude=normpdf(kall,k_0,k_c); % intensity distribution each k value
amplitude=sqrt(amplitude/max(amplitude)); %normalize and transform transform it into amplitude

omegaall=kall*c;
deltaomega=omegaall(2)-omegaall(1);




tin_all_ns=xdtof_irf; % ns

%width_omega=k_c*c; % in mm, ns
%width_t=1/2/pi/width_omega;



%width_t=width_t*10; % we make the widht larger compared to the band limited withdth
%IRF_in=fftshift(fft(amplitude.^2)); % intensity distribution of the band limited IRF
IRF_in = ydtof_irf;
plot(tin_all_ns,IRF_in);
%%IRF_in=
IRF_in=IRF_in/sum(IRF_in); 
proball=rand(1,Nphotons);
IRF_prob=cumsum(IRF_in);
index2=zeros(1,Nphotons);


 

%  pathlengthall=pathlengthall(:,1:Nphotons);
%  Nscatteringall=Nscatteringall(:,1:Nphotons);
 for ii=1:Nphotons
    index2(ii)=find(proball(ii)<IRF_prob,1);
 end 
    
    
    
    t_incident=tin_all_ns(index2);
    t_L=sum(pathlengthall,1)/v;
    t_s=t_L+t_incident;
    
    t_radius=0.1; %widht of the window;
    index3=find(t_s<(t_s0+t_radius) & t_s>(t_s0-t_radius));
    Nphotons2=length(index3);
    
    pathlengthall=pathlengthall(:,index3);
    Nscatteringall=Nscatteringall(:,index3);
    
    Nsmax=max(Nscatteringall(:));
    mask=(zeros(Ntissue*Nphotons2,Nsmax));
    
 for ii=1:Nphotons2
    for jj=1:Ntissue
       mask((ii-1)*Ntissue+jj,1:Nscatteringall(jj,ii))=1; %% I still don't have a good way to assign the elements to the matrix
    end
end 
    
    index=find(mask~=0);
    clear mask
 
    
    

%%
initial_phase=rand(Nphotons2,1)*2*pi; % number of phases 
g=0; % we consider uniform scattering 
Dtissue=[D,D];% potentially we will have different diffusion coefficient in the two tissue types



%Eall_np_nk=zeros(Nphotons2,length(kall)); % number of photons * number of kvectors
Iall_t=zeros(1,length(tall));


xall =(zeros(Ntissue*Nphotons2,Nsmax));
yall =(zeros(Ntissue*Nphotons2,Nsmax));
zall =(zeros(Ntissue*Nphotons2,Nsmax));

xall2 =(zeros(Ntissue*Nphotons2,Nsmax));
yall2 =(zeros(Ntissue*Nphotons2,Nsmax));
zall2 =(zeros(Ntissue*Nphotons2,Nsmax));



qxall =(zeros(Ntissue*Nphotons2,Nsmax));
qyall =(zeros(Ntissue*Nphotons2,Nsmax));
qzall =(zeros(Ntissue*Nphotons2,Nsmax));

%phase_dynamics2=sparse(zeros)



    theta=acos((rand(1,length(index))-0.5)*2); %initialize the q matrix
    stheta=sin(theta);
    ctheta=cos(theta);
    phi=rand(1,length(index))*2*pi;
    sphi=sin(phi);
    cphi=cos(phi);
    ksx=k0.*stheta.*cphi; % the k vector of the scattered light
    ksy=k0.*stheta.*sphi;
    ksz=k0.*ctheta; 
    %ksx.^2+ksy.^2+ksz.^2 %sanity check
    qx=k0-ksx; % we assume the initial direction of the wave vector is in [1,0,0]
    qy=0-ksy;
    qz=0-ksz;
    qxall(index)=qx;
    qyall(index)=qy;
    qzall(index)=qz;
    
    phasetemp=pathlengthall(:)*kall;
    phasenow=zeros(Nphotons2,length(kall));
    
    for jj=1:Ntissue
       phasenow=phasenow+phasetemp(jj:Ntissue:end,:); % the accumulated phase for both tissue types  
    end
    
    clear phasetemp;
    
    
    pathlengthnow=repmat(sum(pathlengthall)',[1,length(kall)]);
    absorptionTerm=sqrt(exp(-mua*pathlengthnow));
    
    phasenow=phasenow+repmat(initial_phase,[1,length(kall)]); % add the initial phase the matrix is in nphotons*nk
    Eall_np_nk=exp(1i*phasenow).*absorptionTerm;
    Iall_t(1)=sum(abs(sum(Eall_np_nk,1)).^2,2);
    
    del=sqrt(2*Dtissue(1)*(tall(2)-tall(1)));% We first use the diffusion coefficient of the first medium and rescale it
    ratio=sqrt(Dtissue(2)/Dtissue(1)); 
    
% xall_sub = xall(index);
% yall_sub = yall(index);
% zall_sub = zall(index);
% qxall_sub = qxall(index);
% qyall_sub = qyall(index);
% qzall_sub = qzall(index);

n_index = length(index);

phase_dynamics_sum = zeros(Nphotons2,1);
phase_dynamics_sum_k = zeros(Nphotons2,length(kall));
amplitudeMatrix=repmat(amplitude,[Nphotons2,1]);



for tt=2:length(tall)
    if mod(tt,10)==0
       display(['Currenalty on tt=',num2str(tt)]) ;
    end
    xall(index) = xall(index)+del*randn(n_index,1); %Ntissue*Nphotons by Nscattering events
    yall(index) = yall(index)+del*randn(n_index,1); % pulling out the elements is expensive
    zall(index) = zall(index)+del*randn(n_index,1);
%     xall_sub=xall_sub+del*randn(n_index,1); %Ntissue*Nphotons by Nscattering events
%     yall_sub=yall_sub+del*randn(n_index,1); % pulling out the elements is expensive
%     zall_sub=zall_sub+del*randn(n_index,1);
%    
    xall2(1:2:end,:)=xall(1:2:end,:);
    xall2(2:2:end,:)=xall(2:2:end,:)*ratio;
    
    yall2(1:2:end,:)=yall(1:2:end,:);
    yall2(2:2:end,:)=yall(2:2:end,:)*ratio;
    
    zall2(1:2:end,:)=zall(1:2:end,:);
    zall2(2:2:end,:)=zall(2:2:end,:)*ratio;
    
        
    phase_dynamics=xall2.*qxall+yall2.*qyall+zall2.*qzall;
    
%     if Ntissue>1 % use this when the medium is not homogeneous
%     for jj=2:Ntissue
%         phase_dynamics(jj:Ntissue:end)=phase_dynamics(jj:Ntissue:end)*sqrt(Dtissue(jj)/Dtissue(1));
%     end
%     end

%    phase_dynamics_sum=sparse(zeros(Nphotons2,Nsmax));
    phase_dynamics_sum = sum(phase_dynamics(1:Ntissue:end,:),2);
    for jj=2:Ntissue % sum over the phase from both tissue types
        phase_dynamics_sum = phase_dynamics_sum + sum(phase_dynamics(jj:Ntissue:end,:),2);
    end
%    clear phase_dynamics % PRE ALLOCATE THIS VARIABLE BEFORE THE tt loop and DO NOT CLEAR
%    phase_dynamics_sum=sum(phase_dynamics_sum,2);
    
    phase_dynamics_sum_k = phase_dynamics_sum * kall/k0; % ALLOCATE NEW VARIABLE phase_dynamics_sum_k
    
    phase_dynamics_sum_k = phase_dynamics_sum_k + phasenow;
    
    Eall_np_nk = exp(1i*phase_dynamics_sum_k) .* amplitudeMatrix.*absorptionTerm;
    Iall_t(tt) = sum(abs(sum(Eall_np_nk,1)).^2,2);
      
    
end 

%mean(Iall_t.^2)/mean(Iall_t).^2;

save(['Iallt_multiWavelength_lc_',num2str(l_coherence),'ts_',num2str(t_s0),'_tw_',num2str(2*t_radius),'_config',num2str(config),'.mat'],'Iall_t','tall');

    end
    
    %end

%toc