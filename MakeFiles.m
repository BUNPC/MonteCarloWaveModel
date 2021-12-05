%% Inputfile Parameters


clc
clear

%cd('C:\Users\cxjph\Documents\CProgrammingUdemy\tMCimgTest_TDDCS\')

filenm='Tissue1';
xmax=80; %mm size in x
ymax=60; %mm
zmax=50; %mm

voxelsize=1; %mm
xvoxels=round(xmax/voxelsize);
yvoxels=round(ymax/voxelsize);
zvoxels=round(zmax/voxelsize);

Nphotons=10^8; % number of photons
randomSeed=1;

xi=25;%mm
yi=30;%mm
zi=voxelsize*2;%mm
NA=0;% we don't consider the numerical aperture
flagFocus=0; % we are not focusing our beam here
cxi=0; cyi=0; czi=1; %initial direction cosines of photons, now 

minT=0;maxT=5*10^(-9); stepT=maxT*2; %recording time span. Make stepT >maxT when not doing temporal measurements % get rid of the number of scattering events, 5 ns for the pathlength
nA1step=1; nA3step=1; % number of angular steps. Set to be 0 if not doing angular measurements;

xstep=voxelsize; nxstep=xvoxels; Ixmin=10; Ixmax=70; %region for recording fluence data
ystep=voxelsize; nystep=yvoxels; Iymin=28; Iymax=32;
zstep=voxelsize; nzstep=zvoxels; Izmin=1; Izmax=50;

Ntissue=2;       %two tissue types 
tmus=[1, 1]; % rho>ls* , ls*=1 mm
tg=[0.01 0.01]; % Note,decrease g to compare with the theory g=0.01,
tmua=[0.001 0.001];
tn=[1,1];
tProbCell=[0.1,0.1];% probability of a scattering event happens in a cell

nDets=7; %number of detectors
detRad=2; %radius of detectors N
xdetPos=xi+[10:5:40]; ydetPos=30+zeros(1,length(xdetPos));zdetPos=zi+zeros(1,length(xdetPos)); %detector size


%%
fileID = fopen([filenm,'.bin'],'w');

TissueProfile=ones(xvoxels,yvoxels,zvoxels);
TissueProfile(:,:,1)=0;% first layer to be zero

% x_local_center=40; % local region where the activation is
% y_local_center=30;
% z_local_center=30;
% xspan=5;
% yspan=5;
% zspan=2;
% 
% xrange=[round((x_local_center-xspan)/voxelsize):round((x_local_center+xspan)/voxelsize)];
% yrange=[round((y_local_center-yspan)/voxelsize):round((y_local_center+yspan)/voxelsize)];
% zrange=[round((z_local_center-zspan)/voxelsize):round((z_local_center+zspan)/voxelsize)];
zboundary=15; %mm
zboundaryv=round(zboundary/voxelsize);
TissueProfile(:,:,zboundaryv:end)=2;

% totalN=length(TissueProfile(:));
% for ii=1:totalN
%    fprintf(fileID,'%c',TissueProfile(ii)); 
%     
% end

fwrite(fileID,TissueProfile);

fclose(fileID);
%% 
fileID = fopen([filenm,'.inp'],'w');
fprintf(fileID,'%d\n',Nphotons);
fprintf(fileID,'%d\n',randomSeed);
%fprintf(fileID,'%d\n',randomSeed);
fprintf(fileID,'%f %f %f %f %d\n',xi,yi,zi,NA,flagFocus);
fprintf(fileID,'%f %f %f \n',cxi, cyi, czi);
fprintf(fileID,'%f %1.10f %1.10f \n',minT, maxT, stepT);

fprintf(fileID,'%d %d \n',nA1step, nA3step);
fprintf(fileID,'%s.bin \n',filenm);
%ASSERT(fscanf(fp, "%s", segFile)!=1);
fprintf(fileID,'%f %d %d %d \n',xstep, nxstep, Ixmin, Ixmax);
fprintf(fileID,'%f %d %d %d \n',ystep, nystep, Iymin, Iymax);
fprintf(fileID,'%f %d %d %d \n',zstep, nzstep, Izmin, Izmax);
fprintf(fileID,'%d\n',Ntissue);
for ii=1:Ntissue
fprintf(fileID,'%f %f %f %f %f\n',tmus(ii),tg(ii),tmua(ii),tn(ii),tProbCell(ii));
end
fprintf(fileID,'%d %f\n',nDets,detRad);
for ii=1:nDets
    fprintf(fileID,'%f %f %f \n',xdetPos(ii),ydetPos(ii),zdetPos(ii));
end



%"%d %lf", &nDets, &detRad
%detPos[i], detPos[i]+1, detPos[i]+2

fclose(fileID);

%% Run the c program with executable
%clc
%clear
exefile='tMCimgTest';
inputfile='Tissue1';
system([exefile,' ',inputfile]); 


%% record the fluence data
filenm='Tissue1';

fid = fopen([filenm,'.2pt'],'rb');
Io = fread(fid,inf,'float32');
fclose(fid);

Io=reshape(Io,[(Ixmax-Ixmin)+1,(Iymax-Iymin)+1,(Izmax-Izmin)+1]);
imagesc(squeeze(Io(:,3,:)));
% %% read history file
% clc
% clear
% nDets=7;
% Ntissue=2;
% filenm='Tissue1';
% fid = fopen([filenm,'.his'], 'rb');
%     his = fread(fid, inf, 'float');
% totalP=length(his);
% VariableNumber=3;
% his=reshape(his,[1+Ntissue*VariableNumber,totalP/(1+Ntissue*VariableNumber)]);
% % calculate g1 and g2
% 
% 
% 
% tau=10.^[-4:0.1:0];%:(0:20)*10; %ms
% Ntau=length(tau);
% D=1*10^(-7);
% 
% k0=2*pi/(1.2*10^(-3));
% tmua=[0,0];
%  for ii=1:nDets
% 
%      index=find(his(1,:)==ii-1);
%      Nphotons=length(index);
%      NphotonsAll(ii)=Nphotons;
%      g1=ones(Nphotons,Ntau);
%     for jj=1
%         Delta_r2=6*D*tau;
%         
%         mua=tmua(jj);
%         pathlength=his(jj+1,index);
%         momentumTrasfer=his(jj+Ntissue+1,index);
%       Delta_r2=repmat(Delta_r2,[Nphotons,1]);
%       momentumTrasfer=squeeze(repmat(momentumTrasfer',[1,Ntau]));
%       pathlength=squeeze(repmat(pathlength',[1,Ntau]));
%       g1partial=(exp(-1/3.*momentumTrasfer.*Delta_r2*k0^2).*exp(-pathlength*mua));% multiple before taking the mean
%       g1=g1.*g1partial;
% %Pathlength_photon_tissuetype(ii,jj,:)=his((jj+(ii-1)*Ntissue*VariableNumber):Ntissue*nDets*VariableNumber:end);
% %MomentumTransfer_photon_tissuetype(ii,jj,:)=his((jj+(ii-1)*Ntissue*VariableNumber+Ntissue):Ntissue*nDets*VariableNumber:end);
%     end
%     g1=mean(g1); 
%     f=polyfit(sqrt(tau(5:end)),log(g1(5:end)),1);
%     slopeall(ii)=f(1);
%     
%     plot(sqrt(tau),log(g1));hold on
% 
%  end
% legend('10 mm','15 mm', '20 mm', '25 mm', '30 mm','35 mm', '40 mm')
% xlabel('\tau (ms)');
% ylabel('ln(g1)')
