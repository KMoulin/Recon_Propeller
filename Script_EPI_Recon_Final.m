close all;
clear all;

file=uigetfile('.dat');
twix_obj = mapVBVD(file);

%% EXTRACT RAW DATA INFORMATION
ICE=[];
ICE.DwellTime=twix_obj{1, 2}.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1000; % us
ICE.GradRasterTime=10; % us
ICE.dT=1e-6; % s
ICE.Gamma=42575.6; % Hz/mT from 42.576 MHz/T
ICE.dTwoPiPrecis = 6.2831853;

ICE.FOV=twix_obj{1, 2}.hdr.Config.RoFOV/1000; % mm -> m
ICE.NLine=twix_obj{1, 2}.hdr.Config.RawLin;
ICE.NCol=twix_obj{1, 2}.hdr.Config.RawCol;
ICE.NColMeas=twix_obj{1, 2}.hdr.Config.NColMeas;
ICE.EchoPosReadout=twix_obj{1, 2}.image.centerCol;
ICE.EchoPosPhase=twix_obj{1, 2}.image.centerLin;
ICE.Matrix=twix_obj{1, 2}.hdr.Config.BaseResolution;

ICE.dK=1/ICE.FOV;
ICE.Kmax=ICE.Matrix/(2*ICE.FOV);	

ICE.Nbleaves=twix_obj{1, 2}.hdr.Config.RawSlc;

%% GENERATE K-SPACE TRAJECTORY FROM PROTOCOL
[K_ADC]=calc_EPI(ICE); % K_ADC: K-space trajectory "interpolated" to dwell Times, G_ADC: Gradients interpolated to dwell Times

%% RECON PARAMETER
RECON=[];
RECON.MatrixOverSampling=1;
RECON.Delay=[1,1];  % Gradient delay units in points
RECON.Kshift=[0,0];
RECON.Matrix=ICE.Matrix;
RECON.Debug=1; % 1 for display the plot
RECON.Nbleaves=ICE.Nbleaves;

%% Reorder the RAW data
RAW=[];
RAW_PreScan=[];
RAW_tmp=[];
RAW_tmp=squeeze(twix_obj{1,2}.image(:,:,:,1,:,1,1,1,1,1,:));              % RAW Columns Chanels Lines Slices Segm
RAW=zeros(ICE.NColMeas,size(RAW_tmp,2),ICE.Matrix,ICE.Nbleaves);

% Combine Odd and Even echoes
RAW(:,:,1:size(RAW_tmp,3),:)=squeeze(RAW_tmp(:,:,:,:,1)); 
RAW(:,:,2:2:size(RAW_tmp,3),:)=squeeze(RAW_tmp(:,:,2:2:end,:,2));

% Early FID data (three times the central K-space Line [+ - +])
RAW_PreScan(:,:,:,1)=squeeze(twix_obj{1,2}.phasecor(:,:,twix_obj{1, 2}.phasecor.NLin,1,:,1,1,1,1,1,1)); % RAW Prescan Columns Chanels Lines Slices
RAW_PreScan(:,:,:,2)=squeeze(twix_obj{1,2}.phasecor(:,:,twix_obj{1, 2}.phasecor.NLin,1,:,1,1,1,1,1,2)); % RAW Prescan Columns Chanels Lines Slices
RAW_PreScan(:,:,:,3)=squeeze(twix_obj{1,2}.phasecor(:,:,twix_obj{1, 2}.phasecor.NLin,1,:,2,1,1,1,1,1)); % RAW Prescan Columns Chanels Lines Slices


%% EARLY FID & EPI Trajectory corretion 
% [RECON.Kshift(1)] = correction_traj_EPI2(RAW_PreScan,K_ADC(:,ICE.EchoPosPhase,1),ICE); 


%% GENERATE K-SPACE trajectory for every shot (From K_ADC -> K_rot)
for cpt=1:1:RECON.Nbleaves
    ICE.dRotAngle(cpt) = ICE.dTwoPiPrecis *(cpt-1)/(ICE.Nbleaves/2); % Fix Angle version
    ICE.RotMat(:,:,cpt) = [cos(ICE.dRotAngle(cpt)) -sin(ICE.dRotAngle(cpt));sin(ICE.dRotAngle(cpt)) cos(ICE.dRotAngle(cpt))];
    for cpt2=1:1:size(K_ADC,1)
        for cpt3=1:1:size(K_ADC,2)
            K_rot(cpt2,cpt3,cpt,:)=squeeze(K_ADC(cpt2,cpt3,:))'*squeeze(ICE.RotMat(:,:,cpt));
        end
    end
end

%% NORMALIZE K-SPACE trajectory (From K_rot -> K1)

K1=complex(squeeze(K_rot(:,:,:,1)./(2*max(max(max(abs(K_rot(:,:,:,1))))))),squeeze(K_rot(:,:,:,2)./(2*max(max(max(abs(K_rot(:,:,:,2))))))));
K1=reshape(K1,size(K1,1)*size(K1,2),size(K1,3));

RECON.Dcf_s=ones(size(K1,1),1); % Todo: Density function calculation 


%% REORDER the RAW data (From RAW -> RAW4)
RAW1=permute(RAW,[1 3 4 2 5]);
for cpt_shot=1:1:size(RAW1,3)
    for cpt_coil=1:1:size(RAW1,4)
        ImgCar(:,:,cpt_shot,cpt_coil)=abs(fftshift(fft2(RAW1(:,:,cpt_shot,cpt_coil))));       
    end
end    
ImgCar=sqrt(sum(ImgCar.^2,4));

figure, imagescn(abs(ImgCar))

RAW1=reshape(RAW1,size(RAW1,1)*size(RAW1,2),size(RAW1,3),size(RAW1,4));

%% REGRIDDING IMG (From K1 & RAW1 -> Final_Img & Final_Phase)
RAW_tmp=[];
K2=[];
RECON.Dcf=[];
Final_Img=[];
Final_Phase=[];

% reorder the multi-shot data into one format
for cpt_lv=1:1:RECON.Nbleaves
    RAW_tmp=[RAW_tmp ;squeeze(RAW1(:,cpt_lv,:))];    
    RECON.Dcf=[RECON.Dcf ;RECON.Dcf_s];
    K2=[K2 ;K1(:,cpt_lv)];
end

% Regridding and Coil combination
disp('Recon Img');
[Final_Img(:,:) S M]= grid_n_combine(RAW_tmp,K2,RECON,RECON.Matrix);
Final_Img(isnan(Final_Img))=0;
Final_Phase(:,:)=phase_unwrap(angle(Final_Img(:,:)));
     
clear cpt cpt_lv;

%% Final Image
figure
imagesc(abs(Final_Img(:,:,1))),colormap('gray');
figure
imagesc(angle(Final_Img(:,:,1)));



%% ESPIRIT and ADAPTIVE COIL COMBINE 
function [S M]= adaptive_coil(data)
    data_2d(:,:,1,:)=data;
    [Nx,Ny,Nz,Nc] = size(data_2d);
    S = zeros(Nx,Ny,Nz,Nc);
    M = zeros(Nx,Ny,Nz);
    w = 5;
    for i = 1:Nx
        ii = max(i-w,1):min(i+w,Nx);
        for j = 1:Ny
            jj = max(j-w,1):min(j+w,Ny);
            for k = 1:Nz
                kk = max(k-w,1):min(k+w,Nz);
                kernel = reshape(data_2d(ii,jj,kk,:),[],Nc);
                [V,D] = eigs(conj(kernel'*kernel),1);
                S(i,j,k,:) = V*exp(-1j*angle(V(1)));
                M(i,j,k) = sqrt(D);
            end
        end
    end
    S = squeeze(S.*(M>0.01*max(abs(M(:)))));
end

function [Img, S, M] = grid_n_combine(RAW,K,RECON,Matrix)
    
    Grid_K=[];  

     disp('Coil');
     h = waitbar(0,'Coil...');
    for cpt_coil=1:1:size(RAW,2)
        Data=RAW(:,cpt_coil); % Take the acquisition part after grad delay;
        Data=Data(:);
        K=K(:);  
       [dat m] = gridmat(K,Data,RECON.Dcf,Matrix);
        
        Grid_K(:,:,cpt_coil)=dat;
        waitbar(cpt_coil/size(RAW,2),h);
    end   
    close(h)
    
    % FFT (Grid_K -> Grid_Img)
    Grid_Img=fftshift(fft2(fftshift(Grid_K)));
    Grid_Img(isnan(Grid_Img))=0;
    
    % Deconvolution
    Grid_Img=Grid_Img./m; 

    % ADAPTIE COIL COMBINATION
    [S M]= adaptive_coil(Grid_Img);
    
    % COIL COMBINATION (From Coil_Img -> Final_Img): Sum of square of the magnitude
    %S_test(:,:,:,cpt_lv,cpt_diff)=S;
    %Img=sqrt(sum((Grid_Img).^2,3));
    Img=sum(Grid_Img.*conj(S),3)./sum(S.*conj(S),3);
    
end



