close all;
clear all;

file=uigetfile('.dat');
twix_obj = mapVBVD(file);

%% %%%%%%%%%%%%%%  EXTRACT RAW DATA INFORMATION %%%%%%%%%%%%%%

ICE=[];
ICE.DwellTime=twix_obj{1, 2}.hdr.MeasYaps.sRXSPEC.alDwellTime{1}/1000; % us
ICE.GradRasterTime=10; % us
ICE.dT=1e-6; % s
ICE.Gamma=42575.6; % Hz/mT from 42.576 MHz/T
ICE.dTwoPiPrecis = 6.2831853;

ICE.FOV=twix_obj{1, 2}.hdr.Config.RoFOV/1000; % mm -> m
ICE.NLine= twix_obj{1, 2}.hdr.Config.ImageLines;
ICE.NCol=twix_obj{1, 2}.hdr.Config.ImageColumns;%twix_obj{1, 2}.hdr.Config.RawCol;
ICE.NColMeas=ICE.NCol*2;%twix_obj{1, 2}.hdr.Config.NColMeas;
ICE.NLineMeas=twix_obj{1, 2}.hdr.Config.NLinMeas;

ICE.EchoPosPhase=twix_obj{1, 2}.image.centerLin-1;
ICE.EchoPosRead=twix_obj{1, 2}.image.centerCol-1;

ICE.LineIdx=ICE.EchoPosPhase(1)+twix_obj{1, 2}.image.Lin(1:ICE.NLineMeas);
ICE.Matrix=twix_obj{1, 2}.hdr.Config.BaseResolution;
ICE.dK=1/ICE.FOV;
ICE.Kmax=ICE.Matrix/(2*ICE.FOV);	
ICE.Nbleaves=twix_obj{1, 2}.hdr.Config.NSlc;

%% %%%%%%%%%%%%%% GENERATE K-SPACE TRAJECTORY FROM PROTOCOL %%%%%%%%%%%%%%
[K_ADC]=calc_EPI(ICE); % K_ADC: K-space trajectory "interpolated" to dwell Times, G_ADC: Gradients interpolated to dwell Times

%% %%%%%%%%%%%%%% RECON PARAMETER %%%%%%%%%%%%%%

RECON=[];
RECON.MatrixOverSampling=1;
RECON.Delay=[1,1];  % Gradient delay units in points
RECON.Kshift=[];
RECON.Matrix=ICE.Matrix;
RECON.Debug=1; % 1 for display the plot
RECON.Nbleaves=ICE.Nbleaves;

%% %%%%%%%%%%%%%% Reorder the RAW data %%%%%%%%%%%%%%

RAW=[];
RAW_tmp=[];

RAW_EarlyFid=[];

RAW_Prescan=[];
RAW_Prescan_tmp=[];

RAW_tmp=squeeze(twix_obj{1,2}.image(:,:,:,1,:,1,1,1,1,1,:,1));              % RAW Columns Chanels Lines Slices Segm IDA

% Combine Odd and Even echoes
RAW=zeros(ICE.NColMeas,size(RAW_tmp,2),ICE.Matrix,ICE.Nbleaves);
RAW_test=zeros(ICE.NColMeas,size(RAW_tmp,2),ICE.Matrix,ICE.Nbleaves);
RAW(:,:,ICE.LineIdx,:)=squeeze(RAW_tmp(:,:,:,:,1)+squeeze(RAW_tmp(:,:,:,:,2))); 


% Early FID data (three times the central K-space Line [+ - +])
RAW_EarlyFid(:,:,:,1)=squeeze(twix_obj{1,2}.phasecor(:,:,twix_obj{1, 2}.phasecor.NLin,1,:,1,1,1,1,1,1)); % RAW Prescan Columns Chanels Lines Slices
RAW_EarlyFid(:,:,:,2)=squeeze(twix_obj{1,2}.phasecor(:,:,twix_obj{1, 2}.phasecor.NLin,1,:,1,1,1,1,1,2)); % RAW Prescan Columns Chanels Lines Slices
RAW_EarlyFid(:,:,:,3)=squeeze(twix_obj{1,2}.phasecor(:,:,twix_obj{1, 2}.phasecor.NLin,1,:,2,1,1,1,1,1)); % RAW Prescan Columns Chanels Lines Slices

RAW_Prescan_tmp=twix_obj{1,2}.image{''};
RAW_Prescan=zeros(ICE.NColMeas,size(RAW_tmp,2),ICE.Matrix,ICE.Nbleaves,4);
RAW_Prescan(:,:,ICE.LineIdx,:,1)=squeeze(RAW_Prescan_tmp(:,:,:,:,1,1,2))+squeeze(RAW_Prescan_tmp(:,:,:,:,1,2,2)); % RAW Columns Chanels Lines Slices Segm Rep IDA
RAW_Prescan(:,:,ICE.LineIdx,:,2)=squeeze(RAW_Prescan_tmp(:,:,:,:,2,1,3))+squeeze(RAW_Prescan_tmp(:,:,:,:,2,2,3));
RAW_Prescan(:,:,ICE.LineIdx,:,3)=squeeze(RAW_Prescan_tmp(:,:,:,:,3,1,4))+squeeze(RAW_Prescan_tmp(:,:,:,:,3,2,4));
RAW_Prescan(:,:,ICE.LineIdx,:,4)=squeeze(RAW_Prescan_tmp(:,:,:,:,4,1,5))+squeeze(RAW_Prescan_tmp(:,:,:,:,4,2,5));

RAW=permute(RAW,[1 3 2 4]);                      % RAW Columns Lines Chanels Shots
RAW_EarlyFid=permute(RAW_EarlyFid,[1 3 2 4 5]);  % RAW Columns Lines Chanels Shots Scan
RAW_Prescan=permute(RAW_Prescan,[1 3 2 4 5]);    % RAW Columns Lines Chanels Shots Scan

clear RAW_tmp RAW_Prescan_tmp

%% %%%%%%%%%%%%%% (WIP) ESTIMATE EPI TRAJECTORY CORRECTION FROM PRESCAN AND APPLY IT (RAW -> RAW2) %%%%%%%%%%%%%%

for cpt=1:1:RECON.Nbleaves
    [RECON.Kshift(cpt,:)] = correction_traj_EPI2(squeeze(RAW_Prescan(:,:,:,cpt,2)),K_ADC(:,:,1),ICE); 
    [RAW2(:,:,:,cpt)] = correction_traj_shift(squeeze(RAW(:,:,:,cpt)),RECON.Kshift(cpt,:))  ;
end
RAW2=RAW; % As it's still WIP we do not apply it for now.

%% %%%%%%%%%%%%%% GENERATE K-SPACE trajectory for every shot (From K_ADC -> K_rot) %%%%%%%%%%%%%%

for cpt=1:1:RECON.Nbleaves
    ICE.dRotAngle(cpt) = ICE.dTwoPiPrecis *(cpt-1)/(ICE.Nbleaves); % Fix Angle version
    ICE.RotMat(:,:,cpt) = [cos(ICE.dRotAngle(cpt)) -sin(ICE.dRotAngle(cpt));sin(ICE.dRotAngle(cpt)) cos(ICE.dRotAngle(cpt))];
    for cpt2=1:1:size(K_ADC,1)
        for cpt3=1:1:size(K_ADC,2)
            K_rot(cpt2,cpt3,cpt,:)=squeeze(K_ADC(cpt2,cpt3,:))'*squeeze(ICE.RotMat(:,:,cpt));
        end
    end
end

%% %%%%%%%%%%%%%% (WIP) PARTIAL FOURRIER POCS (From RAW2 -> RAWPocs) %%%%%%%%%%%%%%

RAWPocs=[];
% RAW_tmp=permute(RAW,[3 1 2 4]); 
%  for cpt_shot=1:1:size(RAW_tmp,4)
%      [~, RAWPocs(:,:,:,cpt_shot)]= pocs( RAW_tmp(:,:,:,cpt_shot), 30, true );
%  end
% RAWPocs= permute(RAWPocs,[2 3 1 4]);
RAWPocs=RAW2; % As it's still WIP we do not apply it for now.

%% %%%%%%%%%%%%%% PHASE CORRECTION (From RAWPocs -> RAW3) %%%%%%%%%%%%%% 

[RAW3]= phase_correction(RAWPocs);

%% %%%%%%%%%%%%%%  (WIP) KSPACE CORRECTION of ROTATION AND TRANSLATION (From RAW3 & K_rot -> RAW4 & K_rot2 ) %%%%%%%%%%%%%% 

[RAW4 K_rot2] = correction_rot_trans(RAW3,K_rot);  
K_rot2=K_rot; % As it's still WIP we do not apply it for now.
RAW4=RAW3;    % As it's still WIP we do not apply it for now.

%% %%%%%%%%%%%%%% NORMALIZE K-SPACE trajectory & FINAL RAW DATA (From K_rot2 & RAW4 -> K1 & RAW_final) %%%%%%%%%%%%%%

K1=complex(squeeze(K_rot2(:,:,:,1)./(2*max(max(max(abs(K_rot2(:,:,:,1))))))),squeeze(K_rot2(:,:,:,2)./(2*max(max(max(abs(K_rot2(:,:,:,2))))))));
K1=reshape(K1,size(K1,1)*size(K1,2),size(K1,3));
RAW_final=reshape(RAW4,size(RAW4,1)*size(RAW4,2),size(RAW4,3),size(RAW4,4));

RECON.Dcf_s=ones(size(K1,1),1); % Todo: Density function calculation 

%% %%%%%%%%%%%%%% REGRIDDING IMG (From K1 & RAW_final -> Final_Img & Final_Phase) %%%%%%%%%%%%%%

RAW_tmp=[];
K_tmp=[];
RECON.Dcf=[];
Final_Img=[];
Final_Phase=[];

% reorder the multi-shot data into one format
for cpt_lv=1:1:RECON.Nbleaves
    RAW_tmp=[RAW_tmp ;squeeze(RAW_final(:,:,cpt_lv))];    
    RECON.Dcf=[RECON.Dcf ;RECON.Dcf_s];
    K_tmp=[K_tmp ;K1(:,cpt_lv)];
end

% Regridding and Coil combination
disp('Recon Img');
[Final_Img(:,:) S]= nufft_n_combine(RAW_tmp,K_tmp,RECON,RECON.Matrix,2);
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


function [Img, S] = grid_n_combine(RAW,K,RECON,Matrix,CoilMode)
    
    Grid_K=[];  
    S=[];
     
     %% Regridding 
     disp('Coil');
     h = waitbar(0,'Coil...');
    for cpt_coil=1:1:size(RAW,2)
        Data=RAW(:,cpt_coil); % Take the acquisition part after grad delay;
        Data=Data(:);
        K=K(:);
        
       [tmp_dat m] = gridmat(K,Data,RECON.Dcf,Matrix);
        Grid_K(:,:,cpt_coil)=tmp_dat;
        waitbar(cpt_coil/size(RAW,2),h);
    end   
    close(h)
    
    %% FFT (Grid_K -> Grid_Img)
    Grid_Img=fftshift(fft2(fftshift(Grid_K)));
    Grid_Img(isnan(Grid_Img))=0;
    
    %% Deconvolution
    Grid_Img=Grid_Img./m; 
        
    %% % COIL COMBINATION (From Grid_Img -> Final_Img): Sum of square of the magnitude
    if (CoilMode==2) % Adaptive coil combine (eSPIRIT), works well with a pre-coil selection (WIP)
        [val index]=maxk(squeeze(median(median(abs(Grid_Img(:,:,:))))),12);
        S = adaptive_coil(Grid_Img(:,:,index));
        Img=sum(Grid_Img(:,:,index).*conj(S),3)./sum(S.*conj(S),3);
    else% Sum of square
        Img=sqrt(sum((Grid_Img).^2,3));   
    end   
    
end

function [Img, S] = nufft_n_combine(RAW,K,RECON,Matrix,CoilMode)  
     
     Grid_Img=[];
     S=[];
     
     %% Nufft
     st = nufft_init([2*pi*real(K) 2*pi*imag(K)] , [Matrix Matrix], [5 5], [Matrix*2 Matrix*2],[Matrix/2 Matrix/2], 'minmax:kb');  % Basic parameters
     disp('Coil');
     h = waitbar(0,'Coil...');
    for cpt_coil=1:1:size(RAW,2)
        Data=complex(double(real(RAW(:,cpt_coil))), double(imag(RAW(:,cpt_coil)))).*RECON.Dcf ; % Take the acquisition part after grad delay;    
        Grid_Img(:,:,cpt_coil)  = nufft_adj(Data,st); 
        waitbar(cpt_coil/size(RAW,2),h);
    end   
    close(h)
    
    %% % COIL COMBINATION (From Grid_Img -> Final_Img): Sum of square of the magnitude
    if (CoilMode==2) % Adaptive coil combine (eSPIRIT), works well with a pre-coil selection (WIP)
        [val index]=maxk(squeeze(median(median(abs(Grid_Img(:,:,:))))),12);
        S = adaptive_coil(Grid_Img(:,:,index));
        Img=sum(Grid_Img(:,:,index).*conj(S),3)./sum(S.*conj(S),3);
    else% Sum of square
        Img=sqrt(sum((Grid_Img).^2,3));   
    end
end

function [RAW2]= phase_correction(RAW)
    w=window2_KM(size(RAW,1),size(RAW,2),@triang,1/4);
    RAW2=[];
    RAW_tmp=[];
    for cpt_shot=1:1:size(RAW,4)   
         for cpt_coil=1:1:size(RAW,3)       
            RAW_tmp=squeeze(RAW(:,:,cpt_coil,cpt_shot)).*w;

            IMG_triangle=(fft2(ifftshift(RAW_tmp)));

            IMG=(fft2(ifftshift(squeeze(RAW(:,:,cpt_coil,cpt_shot)))));

            IMG2=abs(IMG).*exp(i*( angle(IMG)-angle(IMG_triangle)));

            RAW2(:,:,cpt_coil,cpt_shot)=fftshift(ifft2((IMG2)));
         end 
    end   
end

function [RAW2] = correction_traj_shift(RAW,Kshift)    
    %RAW2=RAW;
    vectX=1:size(RAW,1);
        for cpt=1:1:size(RAW,2) 
                    R_tmp=RAW(:,cpt,:);           
                    RAW2(:,cpt,:)=ifft( ( fft(R_tmp).* repmat( exp(1i*pi*Kshift(cpt)*0.0005.*vectX'),1,1,size(RAW,3) ) ) );
        end
end

function [RAW2 K2] = correction_rot_trans(RAW,K)    
    RAW2=RAW;
    RAWtmp=squeeze(max(RAW(1:2:end,:,:,:),[],3));
    K2=K;
    for cpt=1:1:size(RAWtmp,3)
        dRotAngle(cpt) = 2*pi *(cpt-1)/(size(RAWtmp,3)); % Fix Angle version
        RAWtmp2(:,:,cpt) = imrotate(abs(RAWtmp(:,:,cpt)),180*dRotAngle(cpt)/pi,'bilinear','crop');
    end
    L=floor(65/sqrt(2)/2);
    RAW4=squeeze((RAWtmp2(end/2-L:1:end/2+L,end/2-L:end/2+L,:,:)));
     [optimizer, metric] = imregconfig('multimodal');
    RAW5=RAW4; 
    

    Mt(:,:,1)=[1 0 0; 0 1 0; 0 0 1];
    for cpt=2:1:size(RAW4,3)
        % Rz = [cos(t) -sin(t) 0; sin(t) cos(t) 0; 0 0 1]
        %tform=affine2d([cos(dRotAngle(cpt)) -sin(dRotAngle(cpt)) 0;sin(dRotAngle(cpt)) cos(dRotAngle(cpt)) 0 ;0 0 1]);  
       [RAW5(:,:,cpt) tt]=imregister(abs(RAW4(:,:,cpt)), abs(RAW4(:,:,1)), 'rigid', optimizer, metric);
       tt = imregtform(abs(RAW4(:,:,cpt)), abs(RAW4(:,:,1)), 'rigid', optimizer, metric);
       Mt(:,:,cpt)=tt.T;
    end
    Mt(:,:,cpt)
    for cpt2=1:1:size(K,1)
        for cpt3=1:1:size(K,2)
            for cpt=1:1:size(RAW4,3)
                K2(cpt2,cpt3,cpt,:)=squeeze(K(cpt2,cpt3,cpt,:))'*squeeze(Mt(1:2,1:2,cpt));
            end
        end
    end
%     tform = affine2d([2 0.33 0; 0 1 0; 0 0 1])
%     'InitialTransformation'
end

function FFT_n_Display(RAW)
    for cpt_shot=1:1:size(RAW,3)
        for cpt_coil=1:1:size(RAW,4)
            ImgCar(:,:,cpt_shot,cpt_coil)=fftshift(ifft2((RAW(:,:,cpt_shot,cpt_coil))));       
        end
    end    
     ImgCar=sqrt(sum(ImgCar.^2,3));
     figure, imagesc(abs(squeeze(ImgCar(:,:,:)))) 
end
