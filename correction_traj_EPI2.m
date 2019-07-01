function [dKshift] = correction_traj_EPI2(RAW,K,ICE)
 

      % RAW2=RAW;
      RAW_tmp=max(RAW,[],3); % Take the max of each channel  
      [Echo Echoi]= max(abs(RAW_tmp(:,:)));
   
     % figure,plot(K(:,1),abs(RAW_tmp))
      
      dKshift=ICE.EchoPosRead(1)-Echoi;
      [RAW2] = correction_traj_shift_local(RAW_tmp,dKshift)  ;
      [Echo2 Echoi2]= max(abs(RAW2(:,:)));
     % figure,plot(K(:,1),abs(RAW2))
%       for cpt_coil=1:1:size(RAW,2)
%             RAW2(:,cpt_coil,cpt_line,1)=RAW;
%             for cpt_line=1:1:size(K,3)   % For every Line. 
%                 dK=0.001*cpt_dK.*K(:,cpt_line);
%                 sign=(-1)^cpt_line;
%                 R_tmp=RAW(:,cpt_coil,cpt_line,1);           
%                 RAW2(:,cpt_coil,cpt_line,1)=ifft(((fft(R_tmp)).*exp(sign*1i*2*pi*dK)));
%             end
%        end
    
end

function [RAW2] = correction_traj_shift_local(RAW,Kshift)    
    %RAW2=RAW;
    vectX=1:size(RAW,1);
    for cpt=1:1:size(RAW,2) 
                R_tmp=RAW(:,cpt);           
                RAW2(:,cpt)=ifft( ( fft(R_tmp).*exp(1i*pi*Kshift(cpt)*0.0005.*vectX') ) );
    end
end