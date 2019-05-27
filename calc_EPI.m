
function [K_ADC,G_ADC]=calc_EPI(ICE)
    
	% Calculate the k-space path
    dKy=linspace(-ICE.Kmax,ICE.Kmax, ICE.Matrix);
    dKx=linspace(-ICE.Kmax,ICE.Kmax, ICE.NColMeas);
    
    [X,Y] = meshgrid(dKy,dKx);
    
    K_ADC(:,:,1)=Y;
    K_ADC(:,:,2)=X;
    G_ADC(:,:,1)=diff(K_ADC(:,:,1));
    G_ADC(:,:,2)=diff(K_ADC(:,:,2));
   
end