%
%	Function [dat] = gridmat(ksp,kdat,dcf,gridsize)
%
%	function performs gridding in Matlab.  This is designed to
%	be a teaching function, and is not necessarily fast(!), but
%	reasonable for 2D gridding.
%
%	INPUT:
%         ksp = k-space locations, kx+i*ky, normalized to |k|<=0.5, 1/pixel.
%         dat = complex data samples at k-space locations.
%	  dcf = density compensation factors at k-space locations.
%	  gridsize = size of grid (default is 256)
%
%	OUTPUT:
%	  dat = matrix of gridded data samples
%
%	Many parameters come from Jackson 1991 here.

%	B.Hargreaves - October 2014.

function [dat m] = gridmat(KSpace,KData,dcf,MatrixSize)

KSpace=KSpace(:);		% Translate to 1D vector in case this is 2D
dcf=dcf(:);
KData=KData(:);


% -- Use a lookup-table for grid kernel
Kernel.Width = 3;		% Convolution kernel width
Kernel.Beta = 4.2;		% Kaiser-bessel Beta
Kernel.Dr = 0.01;			 
Kernel.R = [0:Kernel.Dr:2*Kernel.Width];		% input to kaiser-bessel
Kernel.Tab = kb(Kernel.R,Kernel.Width,Kernel.Beta);	% kaiser-bessel function
%Kernel.Tab = kb2(Kernel.R,Kernel.Width,Kernel.Beta);	% kaiser-bessel function


% -- Allocate grid to accumulate data
PadMatrixSize = MatrixSize+4*Kernel.Width;			% Padded grid to avoid errors
PadMatrix = zeros(PadMatrixSize,PadMatrixSize);	% Padded grid.

% -- Sample Density correction
KData = KData .* dcf;		% Density correct (Simple!).



% -- Scale kspace points to grid units
kxi = real(KSpace)*(MatrixSize-1)+PadMatrixSize/2;	% Scale kx to grid units
kyi = imag(KSpace)*(MatrixSize-1)+PadMatrixSize/2;	% Scale ky to grid units

% -- Round to nearest grid point
ikxi = round(kxi);	% Closest integer
ikyi = round(kyi);	% Closest integer

% -- Make a small matrix, that completely encompasses the grid points
% 	within the convolution kernel around a k-space sample --

sgridext = ceil(Kernel.Width/2)+1;	% Size of small grid around sample
[smaty,smatx] = meshgrid([-sgridext :sgridext ],[-sgridext :sgridext ]);
sgrid = 0*smatx;		% allocate


% -- Go through k-space samples to do convolution

for p = 1:length(KSpace)
  

  gridx = smatx + ikxi(p);	% grid of 'closest' integer pts
  gridy = smaty + ikyi(p);	% same in y	
			
  % -- Distance index (array), to use for kernel lookup	 
  %    Just calculating kernel value between ksp(p) and every point
  %    in this small grid.
  dist = round(sqrt((gridx - kxi(p)).^2 + (gridy -kyi(p)).^2)/Kernel.Dr)+1;

  sgrid(:) = KData(p) * Kernel.Tab(dist(:)); 	% Convolve sample w/ kernel

  % -- Add the 'small-grid' into the padded grid
  PadMatrix(ikxi(p)-sgridext:ikxi(p)+sgridext,ikyi(p)- sgridext:ikyi(p)+sgridext) = PadMatrix(ikxi(p)-sgridext:ikxi(p)+sgridext,ikyi(p)-sgridext:ikyi(p)+sgridext) + sgrid;
   
end;

% -- Extract the main grid from the padded grid.
dat = PadMatrix(2*Kernel.Width+1:2*Kernel.Width+MatrixSize, 2*Kernel.Width+1:2*Kernel.Width+MatrixSize);
m = kb_ft(MatrixSize,MatrixSize,Kernel.Width,Kernel.Beta);

end



% Kaiser-bessel Kernel


%
%	function y = kb(u,w,beta)
%
%	Computes the Kaiser-Bessel function used for gridding, namely
%
%	y = f(u,w,beta) = I0 [ beta*sqrt(1-(2u/w)^2) ]/w
%
%	where I0 is the zero-order modified Bessel function
%		of the first kind.
%
%	INPUT:
%		u = vector of k-space locations for calculation.
%		w = width parameter - see Jackson et al.
%		beta = beta parameter - see Jackson et al.
%
%	OUTPUT:
%		y = vector of Kaiser-Bessel values.
%
%	SEE ALSO:
%		kbk2x.m

%	B. Hargreaves	Oct, 2003.




function y = kb(u,w,beta)

if (nargin < 3)
	error('Not enough arguments -- 3 arguments required. ');
end;


if (length(w) > 1)
	error('w should be a single scalar value.');
end;


y = 0*u;				% Allocate space.
uz = find(abs(u)< w/2);			% Indices where u<w/2.

if (length(uz) > 0)				% Calculate y at indices uz.
	x = beta*sqrt(1-(2*u(uz)/w).^2);	% Argument - see Jackson '91.
	y(uz) = besseli(0,x)./w;
end;

y = real(y);		% Force to be real.

% y fft ->

end

function y = kb2(u,w,beta)
% function y = kb(u,w,beta)
%
% calculate the kaiser-bessel kernel without normalization
%
% (c) Xue Feng 2012
% University of Virginia

x = beta*sqrt(1-(2*u/w).^2);	% Argument - see Jackson '91.

nx = length(x);
x = abs(x);
for i = 1:nx,
    if x(i) < 3.75,
        temp = (x(i)/3.75)^2;
        y(i) = 1.0+temp*(3.5156229+temp*(3.0899424+temp*(1.2067492+temp*(0.2659732+temp*(0.360768e-1+temp*0.45813e-2)))));
    else
        temp = 3.75/x(i);
        y(i) = (exp(x(i))/sqrt(x(i)))*(0.39894228+temp*(0.1328592e-1+temp*(0.225319e-2+temp*(-0.157565e-2+temp*(0.916281e-2+temp*(-0.2057706e-1+temp*(0.2635537e-1+temp*(-0.1647633e-1+temp*0.392377e-2))))))));
    end
end
y = y/w;

end


function m = kb_ft(gridsize,imagesize,kwidth,beta)
% function m = kb_ft(gridsize,imagesize,kwidth,beta)
% 
% Fourier Transform of Kaiser window function
% returns the deapodization value for imagesize*imagesize
% imagesize = gridsize/overgridfactor
%     gridsize -- image size*overgridfactor
%     kwidth -- length of the window of Kaiser window function.
% parameters for the Kaiser window function, i.e, the kernel function
%
% (c) Xue Feng 2012
% University of Virginia

L = kwidth;
n = gridsize;
kbmax = sinh(sqrt(beta^2))/sqrt(beta^2);
for i = 1:imagesize,
    mi = sqrt(beta^2-(pi*L*(i-imagesize/2-1)/n)^2);
    if mi == 0,
        deconx(i) = 1/kbmax;
    else
        deconx(i) = sinh(mi)/mi/kbmax;
    end
end
for j = 1:imagesize,
    mj = sqrt(beta^2-(pi*L*(j-imagesize/2-1)/n)^2);
    if mj == 0,
        decony(j) = 1/kbmax;
    else
        decony(j) = sinh(mj)/mj/kbmax;
    end
end

m = decony'*deconx;

end