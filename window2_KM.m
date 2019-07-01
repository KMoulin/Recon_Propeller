function w=window2_KM(N,M,w_func,width)

w=zeros(N,M);
wc=zeros(N,1);
wr=zeros(M,1);
wc(N/2-N*width/2+1:1:N/2+N*width/2)=window(w_func,N*width);
wr(M/2-M*width/2+1:1:M/2+M*width/2)=window(w_func,M*width);

[maskr,maskc]=meshgrid(wr,wc);

%maskc=repmat(wc,1,M); Old version
%maskr=repmat(wr',N,1);

w=maskr.*maskc;

end