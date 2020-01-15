
%%
% Initializing parameters
 function [papr,snr,MM,r] = weight_ofdm() 
OFDM.N=input('Size of OFDM Symbol N = ');
OFDM.m=input('Number of OFDM symbols to be simulated m = ');
OFDM.M=input('Size of Alphabet M = ');
OFDM.L=input('Up-sampling factor L = ');
OFDM.PoQ=input('Type of Mapping (1 for PSK) and (2 for QAM) = ');
OFDM.Phase_Offset=input('constellation phase offset = ');
OFDM.Symbol_Order=input('constellation Symbol Order (1 for Binary) and (2 for Gray) = ');
OFDM.Ncp=input('size of cyclic prefix samples Ncp = ');
MM=1:OFDM.m;
%% Transmitter
% Creating Baseband modems Tx/Rx
if OFDM.Symbol_Order == 1
    OFDM.Symbol_Order = 'binary';
else
    OFDM.Symbol_Order = 'gray';
end
if OFDM.PoQ == 1
    hTx = modem.pskmod('M',OFDM.M,'PhaseOffset',OFDM.Phase_Offset,'SymbolOrder',OFDM.Symbol_Order);
    hRx = modem.pskdemod('M',OFDM.M,'PhaseOffset',OFDM.Phase_Offset,'SymbolOrder',OFDM.Symbol_Order);
else
    hTx = modem.qammod('M',OFDM.M,'PhaseOffset',OFDM.Phase_Offset,'SymbolOrder',OFDM.Symbol_Order);
    hRx = modem.qamdemod('M',OFDM.M,'PhaseOffset',OFDM.Phase_Offset,'SymbolOrder',OFDM.Symbol_Order);
end
% data generation
OFDM.DATA=randi([0 OFDM.M-1],OFDM.m,OFDM.N/OFDM.L);
% figure;
% plot(OFDM.DATA);
% Mapping
OFDM.Dmap=modulate(hTx,OFDM.DATA);
% Serial to Parallel
OFDM.parallel=OFDM.Dmap.';
% Oversampling 
OFDM.upsampled=upsample(OFDM.parallel,OFDM.L);
% Amplitude modulation (IDFT using fast version IFFT)
ofdm.am=ifft(OFDM.upsampled,OFDM.N);
% Parallel to serial
ofdm.serial=ofdm.am.';
% Cyclic Prefixing
ofdm.CP_part=ofdm.serial(:,end-OFDM.Ncp+1:end); % this is the Cyclic Prefix part to be appended.
ofdm.cp=[ofdm.CP_part ofdm.serial];
%%  Reciever
% % Adding Noise using AWGN
% SNRstart=0;
% SNRincrement=2;
% SNRend=30;
% c=0;
% r=zeros(size(SNRstart:SNRincrement:SNRend));
% for snr=SNRstart:SNRincrement:SNRend
%     c=c+1;
%     ofdm.noisy=awgn(ofdm.cp,snr,'measured');
% % Remove cyclic prefix part
%     ofdm.cpr=ofdm.noisy(:,OFDM.Ncp+1:OFDM.N+OFDM.Ncp); %remove the Cyclic prefix
% % serial to parallel
%     ofdm.parallel=ofdm.cpr.';
% % Amplitude demodulation (DFT using fast version FFT)
%     OFDM.amdemod=fft(ofdm.parallel,OFDM.N);
% % Down-Sampling
% OFDM.downsampled=downsample(OFDM.amdemod,OFDM.L);
% % Parallel to serial
%     OFDM.rserial=OFDM.downsampled.';
% % Baseband demodulation (Un-mapping)
%     OFDM.Umap=demodulate(hRx,OFDM.rserial);
% % Calculating the Symbol Error Rate
%     [n, r(c)]=symerr(OFDM.DATA,OFDM.Umap);
%     %disp(['SNR = ',num2str(snr),' step: ',num2str(c),' of ',num2str(length(r))]);
% end
% snr=SNRstart:SNRincrement:SNRend;
% % Plotting SER vs SNR
% semilogy(snr,r,'-ok','linewidth',2,'markerfacecolor','r','markersize',8,'markeredgecolor','b');grid;
% title('OFDM BIT Error Rate vs SNR');
% ylabel('Bit Error Rate');
% xlabel('SNR [dB]');
% legend(['SER N = ', num2str(OFDM.N),' ',num2str(hTx.M),'-',hTx.type]);
%  ccdf1=zeros(1,32);
 %% WEIGHTED OFDM PAPR CALCULATIONS


for i=1:OFDM.m
  ofdmpapr=OFDM.DATA;
    j=ifft(ofdmpapr);
    k=j(:,1:32);
     x2=(abs(k));
     m1=mean(x2(i));
end 

   
 
 %L=OFDM.L;
 %N=OFDM.N;
%Clipping the signals above threshold(here 70% of original value)
for i=1:OFDM.m
 x_mag(i)=x2(i);
 x_max=0.9*max(x_mag);
    
    if(x_mag(i)>x_max)
    x_mag(i)=x_max;
    else
       x_mag(i)=x_mag(i);

    end
    clip=fft(x_mag,OFDM.N);
    z=real(clip);
 p(i)=10*log10(double(x_mag(i)/m1));
 
end; 
figure;
 semilogy(MM,sort(abs(p),'descend'),'-');
 xlabel('time variations in sec'),ylabel('PAPR')
 grid on;
 hold on;
 
 for i=1:OFDM.m
    %z=pap(i);

    % 
%    for i=1:OFDM.m;
% %z1=[ifftdata12(i)(1:ifftn/2),zeros(1,3*ifftn),ifftdata(i)(ifftn/2+1:ifftn)];
% %4 oversample
     rect12 = rectwin(OFDM.m);
    rect21 = rectwin(OFDM.N);
    rect21x = rect21(1:OFDM.m)';
    rectx = rect12*rect21x;
    rectd = padarray(rectx,[0,OFDM.m],'replicate','post');
    x1 = randi(100,OFDM.m,1);
    xq = padarray(x1,[0,OFDM.N-1],'replicate','post');
    ofdmpapr = OFDM.DATA.*rectd.*xq;
     w1(:,i)=ifft(ofdmpapr(:,i));            
     w1(:,i)=w1(:,i)*OFDM.N;
     x2=(abs(w1(:,i))).^2;
     m1=mean(x2);
% %     v1=max(x2);
    papr(i)=10*log10(double(x2(i)/m1));
    for s1 = 1:length(papr)
    ccdf(s1) = double(length(find(papr >= papr(s1)))/OFDM.m); % # of samples exceeding papr(i)
    end
 end
% % [c,d]=butter(6,0.5);
% % w3=filter(c,d,w2);
% % MM=1:.1:10;
% %  ccdf0=ccdf1./n1;
% %  ccdf3=ccdf2./n1;
% % ccdf1 = ccdf1./OFDM.m;
%figure; 
h=papr;
semilogy(MM,sort(abs(h),'descend'),'-');
 xlabel('time variations in sec'),ylabel('PAPR')
 grid on;
 hold off;  
 legend('clipper=p','weighted=h');
 figure;
 semilogy(MM,sort(ccdf,'descend'),'*');
 xlabel('Time Variations in (sec)'),ylabel('ccdf')
% grid on;
%  hold off; 
%end
new=OFDM.DATA;
r=real(new);
f=r(:,1:15);
l=(abs(f)).^2;
g=z(:,1:15);
u=(abs(g)).^2;

SNRstart=0;
 SNRincrement=2;
 SNRend=30;
 c=0;
 r=zeros(size(SNRstart:SNRincrement:SNRend));
 for snr=SNRstart:SNRincrement:SNRend
  c=c+1;  
   ofdm.noisy=awgn(ofdm.cp,snr,'measured');
% Remove cyclic prefix part
    ofdm.cpr=ofdm.noisy(:,OFDM.Ncp+1:OFDM.N+OFDM.Ncp); %remove the Cyclic prefix
% serial to parallel
    ofdm.parallel=ofdm.cpr.';
% Amplitude demodulation (DFT using fast version FFT)
    OFDM.amdemod=fft(ofdm.parallel,OFDM.N);
% Down-Sampling
OFDM.downsampled=downsample(OFDM.amdemod,OFDM.L);
% Parallel to serial
    OFDM.rserial=OFDM.downsampled.';
% Baseband demodulation (Un-mapping)
    OFDM.Umap=demodulate(hRx,OFDM.rserial);
   
     
% Calculating the Symbol Error Rate
     [n, r(c)]=symerr(OFDM.DATA,OFDM.Umap);
     disp(['SNR = ',num2str(snr),' step: ',num2str(c),' of ',num2str(length(r))]);
     
 
 end

snr=SNRstart:SNRincrement:SNRend;
 % Plotting SER vs SNR
 semilogy(snr,r,'-ok','linewidth',2,'markerfacecolor','r','markersize',8,'markeredgecolor','b');grid;
 title('Clipped OFDM BIT Error Rate vs SNR');
 ylabel('Bit Error Rate');
 xlabel('SNR [dB]');
 %legend(['SER N = ', num2str(OFDM.N),' ',num2str(hTx.M),'-',hTx.type]);
 hold on;

 %%  Reciever
% % Adding Noise using AWGN
SNRstart=0;
SNRincrement=2;
SNRend=30;
c=0;
k=zeros(size(SNRstart:SNRincrement:SNRend));
for snr=SNRstart:SNRincrement:SNRend
    c=c+1;
    ofdm.noisy=awgn(ofdm.cp,snr,'measured');
% Remove cyclic prefix part
    ofdm.cpr=ofdm.noisy(:,OFDM.Ncp+1:OFDM.N+OFDM.Ncp); %remove the Cyclic prefix
% serial to parallel
    ofdm.parallel=ofdm.cpr.';
% Amplitude demodulation (DFT using fast version FFT)
    OFDM.amdemod=fft(ofdm.parallel,OFDM.N);
% Down-Sampling
OFDM.downsampled=downsample(OFDM.amdemod,OFDM.L);
% Parallel to serial
    OFDM.rserial=OFDM.downsampled.';
% Baseband demodulation (Un-mapping)
    OFDM.Umap=demodulate(hRx,OFDM.rserial);
% Calculating the Symbol Error Rate
    [n, k(c)]=biterr(OFDM.DATA,OFDM.Umap);
    %disp(['SNR = ',num2str(snr),' step: ',num2str(c),' of ',num2str(length(r))]);
end
snr=SNRstart:SNRincrement:SNRend;
% Plotting SER vs SNR
%6figure;
semilogy(snr,k,'-ok','linewidth',2,'markerfacecolor','r','markersize',8,'markeredgecolor','b');grid;
title('Weighted OFDM BIT Error Rate vs SNR');
ylabel('Bit Error Rate');
xlabel('SNR [dB]');

legend('CLIPPING OFDM BER VS SNR=r','weighted OFDM BER VS SNR=k');
hold off;