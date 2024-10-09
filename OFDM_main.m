
clc; clear all;
nFFT        = 64; % fft size
nDSC        = 52; % number of data subcarriers
nBitPerSym  = 52; % number of bits per OFDM symbol (same as the number of subcarriers for BPSK)
nSym        = 10^3; % number of symbols

EbN0dB      = [0:35]; % bit to noise ratio
EsN0dB      = EbN0dB + 10*log10(nDSC/nFFT) + 10*log10(64/80); % converting to symbol to noise ratio

for ii = 1:length(EbN0dB)

   % Transmitter
   ipBit1 = rand(1,nBitPerSym*nSym) > 0.5; % random 1's and 0's
   ipMod1 = 2*ipBit1-1; % BPSK modulation 0 --> -1, 1 --> +1
   ipMod1 = reshape(ipMod1,nBitPerSym,nSym).'; % grouping into multiple symbolsa
   ipBit2 = rand(1,nBitPerSym*nSym) > 0.5; % random 1's and 0's
   ipMod2 = 2*ipBit2-1; % BPSK modulation 0 --> -1, 1 --> +1
   ipMod2 = reshape(ipMod2,nBitPerSym,nSym).'; % grouping into multiple symbolsa

   % Assigning modulated symbols to subcarriers from [-26 to -1, +1 to +26]
   xF1 = [zeros(nSym,6) ipMod1(:,[1:nBitPerSym/2]) zeros(nSym,1) ipMod1(:,[nBitPerSym/2+1:nBitPerSym]) zeros(nSym,5)] ;
   xF2 = [zeros(nSym,6) ipMod2(:,[1:nBitPerSym/2]) zeros(nSym,1) ipMod2(:,[nBitPerSym/2+1:nBitPerSym]) zeros(nSym,5)] ;
   
   % Taking FFT, the term (nFFT/sqrt(nDSC)) is for normalizing the power of transmit symbol to 1 
   xt1 = (nFFT/sqrt(nDSC))*ifft(fftshift(xF1.')).';
   xt2 = (nFFT/sqrt(nDSC))*ifft(fftshift(xF2.')).';


   % Appending cylic prefix
   xt1 = [xt1(:,[49:64]) xt1];
   xt2 = [xt2(:,[49:64]) xt2];

   % multipath channel
   nTap = 10;
   ht = 1/sqrt(2)*1/sqrt(nTap)*(randn(nSym,nTap) + j*randn(nSym,nTap));
   
   % computing and storing the frequency response of the channel, for use at recevier
   hF = fftshift(fft(ht,64,2));

   % convolution of each symbol with the random channel
   for jj = 1:nSym
      xht1(jj,:) = conv(ht(jj,:),xt1(jj,:));
      xht2(jj,:) = conv(ht(jj,:),xt2(jj,:));
   end
   xt1 = xht1;
   xt2 = xht2;

   % Concatenating multiple symbols to form a long vector
   xt1 = reshape(xt1.',1,nSym*(80+nTap-1));
   xt2 = reshape(xt2.',1,nSym*(80+nTap-1));
   % Gaussian noise of unit variance, 0 mean
   nt = 1/sqrt(2)*[randn(1,nSym*(80+nTap-1)) + j*randn(1,nSym*(80+nTap-1))];

   % Adding noise, the term sqrt(80/64) is to account for the wasted energy due to cyclic prefix
   yt1 = sqrt(80/64)*xt1 + 10^(-EsN0dB(ii)/20)*nt;
   yt2 = sqrt(80/64)*xt2 + 10^(-EsN0dB(ii)/20)*nt;

   % Receiver
   yt1 = reshape(yt1.',80+nTap-1,nSym).'; % formatting the received vector into symbols
   yt2 = reshape(yt2.',80+nTap-1,nSym).';
   yt1 = yt1(:,[17:80]); % removing cyclic prefix
   yt2 = yt2(:,[17:80]);

   % converting to frequency domain
   yF1 = (sqrt(nDSC)/nFFT)*fftshift(fft(yt1.')).'; 
   yF2 = (sqrt(nDSC)/nFFT)*fftshift(fft(yt2.')).'; 

   % equalization by the known channel frequency response
   yF1 = yF1./hF;
   yF2 = yF2./hF;

   % extracting the required data subcarriers
   yMod1 = yF1(:,[6+[1:nBitPerSym/2] 7+[nBitPerSym/2+1:nBitPerSym] ]);
   yMod2 = yF2(:,[6+[1:nBitPerSym/2] 7+[nBitPerSym/2+1:nBitPerSym] ]);

   % BPSK demodulation
   % +ve value --> 1, -ve value --> -1
   ipModHat1 = 2*floor(real(yMod1/2)) + 1;
   ipModHat1(find(ipModHat1>1)) = +1;
   ipModHat1(find(ipModHat1<-1)) = -1;

   ipModHat2 = 2*floor(real(yMod2/2)) + 1;
   ipModHat2(find(ipModHat2>1)) = +1;
   ipModHat2(find(ipModHat2<-1)) = -1;
   % converting modulated values into bits
   ipBitHat1 = (ipModHat1+1)/2;
   ipBitHat1 = reshape(ipBitHat1.',nBitPerSym*nSym,1).';
   ipBitHat2 = (ipModHat2+1)/2;
   ipBitHat2 = reshape(ipBitHat2.',nBitPerSym*nSym,1).';

   % counting the errors
   nErr1(ii) = size(find(ipBitHat1 - ipBit1),2);
   nErr2(ii) = size(find(ipBitHat2 - ipBit2),2);

end

simBer1 = nErr1/(nSym*nBitPerSym);
simBer2 = nErr2/(nSym*nBitPerSym);
EbN0Lin = 10.^(EbN0dB/10);
theoryBer = 0.5.*(1-sqrt(EbN0Lin./(EbN0Lin+1)));

close all; figure

semilogy(EbN0dB,simBer1,'rx-','LineWidth',1);
hold on;
semilogy(EbN0dB,simBer2,'co-','LineWidth',1);
axis([0 35 10^-5 1])
grid on
legend('User 1','User 2');
xlabel('SNR, dB')
ylabel('BER')
title('BER for BPSK using OFDM in a 10-tap Rayleigh channel')

