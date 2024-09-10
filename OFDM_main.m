clc;
clear all;
close all;


N = 1e5;
SNR = 0:30; % SNR range in dB
snr = db2pow(SNR); % SNR range in linear scale


nfft = 128; % Number of data carriers
M = 16; % Modulation order for 16QAM
cplen = 16; % Cyclic prefix length
nSym = 5;
nt = 1;
dataIn = randi([0 M-1],nfft,nSym,nt);
qamSig = qammod(dataIn,M,'UnitAveragePower',true);
ofdmSig = ofdmmod(qamSig,nfft,cplen);

for i = 1:length(snr)
    awgnSig = awgn(ofdmSig,SNR(i),'measured');
    ofdmdemodSig = ofdmdemod(awgnSig,nfft,cplen);
    qamdemodSig = qamdemod(ofdmdemodSig,M,'UnitAveragePower',true);
    ber1(i) = biterr(dataIn,qamdemodSig)/nfft;
end

semilogy(SNR,ber1,'linewidth',3);
grid on;
xlabel('SNR (dB)');
ylabel('BER');
title('OFDM BER for AWGN');