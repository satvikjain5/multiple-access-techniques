clc; clear all;

% Some global constants %
N = 5.2e4;
SNR = 0:35; % SNR range in dB
snr = db2pow(SNR); % SNR range in linear scale

%%%%% OMA Transmitter %%%%%

% BPSK Signalling % 
x1 = randi([0 1],1,N);
x2 = randi([0 1],1,N);
%fprintf('Data of user 1 is %s\n',sprintf('%d', x1));
%fprintf('Data of user 2 is %s\n',sprintf('%d', x2));

x1bpsk = 2 * x1 - 1;
x2bpsk = 2 * x2 - 1;

% Assigning Power Coefficients to User 1 and 2 %
a1 = 0.75;
a2 = 0.25;

% Superposition coding %
x = sqrt(a1) * x1bpsk + sqrt(a2)*x2bpsk;

%%%%% OMA Receiver %%%%%

for i = 1:length(snr)    
    x_awgn = awgn(x,SNR(i),'measured');
    x1_received = ones(1,N);
    x1_received(x_awgn < 0) = 0;
    x1_received_ASK = ones(1,N);
    x1_received_ASK(x_awgn < 0) = 0;
    x1_received_ASK(x1_received == 0) = -1; %Remodulating BPSK signal
    remaining_signal = x_awgn - sqrt(a1) * x1_received_ASK;
    x2_received = zeros(1,N);
    x2_received(remaining_signal > 0) = 1;
    ber1(i) = biterr(x1,x1_received) / N;
    ber2(i) = biterr(x2,x2_received) / N;
    if (SNR(i) == 10)
        disp(ber1(i))
    end
end


% Plotting %
figure;
semilogy(SNR, ber1,'rx-','linewidth', 1);
hold on;
semilogy(SNR, ber2,'co-', 'linewidth', 1);
grid on;
legend('UE 1 \alpha_1 = 0.75', 'UE 2 \alpha_2 = 0.25');
xlabel('SNR (dB)');
ylabel('BER');
title('BER for BPSK using NOMA in a 10-tap Rayleigh channel');

