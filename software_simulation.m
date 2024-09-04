%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Multiple Access Techniques %%%%
%%%%%%%%%%% OMA and NOMA %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear all; close all

% Some global constants %
N = 1e5;
SNR = -20:50; % SNR range in dB
snr = db2pow(SNR); % SNR range in linear scale

%%%%% OMA Transmitter %%%%%

% BPSK Signalling % 
x1 = randi([0 1],1,N);
x2 = randi([0 1],1,N);
fprintf('Data of user 1 is %s\n',sprintf('%d', x1));
fprintf('Data of user 2 is %s\n',sprintf('%d', x2));

x1bpsk = 2 * x1 - 1;
x2bpsk = 2 * x2 - 1;

% Assigning Power Coefficients to User 1 and 2 %
a1 = 0.75;
a2 = 0.25;

% Superposition coding %
x = sqrt(a1) * x1bpsk + sqrt(a2)*x2bpsk;

% Plotting data signals %
figure(1);
subplot(211);
stem(x1bpsk(1:100),'LineWidth',3);
title('Data of user 1 (x_1)');
grid on;
subplot(212);
stem(x2bpsk(1:100),'LineWidth',3);
title('Data of user 2 (x_2)');
grid on;

% Plotting superposition signal %
figure(2);
stem(x(1:100),'LineWidth',3);
title('Superposition coded signal at transmitter');
grid on;

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
end

% Plotting %
semilogy(SNR, ber1, 'linewidth', 3);
hold on;
semilogy(SNR, ber2, 'linewidth', 3);
grid on;
legend('\alpha_1 = 0.75', '\alpha_2 = 0.25');
xlabel('SNR (dB)');
ylabel('BER');
title('BER for AWGN');