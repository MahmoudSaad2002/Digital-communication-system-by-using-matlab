% BFSK Modulation and Demodulation Simulation
clear all;
clc;

%% Parameters
T = 1; % Bit duration (seconds)
N = 5000 ; % Number of bits
data = randi([0, 1], 1, N); % Generate random binary stream
fs = 1000; % Number of samples per bit
u = fs * N; % Total samples
t = 0:1/fs:(N*T - 1/fs); % Time vector

A=1; % amp of sin
Eb=(A*A)/2; % energy per bit
sqrt_Eb=sqrt(Eb);

% Frequencies for BFSK modulation
fc0 = 5; % Frequency for bit 0
fc1 = 10; % Frequency for bit 1

%% Transmitter
% Baseband modulated data (BFSK)
 
baseband_modulated = [];
for i =1:N
    if data(i)==1
        baseband_modulated=[baseband_modulated ones(1,fs)];
    elseif data(i)==0
        baseband_modulated=[baseband_modulated zeros(1,fs)];
    end
end

% passband
passband_modulated = [];
for i = 1:N
    if data(i) == 0
        passband_modulated = [passband_modulated  A*sqrt(2/T)*sin(2 * pi * fc0 * t((i-1)*fs + 1:i*fs))]; % Frequency for 0
    else
        passband_modulated = [passband_modulated  A*sqrt(2/T)*sin(2 * pi * fc1 * t((i-1)*fs + 1:i*fs))]; % Frequency for 1
    end
end
%}

% Bit stream impulses plot 
subplot(3,1,1);
stem(data);
xlabel('Time (sec)');
ylabel('Amplitude (volts)');
title('Impulses of Bits to be Transmitted');
ylim([-0.5 1.5]);
% output of unipolar NRZ 
grid on;
% baseband ploting 
subplot(3,1,2);
plot(t,baseband_modulated);
xlabel('Time (mel sec)');
ylabel('Amplitude (volts)');
title('generated unipolar NRZ signal');
ylim([-2 2]);
grid on;

% Output of BFSK plotting 
subplot(3,1,3);
plot(t, passband_modulated);
xlabel('Time (s)');
ylabel('Amplitude (volts)');
title('Generated BFSK Signal');
ylim([-1.5 1.5]);
grid on;

% Constellation diagram (for BFSK, 0 and 1 represented as two points)
% Constellation points
constellation = zeros(2, 2); % Preallocate for two points
constellation(1, :) = [sqrt_Eb, 0]; % Point for bit 0
constellation(2, :) = [0, sqrt_Eb]; % Point for bit 1
% Plot constellation diagram
figure;
scatter(constellation(1, :), constellation(2, :), 'filled');
title('Constellation Diagram for BFSK');
xlabel('In-Phase Amplitude');
ylabel('Quadrature Amplitude');
%xlabel('Frequency (Hz)');
%ylabel('Amplitude');
xlim([0, 2]); % Adjust limits as needed
ylim([0, 2]); % Adjust limits as needed
grid on;
text(sqrt_Eb, 0.2, '0', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(0.2, sqrt_Eb , '1', 'HorizontalAlignment', 'center', 'FontSize', 12);

%% Spectrum Analysis
spectrum = fftshift(fft(passband_modulated)); % Spectrum of the BFSK signal
f = -fs/2:fs/length(spectrum):fs/2-fs/length(spectrum); % Frequency vector

figure;
plot(f, abs(spectrum)/N);
title('Spectrum of Passband Modulated Signal (BFSK)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-20 20]);
grid on;

% Bandwidth estimation
[bw, ~, ~, ~] = obw(baseband_modulated, fs);
disp(['Estimated Bandwidth of the Passband Modulated Signal: ', num2str(bw), ' Hz']);

%% Channel: Add AWGN  and receiver 
EbN0_range = 0:5:20; % Eb/N0 range
BER = []; % Placeholder for Bit Error Rate
k = 1;  % bits per symbol 
SNR=[];
for i = 1:length(EbN0_range)
    % Convert Eb/N0 to SNR
    SNR(i) = EbN0_range(i) + 10*log10(k); 
    % Add noise
    received_signal = awgn(baseband_modulated, SNR(i), 'measured'); 
    figure;
    subplot(2,1,1);
    plot(t,received_signal);
    title(['Received Signal at SNR = ' num2str(SNR(i)) 'db']);
    xlabel('Time (sec)');
    ylabel('Amplitude');
    ylim([-2 2]);
    grid on ;
    % detector
    sample_index=1:fs:length(received_signal);
    detected_bit=received_signal(sample_index) > 0.5 ;
          
    subplot(2,1,2);
    stem(detected_bit);
    title('Impulses of Received bits');
    xlabel('Time (seconds)-->');
    ylabel('Amplitude (volts)');
    ylim([-2 2]);
    grid on ;
    % Calculate Bit Error Rate
    errors = sum( data ~= detected_bit);
    BER = [BER (errors/N)+1e-10] ;
end

% Theoretical BER calculation for BFSK
theory_BER = 0.5 * erfc(sqrt(10.^(EbN0_range/10)));

% Plot BER vs Eb/N0
figure;
semilogy(EbN0_range, BER, 'o-', 'DisplayName', 'Simulated BER');
hold on;
semilogy(EbN0_range, theory_BER, 'r--', 'DisplayName', 'Theoretical BER');
title('Bit Error Rate vs Eb/N0 (BFSK)');
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
ylim([1e-10 1]);
grid on;
legend;