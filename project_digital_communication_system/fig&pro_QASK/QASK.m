% 4-ARY ASK Modulation and Demodulation Simulation
clear all;
clc;

%% Parameters
T = 1; % symbol duration (seconds)
N = 10000 ; % Number of symbols
data = randi([0, 3], 1, N); % Generate random 2-bit symbols (0 to 3)
fs = 1000; % Number of samples per bit
u = fs * N; % Total samples
t = 0:1/fs:(N*T - 1/fs); % Time vector
fc=5;

A = 1; % Maximum amplitude
Eb = (A^2)/2; % Energy per bit
sqrt_Eb = sqrt(Eb);

% Amplitude levels for 4-ASK
amplitude_levels = [A/4, A/2, 3*A/4, A]; % Define the four amplitude levels for symbols 0, 1, 2, 3

Amp=amplitude_levels; % amp of sin
Eb=(Amp.*Amp)/2; % energy per bit
sqrt_Eb=sqrt(Eb);

%% Transmitter
% Baseband modulated data (4-ASK)
baseband_modulated = [];
for i = 1:N
    baseband_modulated = [baseband_modulated amplitude_levels(data(i) + 1) * ones(1, fs)];
end


% Passband modulation
carrier = (sqrt(2/T))*cos(2*pi*fc*t); % Carrier signal
passband_modulated = baseband_modulated .* carrier; % Passband modulation


% symbols plotting 

subplot(3,1,1)
stem(data);
xlabel('Symbol Index');
ylabel('Amplitude (volts)');
title('Symbols to be Transmitted');
ylim([-0.5 3.5]);
grid on;

% Output of base_band 4-ASK plotting 
subplot(3,1,2)
plot(t, baseband_modulated);
xlabel('Time (s)');
ylabel('Amplitude (volts)');
title('Generated baseband 4-ASK Signal');
ylim([-0.5 1.5]);
grid on;

%output modulation 
subplot(3,1,3)
plot(t,passband_modulated);
title('Passband Modulated Signal');
xlabel('Time');
ylabel('Amplitude');
ylim([-2 2]);
grid on;

% Create the constellation points
constellation = zeros(2, 4); % Preallocate for four points
constellation(1, :) = sqrt_Eb; % Amplitude levels on the x-axis
constellation(2, :) = zeros(1, 4); % All points on the y-axis (0 for ASK)

% Plot the constellation diagram
figure;
scatter(constellation(1, :), constellation(2, :), 'filled');
title('Constellation Diagram for 4-ARY ASK');
xlabel('Energy');
% ylabel('Quadrature Amplitude');
xlim([-0.5, A + 0.5]); % Adjust limits as needed
ylim([-0.5, 1]); % Adjust limits as needed
grid on;

% Label the constellation points
text(constellation(1, 1), 0.1, '0', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(constellation(1, 2), 0.1, '1', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(constellation(1, 3), 0.1, '2', 'HorizontalAlignment', 'center', 'FontSize', 12);
text(constellation(1, 4), 0.1, '3', 'HorizontalAlignment', 'center', 'FontSize', 12);



%% Spectrum Analysis
spectrum = fftshift(fft(passband_modulated)); % Spectrum of the 4-ASK signal
f = -fs/2:fs/length(spectrum):fs/2-fs/length(spectrum); % Frequency vector

figure(2);
plot(f, abs(spectrum)/N);
title('Spectrum of Passband Modulated Signal (4-ASK)');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-10 10]);
grid on;

% Bandwidth estimation
[bw, ~, ~, ~] = obw(baseband_modulated, fs);
disp(['Estimated Bandwidth of the 4-ASK Signal: ', num2str(bw), ' Hz']);

%% Channel: Add AWGN and receiver 
EbN0_range = 0:5:20; % Eb/N0 range
BER = []; % Placeholder for Bit Error Rate
k = 2; % bits per symbol (4-ASK)
SNR=[];
for i = 1:length(EbN0_range)
    % Convert Eb/N0 to SNR
    SNR(i) = EbN0_range(i) + 10*log10(k); 
    % Add noise
    received_signal = awgn(baseband_modulated, SNR(i), 'measured');
    figure;
    subplot(2,1,1);
    plot(t, received_signal);
    title(['Received Signal at SNR = ' num2str(SNR(i)) 'db']);
    xlabel('Time (s)');
    ylabel('Amplitude');
    ylim([-0.5 1.5]);
    % grid on ;
    
    
    
    % Detector
    detected_symbols = zeros(1, N);
    for j = 1:N
        % Sample the received signal at the symbol intervals
        sample = received_signal((j-1)*fs + 1);
        % Determine the detected symbol based on the amplitude
        [~, detected_symbols(j)] = min(abs(sample - amplitude_levels));
        detected_symbols(j) = detected_symbols(j) - 1; % Adjust index (0-3)
    end
    subplot(2,1,2);
    stem(detected_symbols);
    title('Impulses of Received bits');
    xlabel('Time (seconds)-->');
    ylabel('Amplitude (volts)');
    ylim([-1 4]);
   % grid on ;
    
    % Calculate Bit Error Rate
    errors = sum(data ~= detected_symbols);
    BER = [BER (errors/N)+1e-30];
end

% Theoretical BER calculation for 4-ASK
%theory_BER = 0.5 * erfc(sqrt(10.^(EbN0_range/10)));
theory_BER = 3 * qfunc(sqrt(3 * 10.^(EbN0_range / 20)));

% Plot BER vs Eb/N0
figure;
semilogy(EbN0_range, BER, 'o-', 'DisplayName', 'Simulated BER');
hold on;
semilogy(EbN0_range, theory_BER, 'r--', 'DisplayName', 'Theoretical BER');
title('Bit Error Rate vs Eb/N0 (4-ASK)');
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
ylim([1e-10 1]);
grid on;
legend;