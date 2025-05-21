%BPSK Modulation and Demodulation Simulation

clear all 
clc
%% parameters 
T=1; % assume one bit /sec 
N = 10000; % Number of bits
data = randi([0,1],1,N); % Generate random binary stream [ generata integer values from [0,1]  in matrex 1*N ]
vp=1;% vp is peak voltage of NRZ 
fs=1000 ;% number of samples per bit
u= fs*N; % total samples 
t = 0:1/(fs):(N-(1/(fs)));   % Time vector t = 0 : 0.001 : N-0.001
fc = 5 ; % Carrier frequency
f = -fs/2:1/N:fs/2-1/N; % freq vector samples 
k=1;  % bit per symbol 
%% Transmitter

% baseband [ output NRZ 0 --> -1 and 1 --> 1 ]
baseband_modulated=[];
for i =1:N
    if data(i)==1
        baseband_modulated=[baseband_modulated ones(1,fs)*(vp)];
    elseif data(i)==0
        baseband_modulated=[baseband_modulated ones(1,fs)*(-vp)];
    end
end

% Passband modulation
carrier = (sqrt(2/T))*cos(2*pi*fc*t); % Carrier signal
passband_modulated = baseband_modulated .* carrier; % Passband modulation

% bit stream impulses plot 
subplot(3,1,1);
stem(data);
xlabel('Time (sec)');
ylabel('Amplitude (volts)');
title('impulses of bits to be transmited ');
ylim([-2 2]);
grid on;
% output of NRZ ploting 

subplot(3,1,2);
plot(baseband_modulated);
xlabel('Time (mel sec)');
ylabel('Amplitude (volts)');
title('generated NRZ signal');
ylim([-2 2]);
grid on;

% Plot passband modulated signal and baseband modulated data

subplot(3,1,3);
plot(t,passband_modulated);
title('Passband Modulated Signal');
xlabel('Time(sec)');
ylabel('Amplitude');
ylim([-2 2]);
grid on;


% Constellation diagram of baseband modulated data
scatterplot(baseband_modulated);
title('Constellation Diagram of Baseband Modulated Data for BPSK');
xlabel('In-Phase Amplitude');
ylabel('Quadrature Amplitude');
grid on;

%{
subplot(2,1,2);
plot(t,baseband_modulated);
title('Baseband Modulated Data');
xlabel('Time');
ylabel('Amplitude');
ylim([-2 2]);
grid on;
%}
% Spectrum of passband modulated signal
figure;
spectrum = fftshift(fft(passband_modulated));
plot(f,abs(spectrum)/N);
title('Spectrum of Passband Modulated Signal');
xlabel('Frequency (Hz)');
ylabel('Magnitude');
xlim([-20 20]);
grid on;

% Bandwidth estimation
[bw,flo,fhi,power] = obw(passband_modulated,fs);
%bw = obw(passband_modulated,fs);
bw_estimate = bw;
disp(['Estimated Bandwidth of the Passband Modulated Signal: ', num2str(bw_estimate), ' Hz']);
%% channel and receiver 
% channel
EbN0_range = 0:5:20; % Eb/N0 range
%BER = zeros(1,length(EbN0_range)); % Placeholder for Bit Error Rate
BER =[];
SNR=[];
for i = 1:length(EbN0_range)
    % Add AWGN
    SNR(i) = EbN0_range(i) + 10*log10(k); % Convert Eb/N0 to SNR    
    received_signal = awgn(baseband_modulated,SNR(i));     
    figure;
    subplot(2,1,1);
    plot(t,received_signal);
    title(['Received Signal at SNR = ' num2str(SNR(i)) 'db']);
    xlabel('Time (sec)');
    ylabel('Amplitude');
    ylim([-2 2]);
 
    % detector    
    detected=sign(received_signal);
    sample_index=1:fs:length(received_signal);
    detected_bit=detected(sample_index) > 0 ;   
    %figure;
    subplot(2,1,2);
    stem(detected_bit);
    title('Impulses of Received bits');
    xlabel('Time (seconds)-->');
    ylabel('Amplitude (volts)');
    ylim([-2 2]);
    % Calculate Bit Error Rate
    errors = sum( data ~= detected_bit);
    BER = [BER (errors/N)+10e-10] ;
end
%%  plot Theoretical and BER with EbN0
% Theoretical BER calculation
theory_BER = 0.5*erfc(sqrt(10.^(EbN0_range/10)));
% Plot BER vs Eb/N0
figure;
semilogy(EbN0_range, BER, 'o-', 'DisplayName', 'Simulated BER');
hold on;
semilogy(EbN0_range, theory_BER, 'r--', 'DisplayName', 'Theoretical BER');
title('Bit Error Rate vs Eb/N0');
xlabel('Eb/N0 (dB)');
ylabel('Bit Error Rate');
ylim([10e-10  10e1])
grid on;
legend;
