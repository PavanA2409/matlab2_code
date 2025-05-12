% Digital Modulation Techniques Visualization
clc;
clear all; 
close all; 

%% Simulation Parameters
fs = 10000;         % Sampling frequency (Hz)
fc = 1000;          % Carrier frequency (Hz)
T = 0.1;            % Total simulation time (s)
N = fs * T;         % Number of samples
t = (0:N-1)/fs;     % Time vector
bits = 10;          % Number of bits to transmit

% Generate random bit stream
data = randi([0 1], 1, bits);
samples_per_bit = round(fs / (bits / T));

%% 1. Amplitude Shift Keying (ASK)
figure('Name', 'ASK Modulation', 'Position', [100 100 800 600]);
nrz = repelem(data, samples_per_bit);
ask_signal = nrz .* sin(2*pi*fc*t(1:length(nrz)));

subplot(3,1,1); plot(t(1:length(nrz)), nrz); title('Original Binary Data'); ylim([-0.2 1.2]); grid on;
subplot(3,1,2); plot(t, sin(2*pi*fc*t)); title('Carrier Signal'); grid on;
subplot(3,1,3); plot(t(1:length(ask_signal)), ask_signal); title('ASK Modulated Signal'); grid on;

%% 2. Frequency Shift Keying (FSK)
figure('Name', 'FSK Modulation', 'Position', [200 200 800 600]);
f0 = fc; f1 = 2*fc; fsk_signal = zeros(1, length(nrz));
for i = 1:length(data)
    idx = (i-1)*samples_per_bit+1 : min(i*samples_per_bit, length(t));
    fsk_signal(idx) = sin(2*pi*(f0 + data(i)*(f1-f0))*t(idx));
end

subplot(3,1,1); plot(t(1:length(nrz)), nrz); title('Original Binary Data'); ylim([-0.2 1.2]); grid on;
subplot(3,1,2); plot(t(1:length(nrz)), f0 + data(ceil((1:length(nrz))/samples_per_bit))*(f1-f0)); title('Frequency Variation'); grid on;
subplot(3,1,3); plot(t(1:length(fsk_signal)), fsk_signal); title('FSK Modulated Signal'); grid on;

%% 3. Binary Phase Shift Keying (BPSK)
figure('Name', 'BPSK Modulation', 'Position', [300 300 800 600]);
bpsk_signal = (2*nrz - 1) .* sin(2*pi*fc*t(1:length(nrz)));

subplot(3,1,1); plot(t(1:length(nrz)), nrz); title('Original Binary Data'); ylim([-1.2 1.2]); grid on;
subplot(3,1,2); plot(t(1:length(bpsk_signal)), (2*nrz - 1)); title('Bipolar NRZ Encoding'); ylim([-1.2 1.2]); grid on;
subplot(3,1,3); plot(t(1:length(bpsk_signal)), bpsk_signal); title('BPSK Modulated Signal'); grid on;

%% 4. Quadrature Phase Shift Keying (QPSK)
figure('Name', 'QPSK Modulation', 'Position', [400 400 800 600]);
qpsk_data = randi([0 3], 1, bits/2);
qpsk_signal = pskmod(qpsk_data, 4, pi/4);
scatterplot(qpsk_signal); title('QPSK Constellation Diagram'); grid on;

%% 5. Quadrature Amplitude Modulation (16-QAM)
figure('Name', '16-QAM Modulation', 'Position', [500 500 800 600]);
qam_data = randi([0 15], 1, bits/4);
qam_signal = qammod(qam_data, 16);
scatterplot(qam_signal); title('16-QAM Constellation Diagram'); grid on;