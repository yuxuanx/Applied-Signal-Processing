%% Version:
% Using cross correlation to automatically detect where signal start while
% receiving in real channel
clc; clear all; close all
load('s_training.mat'); % 128 training symbols
%% Generate bits
N = 128;
load('b.mat');
% b = randsrc(1,2*N,[-1 1]);
%% Bits2Symbols Using QPSK
M = 4; % Number of symbols in QPSK
m = log2(M); % Bits per Symbol
s_QPSK = [(1 + 1i) (1 - 1i) (-1 -1i) (-1 + 1i)]/sqrt(2); % QPSK Symbols
b_buffer = buffer(b, m)'; % Group bits into bits per symbol
s = zeros(N,1);
% Look up symbols using the indices
for k=1:N % 4 extra training symbols
    if b_buffer(k,:) == [1 1]
        s(k) = 1 + 1i;
    elseif b_buffer(k,:) == [1 -1]
        s(k) = 1 - 1i;
    elseif b_buffer(k,:) == [-1 -1]
        s(k) = -1 - 1i;
    else
        s(k) = -1 + 1i;
    end
end
%% Generate OFDM Seuqence
zz_training = ifft(s_training); % First block used for channel estimation
zz = ifft(s); % Second block for transmitting information
%% Add Cyclic Prefix
len_cp = 80; % Length of cyclic prefix
cyclic_prefix = zz(end-len_cp+1:end); % Cyclic prefix
zz = [cyclic_prefix;zz]; % Add cyclic prefix to the front
cyclic_prefix_training = zz_training(end-len_cp+1:end);
zz_training = [cyclic_prefix_training;zz_training];
ofdm_package = [zz_training;zz]; % Package for transmitting
%% Interpolation
R = 10; % Upsample factor
ofdm_interpo = upsample(ofdm_package,R);
NN = 2^14; % Number of frequency grid points
f = (0:NN-1)/NN; % Normolized frequency
figure;semilogy(f,abs(fft(ofdm_package,NN))) % Check transform
figure;semilogy(f,abs(fft(ofdm_interpo,NN)))
%% Lowpass filter
B = firpm(32,2*[0 0.5/R*0.9 0.5/R*1.6 1/2],[1 1 0 0]);
ofdmSignal = conv(ofdm_interpo,B);
ofdmSignal = ofdmSignal(1:length(ofdm_interpo));
figure;
semilogy(f,abs([fft(ofdmSignal,NN) fft(B.',NN)]) ) % Check transforms
legend('Interpolated after LP filtering','LP-filter')
xlabel('relative frequency f/fs');
%% Modulation
fs = 22050; % Sampling frequency
fc = 4000; % Carrier frequency
F = (0:NN-1)/NN*fs;
n = (0:length(ofdmSignal)-1)';
ofdmModu = ofdmSignal.*exp(1i*2*pi*fc/fs*n);
semilogy(F,abs([fft(ofdmSignal,NN) fft(ofdmModu,NN) ]) ) % Check transforms
legend('Interpolated','Modulated')
xlabel('Frequency (Hz)');
%% Make signal real
ofdmReal = real(ofdmModu);
semilogy(F,abs([fft(ofdmSignal,NN) fft(ofdmModu,NN) fft(ofdmReal,NN) ]) ) % Check transforms
legend('Interpolated','Modulated','Real and modulated')
xlabel('Frequency (Hz)');
%% Transmit through channel
% Start with scaling from -1 to 1
maxAmp = max(abs(ofdmReal));
ofdmReal = ofdmReal/maxAmp;
ofdmReal = [ofdmReal zeros(length(ofdmReal),1) ]; % Use only left audio channel
sound(ofdmReal,fs);
% audiowrite('trial.wav',ofdmReal,fs); % Use to test on our own lap