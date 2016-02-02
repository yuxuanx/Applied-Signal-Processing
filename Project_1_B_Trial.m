%% Version:
% Testing using simulated channel (simulate_audio_channel)
clc; clear all; close all
load('s_training.mat'); % 128 training symbols
%% Generate bits
N = 128;
b = randsrc(1,2*N,[-1 1]);
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
R =10; % Upsample factor
ofdm_interpo = upsample(ofdm_package,R);
NN = 2^14; % Number of frequency grid points
f = (0:NN-1)/NN; % Normolized frequency
figure;semilogy(f,abs(fft(ofdm_package,NN))) % Check transform
xlabel('relative frequency f/fs');
ylabel('Amplitude(dB)');
title('Original Signal');
hold on;semilogy(f,abs(fft(ofdm_interpo,NN)))
xlabel('relative frequency f/fs');
ylabel('Amplitude(dB)');
title('Signal after upsampling');
%% Lowpass filter
B = firpm(32,2*[0 0.5/R*0.9 0.5/R*1.6 1/2],[1 1 0 0]);
ofdmSignal = conv(ofdm_interpo,B);
ofdmSignal = ofdmSignal(1:length(ofdm_interpo));
figure;
semilogy(f,abs([fft(ofdmSignal,NN) fft(B.',NN)]) ) % Check transforms
legend('Interpolated after LP filtering','LP-filter')
xlabel('relative frequency f/fs');
xlabel('relative frequency f/fs');
ylabel('Amplitude(dB)');
%% Modulation
fs = 22050; % Sampling frequency
fc = 4000; % Carrier frequency
F = (0:NN-1)/NN*fs;
n = (0:length(ofdmSignal)-1)';
ofdmModu = ofdmSignal.*exp(1i*2*pi*fc/fs*n);
semilogy(F,abs([fft(ofdmSignal,NN) fft(ofdmModu,NN) ]) ) % Check transforms
legend('Interpolated','Modulated')
xlabel('Frequency (Hz)');
xlabel('relative frequency f/fs');
ylabel('Amplitude(dB)');
%% Make signal real
ofdmReal = real(ofdmModu);
semilogy(F,abs([fft(ofdmSignal,NN) fft(ofdmModu,NN) fft(ofdmReal,NN) ]) ) % Check transforms
legend('Interpolated','Modulated','Real and modulated')
xlabel('Frequency (Hz)');
xlabel('relative frequency f/fs');
ylabel('Amplitude(dB)');
%% Transmit through channel
% Start with scaling from -1 to 1
% maxAmp = max(abs(ofdmReal));
% ofdmReal = ofdmReal/maxAmp;
% ofdmReal = [ofdmReal zeros(length(ofdmReal),1) ]; % Use only left audio channel
% sound(ofdmReal,fs);
% audiowrite('trial.wav',ofdmReal,fs)
%% Recording
% recObj = audiorecorder(fs,8,1);
% recordblocking(recObj,5);
% signalRec = getaudiodata(recObj);
% signalRec = signalRec/max(signalRec);
% figure;plot(signalRec);
%% Make it through a simulated channel
signalRec = simulate_audio_channel(ofdmReal,0.05); % Second input is noise level
figure;plot(signalRec);
signalRec = signalRec/max(signalRec);
xlabel('Sample');
ylabel('Amplitude');
title('Signal through simulated channel');
%% Determine where the signal starts
k = 1;
while(abs(signalRec(k))^2<0.05) % Energy detection
    k = k+1;
end
signal = signalRec(k:k+length(ofdmReal)-1);
figure;plot(signal);xlabel('Sample');ylabel('Normalized Amplitude');
title('Signal Frame');
%% Demodulation
% signal = ofdmReal; % Used for testing only

fs = 22050;
fc = 4000;
F = (0:NN-1)/NN*fs;
n = (0:length(ofdmSignal)-1)';
signalDemo = signal.*exp(-1i*2*pi*fc/fs*n);
figure;
semilogy(F,abs([fft(signal,NN) fft(signalDemo,NN)]) ) % Check transforms
legend('Modulated','Demodulated')
xlabel('Frequency (Hz)');
ylabel('Amplitude(dB)');
%% Design a LP decimation filter
B = firpm(32,2*[0 0.5/R*0.9 0.5/R*1.6 1/2],[1 1 0 0]);
y = conv(signalDemo,B); % y is now the interpolated signal
y = y(length(B):end);
figure;
semilogy(f,abs([fft(signalDemo,NN) fft(y,NN) fft(B.',NN)]) ) % Check transforms
legend('Demodulated','after LP filtering','LP-filter')
xlabel('relative frequency f/fs');
ylabel('Amplitude(dB)');
%% Down-sampling
y = downsample(y,R);
figure;
semilogy(f,abs(fft(y,NN)) ) % Check transforms
xlabel('relative frequency f/fs');
ylabel('Amplitude (dB)');
title('Signal after downsampling');
%% Channel Estimation
y_training = y(1:length(y)/2); % The half first block is used for estimation
y_training = y_training(len_cp+1:end); % Remove cyclic prefix
r_training = fft(y_training); % OFDM decoding
H = r_training./s_training;
fx = eps:1/(length(H)-1):1;
figure;semilogy(fx,abs(H));
xlabel('relative frequency');
ylabel('Amplitude (dB)');
figure;plot(real(ifft(fftshift(H))));
xlabel('Sample');ylabel('Amplitude');title('Channel estimation in time domain');
%% Channel Equalization
y_info = y(length(y)/2+1:end);
y_info = y_info(len_cp+1:end); % Compensate for the FIR filter delay
r = fft(y_info);
semilogy(abs(y));
ss = sign(real(r.*conj(H)))+1j*sign(imag(r.*conj(H)));
%% Symbols2Bits
bb = zeros(1,2*N); % Bits received
for k=1:N
    bb(2*k-1) = real(ss(k));
    bb(2*k) = imag(ss(k));
end
%% Bits error calculating
diff=b-bb;
error=find(diff~=0);
errorrate=length(error)/(2*N);