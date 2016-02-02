% Version: Channel equalization using known H(k)
clc; clear all; close all
%% Generate bits
N = 128;
load('b.mat');
%% Bits2Symbols Using QPSK
M = 4; % Number of symbols in QPSK
m = log2(M); % Bits per Symbol
s_QPSK = [(1 + 1i) (1 - 1i) (-1 -1i) (-1 + 1i)]/sqrt(2); % QPSK Symbols
b_buffer = buffer(b, m)'; % Group bits into bits per symbol
s = zeros(N,1);
% Look up symbols using the indices
for k=1:N
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
zz = ifft(s);
%% Channel Description
h = Channel(2);
H = fft(h,N); % DTFT
%% Decide and Add Cyclic Prefix
beta = 1.1; % Parameter deciding how long the cyclic prefix
len_cp_c = ceil(length(h)*beta); % Length of cyclic prefix
cyclic_prefix_c = zz(end-len_cp_c+1:end); % Cyclic prefix
zz = [cyclic_prefix_c;zz]; % Add cyclic prefix to the front
y_len = length(zz)+length(h)-1;
%% Generate Noise Sequence
rate = zeros(1,1);
j = 1;
for sigma=0:0.001:0.05
    errorrate = 0;
for i=1:1000
w = 1/sqrt(2)*sigma*(randn(y_len,1) + 1i*randn(y_len,1)); % AWGN Channel1
%% Generate Received Signal
y = conv(h,zz) + w;
y = y(1:length(zz)); % Remove convolution redundancy
y_rec = y(len_cp_c+1:end); % Remove cyclic prefix
%% Simulate Negative Sync Error
ne = 0; % Used for simulate sync error, ne should be greater than 0
y_rec = y(len_cp_c+1-ne:end-ne); % Use some samples from cyclic prefix
%% OFDM Decoding 
r = fft(y_rec);
%% Channel Equalization With known H(k)
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
end
errorrate = errorrate/1000;
rate(j) = errorrate;
j = j+1;
end
%% Some Plots
figure;
stem(real(ss));hold on;
stem(real(s));
xlabel('k');
ylabel('amplitude');
title('Comparison between real s(k) and estimated s(k)');
legend('estimated s(k)','real s(k)');
figure;
plot(0:0.001:0.05,rate);
xlabel('noise level');
ylabel('errorrate');