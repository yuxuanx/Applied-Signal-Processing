% Version: Channel equalization using known H(k)
% Using a for loop to implement Monte Carlo method to run the program for
% 1000 times to get the average errorrate.
clc; clear all; close all
%% Generate bits
N = 128;
% b = randsrc(1,2*N,[-1 1]);
load('b.mat'); % Load the random data generated before to avoid generate a 
% new one eachtime, otherwise the errorrate we calculated can be very
% unstable when we simulate different synchronization errors.
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
h = Channel(1);
H = fft(h,N); % DTFT
%% Decide and Add Cyclic Prefix
beta = 1.1; % Parameter deciding how long the cyclic prefix
len_cp_c = ceil(length(h)*beta); % Length of cyclic prefix
cyclic_prefix_c = zz(end-len_cp_c+1:end); % Cyclic prefix
zz = [cyclic_prefix_c;zz]; % Add cyclic prefix to the front
y_len = length(zz)+length(h)-1;
%% Generate Noise Sequence
errorrate = 0;
for i=1:1000 % Monte Carlo method starts
sigma = 0; % Noise level
w = 1/sqrt(2)*sigma*(randn(y_len,1) + 1i*randn(y_len,1)); % AWGN Channel1
%% Generate Received Signal
y = conv(h,zz) + w;
y = y(1:length(zz)); % Remove convolution redundancy
y_rec = y(len_cp_c+1:end); % Remove cyclic prefix
% Note that when simulate synchronization error, only once a scenario,
% which means you should comment another one
    %% Simulate Negative Sync Error
    ne = 0; % Used for simulate sync error, ne should be greater than 0
    y_rec = y(len_cp_c+1-ne:end-ne); % Use some samples from cyclic prefix
    %% Simulate Positive Sync Error
    %po = 0; % Used for simulate sync error, po should be greater than 0
    %y_rec = [y(len_cp_c+1+po:end);cyclic_prefix_c(1:po)]; % Use samples from neighbor frame
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
errorrate=length(error)/(2*N) + errorrate;
end
errorrate = errorrate/1000;
%% Some Plots
figure;
plot(real(r));hold on;
plot(real(ss));hold on;
plot(real(s));hold on;
xlabel('k');
ylabel('amplitude');
title('Comparison between r(k),real s(k) and estimated s(k)');
legend('r(k)','estimated s(k)','real s(k)');