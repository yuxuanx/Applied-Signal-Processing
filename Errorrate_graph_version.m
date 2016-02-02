% Version: Channel equalization using unknown H(k)
% Using decision feedback for estimation, 4 training symbols distributed
% uniformly in the symbols transmitted.
clc; clear all; close all
%% Generate bits
N = 128;
load('b.mat');
load('rate_knownH.mat');
%% Add training Bits
b = [1,1,b(1:64),1,1,b(65:128),1,1,b(129:192),1,1,b(193:end)];
%% Bits2Symbols Using QPSK
M = 4; % Number of symbols in QPSK
m = log2(M); % Bits per Symbol
s_QPSK = [(1 + 1i) (1 - 1i) (-1 -1i) (-1 + 1i)]/sqrt(2); % QPSK Symbols
b_buffer = buffer(b, m)'; % Group bits into bits per symbol
s = zeros(N+4,1);
% Look up symbols using the indices
for k=1:N+4 % 4 extra training symbols
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
H = fft(h,N+4); % DTFT
%% Decide and Add Cyclic Prefix
beta = 1.1; % Parameter deciding how long the cyclic prefix
len_cp_c = ceil(length(h)*beta); % Length of cyclic prefix
cyclic_prefix_c = zz(end-len_cp_c+1:end); % Cyclic prefix
zz = [cyclic_prefix_c;zz]; % Add cyclic prefix to the front
y_len = length(zz)+length(h)-1;
%%
rate = zeros(1,1);
j = 1;
for sigma=0:0.001:0.05
    errorrate = 0;
for i=1:1000
%     sigma = 0.0145; % Noise level
    w = 1/sqrt(2)*sigma*(randn(y_len,1) + 1i*randn(y_len,1)); % AWGN Channel1
    y = conv(h,zz) + w;
    y = y(1:length(zz)); % Remove convolution redundancy
    y_rec = y(len_cp_c+1:end); % Remove cyclic prefix
    %% Simulate Negative Sync Error
    %ne = 0; % Used for simulate sync error, ne should be greater than 0
    %y_rec = y(len_cp_c+1-ne:end-ne); % Use some samples from cyclic prefix
    %% Simulate Positive Sync Error
    %po = 7; % Used for simulate sync error, po should be greater than 0
    %y_rec = [y(len_cp_c+1+po:end);cyclic_prefix_c(1:po)]; % Use samples from neighbor frame
    r = fft(y_rec);
    ss = zeros(N+4,1);
    HH = zeros(N+4,1); % Guessed H(k)
    trainingSymbol = 1 + 1i;

    ss(1) = trainingSymbol;
    HH(1)=r(1)/ss(1);
for k=1:32 % For loop used for decision feedback
    HH(k+1)=r(k)/ss(k);
    ss(k+1) = sign(real(r(k+1).*conj(HH(k+1))))+1j*sign(imag(r(k+1).*conj(HH(k+1))));
end

ss(34) = trainingSymbol;
HH(34)=r(34)/ss(34);
for k=34:65 % For loop used for decision feedback
    HH(k+1)=r(k)/ss(k);
    ss(k+1) = sign(real(r(k+1).*conj(HH(k+1))))+1j*sign(imag(r(k+1).*conj(HH(k+1))));
end

ss(67) = trainingSymbol;
HH(67)=r(67)/ss(67);
for k=67:98 % For loop used for decision feedback
    HH(k+1)=r(k)/ss(k);
    ss(k+1) = sign(real(r(k+1).*conj(HH(k+1))))+1j*sign(imag(r(k+1).*conj(HH(k+1))));
end

ss(100) = trainingSymbol;
HH(100)=r(100)/ss(100);
for k=100:131 % For loop used for decision feedback
    HH(k+1)=r(k)/ss(k);
    ss(k+1) = sign(real(r(k+1).*conj(HH(k+1))))+1j*sign(imag(r(k+1).*conj(HH(k+1))));
end
bb = zeros(1,2*(N+4)); % Bits received
for k=1:N+4
    bb(2*k-1) = real(ss(k));
    bb(2*k) = imag(ss(k));
end
diff=b-bb;
error=find(diff~=0);
errorrate=length(error)/(2*N) + errorrate;
end
errorrate = errorrate/1000;
rate(j) = errorrate;
j = j+1;
end
%% Some plots
figure;
plot(real(r));hold on;
plot(real(ss));hold on;
plot(real(s));hold on;
xlabel('k');
ylabel('amplitude');
title('Comparison between r(k),real s(k) and estimated s(k)');
legend('r(k)','estimated s(k)','real s(k)');
figure;
plot(real(H));hold on;
plot(real(HH));hold on;
xlabel('k');
ylabel('amplitude');
title('Comparison between real H(k) and estimated H(k)');
legend('real H(k)','estimated H(k)');
figure;
plot(0:0.001:0.05,rate);hold on
plot(0:0.001:0.05,rate2)
xlabel('noise level (dB)');
ylabel('errorrate');
legend('estimated H(k)','known H(k)');
title('Error rate versus noise level');