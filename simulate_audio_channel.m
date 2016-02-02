function [ yrec ] = simulate_audio_channel( zupmr, sigma )
%SIMULATE_AUDIO_CHANNEL 
% 
% Simulates tha audio media including a channel model, additive white noise and 
% a random delay respresenting the asyncronous timing between 
% receiver and transmitter.
%
% zupmr = signal to transmit over channel sampled at a rate of 22050 Hz
% sigma = standar deviation of additive noise at the receiver.
%         if sigma is not supplied or left empty the noise level is set to 
%         zero.
% 
% yrec = simulated received signal with a length corresponding to a
% duration of 5 s (assuming fs = 22050
%
% 20131116  Tomas McKelvey
%
%%

fs = 22050;
f0 = 4000;

if nargin<2 | isempty(sigma),
    sigma=0;
end

if ~isreal(zupmr);
    error('Input signal cannot be complex');
end


[nr,nc]=size(zupmr);
if nr<nc,
    zupmr = zupmr';
    nr = nc;
end

maxz=max(abs(zupmr));
zupmr = zupmr/maxz;  %Limit magnitude to 1.

zero_len= 6*22050-nr;
if mod(zero_len,2)
    % odd
    z1_len = (zero_len-1)/2;
    z2_len = (zero_len-1)/2+1;
else
    % even
    z1_len = zero_len/2;
    z2_len = zero_len/2;
end

% z1_len+z2_len + nr

zupmr_zp = [zeros(z1_len,1); zupmr; zeros(z2_len,1)];

z0 = 0.9*exp(j*2*pi*f0/fs);
a_chan = conv([ 1 -z0],[1 -conj(z0)]);
b_chan = conv([1 -1],[ 1 1]);

% [h,w] = freqz(b_chan,a_chan);
% plot(w/2/pi*fs,20*log10(abs(h)));

y = filter(b_chan,a_chan,zupmr_zp);



% Simulate asyncronuous reception
x = resample(y,4,1);
x = x(ceil(200*rand(1)):end);
yrec = resample(x,1,4);
yrec = yrec(1:5*22050)+ sigma*randn(5*22050,1); % Add noise

end

