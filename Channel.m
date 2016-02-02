function [ h ] = Channel( channelIndex )
% Channel impulse response choosen
if channelIndex == 1
    h = zeros(60,1); % Channel 1
    for n=1:60
        h(n) = 0.8^(n-1);
    end
elseif channelIndex == 2
    h = [0.5;zeros(7,1);0.5]; %Channel 2
end

end

