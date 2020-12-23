function [faxis, pow] = get_PowerSpec(signal, SamplingRate, DecibelsFlag ,plotFlag)
% For example for a 10sec segment sampled at millisec resolution use:
% get_PowerSpec(signal, SamplingRate, DecibelsFlag ,plotFlag)
signal = Testep;
SamplingRate = 1/250;

plotFlag = 1; 
DecibelsFlag = 1; 
% get power
pow = (abs(fft(signal)).^2)/length(signal);
% convert to decibels
if DecibelsFlag==1
    pow = 10*log10(pow/max(pow));
end
% first half of data without negative frequencies
pow = pow(1:min(floor(length(signal)/2)+1,length(signal)));
% define df and fNQ
df = 1/(length(signal)/SamplingRate);
fNQ = SamplingRate/2;

faxis = (0:df:fNQ);
if (plotFlag)
    figure;
    plot(faxis, pow);
    xlabel('Frequency [Hz]')
    ylabel('Power [dB]')
end
end

