clear


% [data,Fs] = audioread('SDRSharp_20210204_232546Z_915675010Hz_IQ.wav');
[data,Fs] = audioread('SDRSharp_20210204_232546Z_915675010Hz_IQ.wav');
T = 1/Fs;             % Sampling period       
L = size(data,1);             % Length of signal
% t = (0:L-1)*T;        % Time vector
t = 0:T:(L-1)*T;
ourmat = zeros(size(data,1),2);
for i = 1:size(data,1)
    ourmat(i,1) = t(i);
    ourmat(i,2) = abs(data(i,1) + data(i,1)*(1i));
end


% sec1 = data[1:.1*Fs]

% spectrogram(data(1:Fs*10))

%%
[y,fs] = audioread('SDRSharp_20210204_232546Z_915675010Hz_IQ.wav');
y = y(1:int64(round(.9*fs)),:);
ydft = fft(y);
% I'll assume y has even length
ydft = ydft(1:length(y)/2+1);
% create a frequency vector
freq = 0:fs/length(y):fs/2;
% plot magnitude
subplot(211);
plot(freq,abs(ydft));
set(gca, 'XScale', 'log')
% plot phase
subplot(212);
plot(freq,unwrap(angle(ydft))); 
xlabel('Hz');