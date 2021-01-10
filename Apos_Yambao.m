%Apos & Yambao ©2018
%EE15 Project
%Specifications: Band-stop filter at 142Hz<=f<=242Hz

close all; clear all; clc;
start();

function [db, W] = freqzM(num, den)
% script from Digital Signal Processing using MatLab, Ingle and Proakis
    [H, W] = freqz(num, den, 1000, 'whole');
    H = (H(1:1:501));
    W = (W(1:1:501));
    mag = abs(H); 
    db = 20*log10((mag+eps)/max(mag));
end

function start()
disp('142-242 Hz Band-stop Filter by Apos & Yambao');
disp('A sampling rate of 1000±40% Hz is recommended.');

prompt1 = 'Please enter the number of input frequencies: ';
prompt2 = 'Please enter input frequency: ';
prompt3 = 'Please enter the sampling rate: ';

N_F = input(prompt1); %user input for number of frequencies

%user input for frequency
F_I = [];
for k = 1:1:N_F
    F_I(k) = input(prompt2);
    k = k + 1;
end
%%%%%

s_r = input(prompt3); %user input for sampling rate
fs = s_r; %for use in spectrogram setup
S_R = 2*pi./s_r; %for use in zeroes and poles setup iteration 

%setting up of zeroes from 142Hz?f?242Hz
Z = [];
a = 1;
for j = 142:10:242
    Z(a) = exp(S_R*j*1i);
    a = a + 1;
    Z(a) = exp(S_R*j*-1i);
    a = a + 1;
end
%%%%%

%setting up of poles: 141.5Hz to 141.9Hz
P = [];
a = 1;
for j = 141.5:0.04:141.9
    P(a) = 0.58*exp(S_R*j*1i);
    a = a + 1;
    P(a) = 0.58*exp(S_R*j*-1i);
    a = a + 1;
end
%%%%%

%setting up of poles: 241.5Hz to 241.9Hz
for j= 241.5:0.04:241.9
    P(a)=0.6*exp(S_R*j*1i);
    a=a+1;
    P(a)=0.6*exp(S_R*j*-1i);
    a=a+1;
end
%%%%%

NUM = poly(Z); %extract coefficients of the characteristic polynomial of Z
DEN = poly(P);%extract coefficients of the characteristic polynomial of P

HZ = filt(NUM,DEN); %creates a discrete-time transfer function, HZ

[mag w] = freqzM(NUM, DEN); %calls freqzM function (defined); returns magnitude and radian frequency. For use in magnitude spectrum plot (3)

%waveform and spectrogram setup
syms t;
k = [0:fs-1];
x = 0;
for a=1:1:N_F
    x = x + cos(2*pi*F_I(a)*t);
end
xk = eval(subs(x, t, k/fs));
y = filter(NUM, DEN, xk);
%%%%

%1st window plots; (*number), *required plot
figure();
subplot(211);pzplot(HZ); %pole-zero plot
subplot(212);plot(w/pi, mag); grid on; title('Magnitude Spectrum') %magnitude spectrum of the designed filter (3)
xlabel('frequency (Hz)')
ylabel('Magnitude (dB)')

%1st window mod
set(gcf, 'Toolbar', 'none', 'Menu', 'none');
set(gcf, 'Name', 'Apos & Yambao Bandstop Project', 'NumberTitle', 'Off')
%%%%%

%2nd window plots
figure();
subplot(221); plot(xk); title('Input Signal Waveform') %input signal waveform (1)
subplot(222); spectrogram(xk); title('Input Signal Frequency Components') %input signal frequency components (2)
subplot(223); plot(y);title('Filtered Signal Waveform') %filtered signal waveform (4)
subplot(224); spectrogram(y); title('Filtered Signal Frequency Components') %filtered signal frequency components (5)

%2nd window mod
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
set(gcf, 'Toolbar', 'none', 'Menu', 'none');
set(gcf, 'Name', 'Apos & Yambao Bandstop Project', 'NumberTitle', 'Off')
%%%%%

end
