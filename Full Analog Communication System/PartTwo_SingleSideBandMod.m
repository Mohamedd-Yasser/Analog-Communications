close all;
clc;
%==========================================================================
%--------------------------------Experiment 2 -----------------------------
%------------------------------------step1---------------------------------
[signal,Fs] = audioread('eric.wav');
L=length(signal);
dt = 1/Fs;
t = 0:dt:(L*dt)-dt;
F_signal = fftshift(fft(signal))/L;
mag = abs(F_signal);
f = linspace(-Fs/2 , Fs/2 , L);

figure('Name','Original Signal','NumberTitle','off');
subplot(2,1,1);
plot(t,signal);
title('Signal in time domain');

subplot(2,1,2);
plot(f,mag);
title('Magnitude spectrum');

%==========================================================================
%------------------------------------step2 ,3 &4 ---------------------------
%Used Low pass filter
idx_p = find(f >= 4000,1);
idx_n = find(f >= -4000, 1); 
filtered_signal = zeros(1,L);
filtered_signal(idx_n:idx_p) = F_signal(idx_n:idx_p);


mag_filtered = abs(filtered_signal);
figure('Name','Filtered Signal','NumberTitle','off');
subplot(2,1,1);
plot(f,mag_filtered);
title('Magnitude of Filtered Signal in Frequancy domain');
L=length(filtered_signal);
filtered_signal = ifft(ifftshift(filtered_signal))*L;
subplot(2,1,2);

plot(t,filtered_signal);
title('Filtered Signal in time domain');

sound(filtered_signal,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 


%==========================================================================
%----------------------------DSB-TC and DSB-SC-----------------------------
%DSB-Tc:

fc = 100000;
fs = 5*fc;
filtered_signal_after_resampling = resample(filtered_signal,fs,Fs);

dt = 1/fs;
L = length(filtered_signal_after_resampling);
t = 0:dt:(L*dt)-dt;

maxSignal= max(filtered_signal_after_resampling);
carrier = cos(2*pi*fc*t);
DSB_TC = carrier .*((2*maxSignal) + filtered_signal_after_resampling );

DSB_TC = fftshift(fft(DSB_TC))/L;
mag = abs(DSB_TC);
f = linspace(-fs/2 , fs/2 , L);

figure('Name','DSB Signals','NumberTitle','off');
subplot(2,1,1);
plot(f,mag); 
title('DSB-TC');

%DSB-SC:

DSB_SC = carrier .* filtered_signal_after_resampling;
DSB_SC = fftshift(fft(DSB_SC))/L;
mag = abs(DSB_SC);

subplot(2,1,2);
plot(f,mag); 
title('DSB-SC');

%==========================================================================
%------------------------------------step5 --------------------------------
%Used Low pass filter
idx_p = find(f >= 100000,1);
idx_n = find(f >= -100000, 1); 
SSB = zeros(1,L);
SSB(idx_n:idx_p) = DSB_SC(idx_n:idx_p);
mag_filtered1 = abs(SSB);

figure('Name','Lower Side Band SC','NumberTitle','off');
subplot(3,1,1);
plot(f,mag_filtered1);
title('Lower Side Band SC using Ideal Filter');
L=length(SSB);
SSB = ifft(ifftshift(SSB))*L;

%==========================================================================
%------------------------------------step6 --------------------------------
coherent_detector= carrier .*SSB;


%down-sampling the message as a part of the carrier removal:
coherent_detector_withoutLPF = resample(coherent_detector,Fs,fs);
L = length(coherent_detector_withoutLPF);
coherent_detector_withoutLPF_F = fftshift(fft(coherent_detector_withoutLPF))/L;

% %Low pass filter centered at zero or at the new received message center, returning to the original Fs
L=length(coherent_detector_withoutLPF);
f = linspace(-Fs/2 , Fs/2 , L); 
positive_index = find(f >= 4000,1);
negative_index = find(f >= -4000, 1); 
coherent_filtered_F= zeros(1,L);
coherent_filtered_F(negative_index:positive_index) = coherent_detector_withoutLPF_F(negative_index:positive_index);

coherent_filtered = real(ifft(ifftshift(coherent_filtered_F)))*L;

%gain=2, because we took at lec that the received power of message is 0.5*m(t)
coherent_SSB = coherent_filtered.*2;

sound(coherent_SSB,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 

%DSB_SC
Len=length(coherent_SSB);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

subplot(3,1,2);
plot(t,coherent_SSB);
title('LSB-SC using coherent detector in time domain');
grid on


%freq_domain
coherent_SSB = fftshift(fft(coherent_SSB))/L;
mag = abs(coherent_SSB);
f = linspace(-Fs/2 , Fs/2 , L); 

subplot(3,1,3);
plot(f,mag); 
title('LSB-SC using coherent detector in freq domain');
grid on

%==========================================================================
%------------------------------------step7 --------------------------------
L=length(DSB_SC);
[A,B] = butter(4,[2*(fc-4000)/fs 2*(fc)/fs]);
DSB_SC = ifft(ifftshift(DSB_SC))*L;
SSB_Butter = filter(A,B,DSB_SC);
SSB_Butter = fftshift(fft(SSB_Butter))/L;
mag= abs(SSB_Butter);
L2=length(SSB_Butter);
f2 = linspace(-fs/2 , fs/2 , L2);

figure('Name','Lower Side Band SC Butterworth','NumberTitle','off');
subplot(3,1,1);
plot(f2,mag);
title('Lower Side Band using Butterworth filter');
SSB_Butter = ifft(ifftshift(SSB_Butter))*L;


%Coherent detection of LSB after butterworth
coherent_detector= carrier .*SSB_Butter;


%down-sampling the message as a part of the carrier removal:
coherent_detector_withoutLPF = resample(coherent_detector,Fs,fs);
L = length(coherent_detector_withoutLPF);
coherent_detector_withoutLPF_F = fftshift(fft(coherent_detector_withoutLPF))/L;

% %Low pass filter centered at zero or at the new received message center, returning to the original Fs
L=length(coherent_detector_withoutLPF);
f = linspace(-Fs/2 , Fs/2 , L); 
positive_index = find(f >= 4000,1);
negative_index = find(f >= -4000, 1); 
coherent_filtered_F= zeros(1,L);
coherent_filtered_F(negative_index:positive_index) = coherent_detector_withoutLPF_F(negative_index:positive_index);

coherent_filtered = real(ifft(ifftshift(coherent_filtered_F)))*L;

%gain=2, because we took at lec that the received power of message is 0.5*m(t)
coherent_SSB = coherent_filtered.*2;

sound(coherent_SSB,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 

%DSB_SC
Len=length(coherent_SSB);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

subplot(3,1,2);
plot(t,coherent_SSB);
title('LSB-SC using coherent detector in time domain');
grid on

%freq_domain
coherent_SSB = fftshift(fft(coherent_SSB))/L;
mag = abs(coherent_SSB);
f = linspace(-Fs/2 , Fs/2 , L); 

subplot(3,1,3);
plot(f,mag); title('LSB-SC using coherent detector in freq domain');
grid on


%==========================================================================
%------------------------------------step8 --------------------------------
%-----------------------------------SNR 0dB------------------------------
SNR0_SSB = awgn(SSB,0);
coherent_detector= carrier .*SNR0_SSB;


%down-sampling the message as a part of the carrier removal:
coherent_detector_withoutLPF = resample(coherent_detector,Fs,fs);
L = length(coherent_detector_withoutLPF);
coherent_detector_withoutLPF_F = fftshift(fft(coherent_detector_withoutLPF))/L;

% %Low pass filter centered at zero or at the new received message center, returning to the original Fs
L=length(coherent_detector_withoutLPF);
f = linspace(-Fs/2 , Fs/2 , L); 
positive_index = find(f >= 4000,1);
negative_index = find(f >= -4000, 1); 
coherent_filtered_F= zeros(1,L);
coherent_filtered_F(negative_index:positive_index) = coherent_detector_withoutLPF_F(negative_index:positive_index);

coherent_filtered = real(ifft(ifftshift(coherent_filtered_F)))*L;
%gain=2, because we took at lec that the received power of message is 0.5*m(t)
coherent_SSB = coherent_filtered.*2;

sound(coherent_SSB,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 

%DSB_SC
Len=length(coherent_SSB);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

figure('Name','SNR Graphs for SSB-SC','NumberTitle','off');
subplot(3,2,1);
plot(t,coherent_SSB);
title('SSB-SC with 0 dB SNR in time domain');
grid on


%freq_domain
coherent_SSB = fftshift(fft(coherent_SSB))/L;
mag = abs(coherent_SSB);
f = linspace(-Fs/2 , Fs/2 , L); 

subplot(3,2,2);
plot(f,mag); title('SSB-SC with 0 dB SNR in freq domain');
grid on

%-----------------------------------SNR 10dB------------------------------
SNR10_SSB = awgn(SSB,10);
coherent_detector= carrier .*SNR10_SSB;


%down-sampling the message as a part of the carrier removal:
coherent_detector_withoutLPF = resample(coherent_detector,Fs,fs);
L = length(coherent_detector_withoutLPF);
coherent_detector_withoutLPF_F = fftshift(fft(coherent_detector_withoutLPF))/L;

% %Low pass filter centered at zero or at the new received message center, returning to the original Fs
L=length(coherent_detector_withoutLPF);
f = linspace(-Fs/2 , Fs/2 , L); 
positive_index = find(f >= 4000,1);
negative_index = find(f >= -4000, 1); 
coherent_filtered_F= zeros(1,L);
coherent_filtered_F(negative_index:positive_index) = coherent_detector_withoutLPF_F(negative_index:positive_index);

coherent_filtered = real(ifft(ifftshift(coherent_filtered_F)))*L;


%gain=2, because we took at lec that the received power of message is 0.5*m(t)
coherent_SSB = coherent_filtered.*2;

sound(coherent_SSB,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 

%DSB_SC
Len=length(coherent_SSB);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

subplot(3,2,3);
plot(t,coherent_SSB);
title('SSB-SC with 10 dB SNR in time domain');
grid on

%freq_domain
coherent_SSB = fftshift(fft(coherent_SSB))/L;
mag = abs(coherent_SSB);
f = linspace(-Fs/2 , Fs/2 , L); 

subplot(3,2,4);
plot(f,mag); title('SSB-SC with 10 dB SNR in freq domain');
grid on



%-----------------------------------SNR 30dB------------------------------
SNR30_SSB = awgn(SSB,30);
coherent_detector= carrier .*SNR30_SSB;

%down-sampling the message as a part of the carrier removal:
coherent_detector_withoutLPF = resample(coherent_detector,Fs,fs);
L = length(coherent_detector_withoutLPF);
coherent_detector_withoutLPF_F = fftshift(fft(coherent_detector_withoutLPF))/L;

% %Low pass filter centered at zero or at the new received message center, returning to the original Fs
L=length(coherent_detector_withoutLPF);
f = linspace(-Fs/2 , Fs/2 , L); 
positive_index = find(f >= 4000,1);
negative_index = find(f >= -4000, 1); 
coherent_filtered_F= zeros(1,L);
coherent_filtered_F(negative_index:positive_index) = coherent_detector_withoutLPF_F(negative_index:positive_index);

coherent_filtered = real(ifft(ifftshift(coherent_filtered_F)))*L;

%gain=2, because we took at lec that the received power of message is 0.5*m(t)
coherent_SSB = coherent_filtered.*2;

sound(coherent_SSB,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 

%DSB_SC
Len=length(coherent_SSB);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

subplot(3,2,5);
plot(t,coherent_SSB);
title('SSB-SC with 30 dB SNR in time domain');
grid on


%freq_domain
coherent_SSB = fftshift(fft(coherent_SSB))/L;
mag = abs(coherent_SSB);
f = linspace(-Fs/2 , Fs/2 , L); 

subplot(3,2,6);
plot(f,mag); title('SSB-SC with 30 dB SNR in freq domain');
grid on


%==========================================================================
%------------------------------------step9 --------------------------------
%Used Low pass filter
L=length(DSB_TC);
f = linspace(-fs/2 , fs/2 , L);
idx_p = find(f >= 100000,1);
idx_n = find(f >= -100000, 1); 
SSB_TC = zeros(1,L);
SSB_TC(idx_n:idx_p) = DSB_TC(idx_n:idx_p);
mag_filtered1 = abs(SSB_TC);

figure('Name','SSB-TC using ideal filter','NumberTitle','off');
subplot(2,1,1)
plot(f,mag_filtered1);
title('Lower Side Band TC using Ideal Filter');
SSB_TC = ifft(ifftshift(SSB_TC))*L;


f = linspace(-fs/2 , fs/2 , L);
envelope_TC=abs(hilbert(SSB_TC));
envelope_TC_F = fftshift(fft(envelope_TC))/L;
%DC blocking using high pass filter
idx_p = find(f >0, 1);%find 1st sample that is greater than zero 
envelope_TC_F(idx_p-1) =0;
%--------plotting & playing the sound after the envelop receiver 
%backing to the original sampling freq before plotting
%DSB_TC
envelope_TC = real(ifft(ifftshift(envelope_TC_F)))*L;
envelope_TC = resample(envelope_TC,Fs,fs);

Len_TC=length(envelope_TC);
dt = 1/Fs;
t = 0:dt:(Len_TC*dt)-dt;

subplot(2,1,2);
plot(t,envelope_TC);
axis([0.1 8.5 0.1 0.6]);
title('SSB-TC using envelope detector');
grid on

sound(envelope_TC,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 











