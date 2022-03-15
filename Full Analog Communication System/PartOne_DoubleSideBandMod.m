%==========================================================================
%--------------------------------Experiment 1 -----------------------------
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
plot(t,signal);title('Signal in time domain');
subplot(2,1,2);
plot(f,mag);title('Magnitude spectrum');

%==========================================================================
%------------------------------------step2 ,3 &4 ---------------------------
%Ideal Low Pass Filter For illustration
idx_p = find(f >= 4000,1); %find 1st sample that is equal or greater than 4000
idx_n = find(f >= -4000, 1); %find 1st sample that is equal or greater than -4000
lpf = zeros(1,L);
lpf (idx_n:idx_p) = 1;


%Used Low pass filter
idx_p = find(f >= 4000,1);
idx_n = find(f >= -4000, 1); 
filtered_signal_F = zeros(1,L);
filtered_signal_F(idx_n:idx_p) = F_signal(idx_n:idx_p);
mag_filtered = abs(filtered_signal_F);
filtered_signal = real(ifft(ifftshift(filtered_signal_F)))*L;

figure('Name','Filtered Signal','NumberTitle','off');
subplot(2,1,1);
plot(t,filtered_signal);title('Filtered Signal in time domain');
subplot(2,1,2);
plot(f,mag_filtered);title('Magnitude of Filtered Signal in Frequancy domain');

sound(filtered_signal,Fs);

%==========================================================================
%------------------------------------step5 --------------------------------
%DSB-Tc:

fc = 100000;
fs = 5*fc;
filtered_signal_after_resampling = resample(filtered_signal,fs,Fs);

dt = 1/fs;
L = length(filtered_signal_after_resampling);
t = 0:dt:(L*dt)-dt;

[maxSignal, indexOfMax] = max(abs(filtered_signal_after_resampling));
carrier = cos(2*pi*fc*t);
DSB_TC = carrier .*((2*maxSignal) + filtered_signal_after_resampling );

DSB_TC_F = fftshift(fft(DSB_TC))/L;
mag = abs(DSB_TC_F);
f = linspace(-fs/2 , fs/2 , L);

figure('Name','DSBTC','NumberTitle','off');
subplot(2,1,1);
plot(t,DSB_TC); title('DSB_TC_time');
subplot(2,1,2);
plot(f,mag); title('DSB_TC_freq');

%DSB-SC:

DSB_SC = carrier .* filtered_signal_after_resampling;
DSB_SC_F = fftshift(fft(DSB_SC))/L;
mag = abs(DSB_SC_F);
figure('Name','DSBSC','NumberTitle','off');
subplot(2,1,1);
plot(t,DSB_SC); title('DSB_SC_time');
subplot(2,1,2);
plot(f,mag); title('DSB_SC_freq');

%==========================================================================
%------------------------------------step6 & 7-----------------------------
%--------Envelope using hilbert takes the signal in the time domain:
f = linspace(-fs/2 , fs/2 , L);
DSB_TC = carrier .*((2*maxSignal) + filtered_signal_after_resampling );

envelope_TC=abs(hilbert(DSB_TC));
envelope_TC_F = fftshift(fft(envelope_TC))/L;
DSB_SC = carrier .* filtered_signal_after_resampling;
envelope_SC=abs(hilbert(DSB_SC));
envelope_SC_F = fftshift(fft(envelope_SC))/L;
%DC blocking using high pass filter
idx_p = find(f >0, 1);%find 1st sample that is greater than zero 
%envelope_TC_F(idx_p-1) =0;
envelope_SC_F(idx_p-1) =0;
%--------plotting & playing the sound after the envelop receiver 
%backing to the original sampling freq before plotting
%DSB_TC
envelope_TC = real(ifft(ifftshift(envelope_TC_F)))*L;
envelope_TC = resample(envelope_TC,Fs,fs);

Len_TC=length(envelope_TC);
dt = 1/Fs;
t = 0:dt:(Len_TC*dt)-dt;

figure('Name','Received Signals by envelope detector','NumberTitle','off');
subplot(2,1,1);
plot(t,envelope_TC);
title('DSB-TC using envelope detector');
grid on

pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 
sound(envelope_TC,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 



%backing to the original sampling freq before plotting
%DSB_SC
envelope_SC = real(ifft(ifftshift(envelope_SC_F)))*L;
envelope_SC = resample(envelope_SC,Fs,fs);
Len_SC=length(envelope_SC);
dt = 1/Fs;
t = 0:dt:(Len_SC*dt)-dt;

subplot(2,1,2);
plot(t,envelope_SC);
title('DSB-SC using envelope detector');
grid on

sound(envelope_SC,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 

%%observation: the signal is distorted in the case of DSB_SC, so this envelope detector can be used with th DSB_TC demodulation type only 

%==========================================================================
%For DSB-SC---------------------- steps 8->10------------------------------


% step 8:
%--------------------------------SNR 0dB  (worst noise case because signal to Noise ratio =1:1 , as the signal power will almost equal to the noise power so the the Noise signal power is very high here) ):
SNR0_DSB_SC = awgn(DSB_SC,0);
coherent_detector= carrier .*SNR0_DSB_SC;


%down-sampling the message as a part of the carrier removal:
coherent_detector_withoutLPF = resample(coherent_detector,Fs,fs);
L = length(coherent_detector_withoutLPF);
coherent_detector_withoutLPF_F = fftshift(fft(coherent_detector_withoutLPF))/L;

%Low pass filter centered at zero or at the new received message center, returning to the original Fs
L=length(coherent_detector_withoutLPF);
f = linspace(-Fs/2 , Fs/2 , L); 
positive_index = find(f >= 4000,1);
negative_index = find(f >= -4000, 1); 
coherent_filtered_F= zeros(1,L);
coherent_filtered_F(negative_index:positive_index) = coherent_detector_withoutLPF_F(negative_index:positive_index);

coherent_filtered = real(ifft(ifftshift(coherent_filtered_F)))*L;
%gain=2, because we took at lec that the received power of message is 0.5*m(t)
coherent_DSBSC = coherent_filtered.*2;

sound(coherent_DSBSC,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 

%DSB_SC
Len=length(coherent_DSBSC);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

figure('Name','Received Signals by coherent detector','NumberTitle','off');
subplot(3,3,1);
plot(t,coherent_DSBSC);
title('DSB-SC using coherent detector in time domain');
grid on



%DSB_SC_without_LPF
Len=length(coherent_detector_withoutLPF);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

subplot(3,3,2);
plot(t,coherent_detector_withoutLPF);
title('DSB-SC using coherent detector without LPF in time domain');
grid on

%freq_domain
coherent_DSBSC = fftshift(fft(coherent_DSBSC))/L;
mag = abs(coherent_DSBSC);
f = linspace(-Fs/2 , Fs/2 , L); 

subplot(3,3,3);
plot(f,mag); title('DSB-SC using coherent detector in freq domain');
grid on
%-----------------------------------SNR 10dB------------------------------
SNR10_DSB_SC = awgn(DSB_SC,10);
coherent_detector= carrier .*SNR10_DSB_SC;


%down-sampling the message as a part of the carrier removal:
coherent_detector_withoutLPF = resample(coherent_detector,Fs,fs);
L = length(coherent_detector_withoutLPF);
coherent_detector_withoutLPF_F = fftshift(fft(coherent_detector_withoutLPF))/L;


%Low pass filter centered at zero or at the new received message center, returning to the original Fs
L=length(coherent_detector_withoutLPF);
f = linspace(-Fs/2 , Fs/2 , L); 
positive_index = find(f >= 4000,1);
negative_index = find(f >= -4000, 1); 
coherent_filtered_F= zeros(1,L);
coherent_filtered_F(negative_index:positive_index) = coherent_detector_withoutLPF_F(negative_index:positive_index);

coherent_filtered = real(ifft(ifftshift(coherent_filtered_F)))*L;

%gain=2, because we took at lec that the received power of message is 0.5*m(t)
coherent_DSBSC = coherent_filtered.*2;

sound(coherent_DSBSC,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 

%DSB_SC
Len=length(coherent_DSBSC);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

subplot(3,3,4);
plot(t,coherent_DSBSC);
title('DSB-SC with 10 dB noise in time domain');
grid on

%DSB_SC_coherent_detector_withoutLPF
Len=length(coherent_detector_withoutLPF);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

subplot(3,3,5);
plot(t,coherent_detector_withoutLPF);
title('DSB-SC with 10 dB noise and without LPF in time domain');
grid on

%freq_domain
coherent_DSBSC = fftshift(fft(coherent_DSBSC))/L;
mag = abs(coherent_DSBSC);
f = linspace(-Fs/2 , Fs/2 , L); 

subplot(3,3,6);
plot(f,mag); title('DSB-SC using coherent detector in freq domain');
grid on



%-----------------------------------SNR 30dB------------------------------
SNR30_DSB_SC = awgn(DSB_SC,30);
coherent_detector= carrier .*SNR30_DSB_SC;


%down-sampling the message as a part of the carrier removal:
coherent_detector_withoutLPF = resample(coherent_detector,Fs,fs);
L = length(coherent_detector_withoutLPF);
coherent_detector_withoutLPF_F = fftshift(fft(coherent_detector_withoutLPF))/L;

%Low pass filter centered at zero or at the new received message center, returning to the original Fs
L=length(coherent_detector_withoutLPF);
f = linspace(-Fs/2 , Fs/2 , L); 
positive_index = find(f >= 4000,1);
negative_index = find(f >= -4000, 1); 
coherent_filtered_F= zeros(1,L);
coherent_filtered_F(negative_index:positive_index) = coherent_detector_withoutLPF_F(negative_index:positive_index);

coherent_filtered = real(ifft(ifftshift(coherent_filtered_F)))*L;

%gain=2, because we took at lec that the received power of message is 0.5*m(t)
coherent_DSBSC = coherent_filtered.*2;

sound(coherent_DSBSC,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 

%DSB_SC
Len=length(coherent_DSBSC);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

subplot(3,3,7);
plot(t,coherent_DSBSC);
title('DSB-SC with 30 dB noise in time domain');
grid on

%DSB_SC_coherent_detector_withoutLPF
Len=length(coherent_detector_withoutLPF);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

subplot(3,3,8);
plot(t,coherent_detector_withoutLPF);
title('DSB-SC with 30 dB noise and without LPF in time domain');
grid on

%freq_domain
coherent_DSBSC = fftshift(fft(coherent_DSBSC))/L;
mag = abs(coherent_DSBSC);
f = linspace(-Fs/2 , Fs/2 , L); 

subplot(3,3,9);
plot(f,mag); title('DSB-SC using coherent detector in freq domain');
grid on

%==========================================================================
%--- step 9, freq shift:
%freq error=100.1-100= 0.1 kHz:
%this is the "beat effect" phenomenon as a freq shift occurs
errored_fc = 100100;
fs = 500000;
dt = 1/fs;
L = length(filtered_signal_after_resampling);
t = 0:dt:(L*dt)-dt;

carrier = cos(2*pi*errored_fc*t);
coherent_detector= carrier .*DSB_SC;
coherent_detector_F=fftshift(fft(coherent_detector))/L;

%Low pass filter centered at zero or at the new received message center, returning to the original Fs
L=length(coherent_detector);
f = linspace(-Fs/2 , Fs/2 , L); 
positive_index = find(f >= 4000,1);
negative_index = find(f >= -4000, 1); 
coherent_filtered_F= zeros(1,L);
coherent_filtered_F(negative_index:positive_index) = coherent_detector_F(negative_index:positive_index);

coherent_detector = real(ifft(ifftshift(coherent_filtered_F)))*L;

%gain=2, because we took at lec that the received power of message is 0.5*m(t)
coherent_DSBSC = coherent_detector.*2;
%down-sampling the message as a part of the carrier removal:
coherent_DSBSC = resample(coherent_DSBSC,Fs,fs);
L = length(coherent_DSBSC);
coherent_DSBSC_F = fftshift(fft(coherent_DSBSC))/L;

sound(coherent_DSBSC,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 

%DSB_SC
Len=length(coherent_DSBSC);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

figure('Name','Received Signals by coherent detector with freq error','NumberTitle','off');
subplot(2,1,1);
plot(t,coherent_DSBSC);
title('DSB-SC with freq-error detector in time domain');
grid on


%freq_domain
coherent_DSBSC = fftshift(fft(coherent_DSBSC))/L;
mag = abs(coherent_DSBSC);
f = linspace(-Fs/2 , Fs/2 , L); 

subplot(2,1,2);
plot(f,mag); title('DSB-SC with freq-error detector in freq domain');
grid on

%==========================================================================
%--- step 10, phase shift error 20 degree:
fc = 100000;
fs = 5*fc;
dt = 1/fs;
L = length(filtered_signal_after_resampling);
t = 0:dt:(L*dt)-dt;
f = linspace(-Fs/2 , Fs/2 , L);

carrier = cos(2*pi*fc*t + (  20*(pi/180)   )   );
coherent_detector= carrier .*DSB_SC;
coherent_detector_f = abs(fftshift(fft(coherent_detector))/L);

%down-sampling the message as a part of the carrier removal:
coherent_detector_withoutLPF = resample(coherent_detector,Fs,fs);
L = length(coherent_detector_withoutLPF);
coherent_detector_withoutLPF_F = fftshift(fft(coherent_detector_withoutLPF))/L;

%Low pass filter centered at zero or at the new received message center, returning to the original Fs
L=length(coherent_detector_withoutLPF);
f = linspace(-Fs/2 , Fs/2 , L); 
positive_index = find(f >= 4000,1);
negative_index = find(f >= -4000, 1); 
coherent_filtered_F= zeros(1,L);
coherent_filtered_F(negative_index:positive_index) = coherent_detector_withoutLPF_F(negative_index:positive_index);

coherent_filtered = real(ifft(ifftshift(coherent_filtered_F)))*L;
%gain=2, because we took at lec that the received power of message is 0.5*m(t)
coherent_DSBSC = coherent_filtered.*2;

sound(coherent_DSBSC,Fs);
pause(10); %pause 10 sec, in order NOT to overlap two sucssessive voices 

%DSB_SC
Len=length(coherent_DSBSC);
dt = 1/Fs;
t = 0:dt:(Len*dt)-dt;

figure('Name','Received Signals by coherent detector with phase error shift with LPF ','NumberTitle','off');
subplot(2,1,1);
plot(t,coherent_DSBSC);
title('DSB-SC with phase-error detector in time domain');
grid on


%freq_domain
coherent_DSBSC = fftshift(fft(coherent_DSBSC))/L;
mag = abs(coherent_DSBSC);
f = linspace(-Fs/2 , Fs/2 , L); 

subplot(2,1,2);
plot(f,mag); title('DSB-SC with phase-error detector in freq domain');
grid on








