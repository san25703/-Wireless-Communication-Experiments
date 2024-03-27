%Amplitude Modulation
clc;
clear all;
close all;
f = 1;
Fc = 10;
Fs = 100;
t = 0:1/Fs:5;
s = sin(2*pi*f*t); % Original signal
c = sin(2*pi*Fc*t);
subplot(4,1,1)
plot(t,s);
xlabel('time');
ylabel('amplitude')
title('original signal (U21EC080)');
subplot(4,1,2)
plot(t,c);
xlabel('time');
ylabel('amplitude');
title('carrier signal (U21EC080)');

y1 = ammod(s,Fc,Fs); % Modulate.
d1 = amdemod(y1,Fc,Fs); % Demodulate.
subplot(4,1,3);
plot(t,y1);
xlabel('time');
ylabel('amplitude');
title('am modulated (U21EC080)');
subplot(4,1,4);
plot(t,d1);
xlabel('time');
ylabel('amplitude');
title('am demodulated (U21EC080)');

figure;
%Frequency Modulation
t = [0:0.001:1]; %upto 1000 samples
f1=5;
m=cos(2*pi*f1*t);
subplot(2,2,1); %plotting message signal
plot(t,m);
xlabel('Time');
ylabel('Amplitude');
title('Message Signal (U21EC080)');
f2=100;
c=sin(2*pi*f2*t);
subplot(2,2,2); %plotting carrier signal
plot(t,c);
xlabel('Time');
ylabel('Amplitude');
title('Carrier Signal (U21EC080)');
mf=15;
y = sin((2*pi*f2*t)+(mf*sin(2*pi*f1*t)));
subplot(2,2,3);%plotting FM (Frequency Modulated) signal
plot(t,y);
xlabel('Time');
ylabel('Amplitude');
title('FM Signal (U21EC080)');
x= diff(y);
z=abs(x);
[b,a]=butter(10, 0.056);
a1=filter(b,a,z);
subplot(2,2,4);
plot(a1);
xlabel('Time');
ylabel('Amplitude');
title('Demodulated Signal (U21EC080)');
figure;
% Frequency Shift Keying

fc1=10 %1st Sine Wave carrier
fc2=30 %2nd Sine Wave carrier
fp=5 %Periodic Binary pulse (Message)
amp=4 %amplitude (For Both Carrier & Binary Pulse Message)
amp=amp/2;
t=0:0.001:1;
c1=amp.*sin(2*pi*fc1*t);%1st Carrier Sine wave
c2=amp.*sin(2*pi*fc2*t);%2nd Carrier Sine wave
subplot(5,1,1); %Carrier wave
plot(t,c1)
xlabel('Time')
ylabel('Amplitude')
title('Carrier 1 Wave (U21EC080)')
subplot(5,1,2) %Carrier wave
plot(t,c2)
xlabel('Time')
ylabel('Amplitude')
title('Carrier 2 Wave (U21EC080)')
m=amp.*square(2*pi*fp*t)+amp;%Generating Square wave message
% Adds amp to the scaled square wave, effectively shifting the entire waveform upwards.
% This is done to ensure that the signal is always positive.
% as the square function produces values in the range [-1, 1].
subplot(5,1,3) %Square Binary Pulse Message
plot(t,m)
xlabel('Time')
ylabel('Amplitude')
title('Binary Message Pulses (U21EC080)')
for i=0:1000 %generating the modulated wave
    if m(i+1)==0
        mm(i+1)=c1(i+1);
        else
        mm(i+1)=c2(i+1);
    end
end
subplot(5,1,4) %Plotting The Modulated wave
plot(t,mm)
xlabel('Time')
ylabel('Amplitude')
title('Modulated Wave (U21EC080)')




for i=0:1000
    if mm(i+1)==c1(i+1)
        dm(i+1)=0;
    else
       dm(i+1)=1;
    end
end
subplot(5,1,5);
plot(t,dm);
xlabel('Time')
ylabel('Amplitude')
title('Demodulated Wave (U21EC080)')
%Phase Shift Keying
figure;
t=0:.001:1;
fc=60 %frequency of Carrier Sine wave
fm=10 %Message frequency
amp=3 %Carrier & Message Amplitude
c=amp.*sin(2*pi*fc*t);% Carrier Sine wave
subplot(4,1,1)
plot(t,c)
xlabel('Time')
ylabel('Amplitude')
title('Carrier (U21EC080)')
m=square(2*pi*fm*t);%Message signal
subplot(4,1,2)
plot(t,m)
xlabel('time')
ylabel('amplitude')
title('Message Signal')
% Sine wave multiplied with square wave in order to generate PSK
x=c.*m;
subplot(4,1,3) %PSK (Phase Shift Keyed) signal
plot(t,x)
title('PSK Modulated Output');
xlabel('t')
ylabel('y')
de = pskdemod(x,1);
subplot(4,1,4);
plot(t,de)
title('PSK Demodulated Output');
xlabel('t')
ylabel('y')
