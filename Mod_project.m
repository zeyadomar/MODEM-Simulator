%**************************************************************************




clear all
close all
clc


%--------------------- Modulating Signal ---------------------



%first Lets import the speech signal


[speech_signal,f_sampling]=audioread("welcome.wav");

t=[0:1/f_sampling:(length(speech_signal)-1)/f_sampling];

speech_signal=speech_signal(1:length(t),1).';


%testing the signal 

sound(speech_signal)


% time domain representation of the speech signal
title_="Modulating signal in time domain";
plot_figures(1,t,speech_signal,title_,"Time (sec)","amplitude",211)

% Freq. domain representation of the speech signal


Y=fft(speech_signal);

title_="Modulating signal in freq domain";
x=f_sampling*linspace(-1,1,length(t));
plot_figures(1,x,abs(Y),title_,"f (Hz)","amplitude",212)




%**************************************************************************





%------------------------ Modulated Signal ---------------------



% Lets generate the carrier signal 

Fc=2*f_sampling;  % Carrier frequency should be >> signal frequency

Ac=2;     % Carrier amplitude


Ka= 0.95;  % Modulator senestivity 





%------- Conventional AM ------

C=Ac*cos(2*pi*Fc*t);
S_t_AM= (1+Ka*speech_signal).*C;



figure(2)
subplot(2,1,1)
plot(S_t_AM)
xlabel("Time (sec)")
ylabel("amplitude")
title("AM Modulated signal in time domain")





Y=fft(S_t_AM);


title_="AM Modulated signal in frequency domain";
x=f_sampling*linspace(-1,1,length(t));
plot_figures(2,x,abs(Y),title_,"f (Hz)","amplitude",212)
ylim([-1 300])



%-------- DSB --------

S_t_DSB = speech_signal.*C;

figure(3)
subplot(2,1,1)
plot(S_t_DSB)
title('DSB Modulated signal in time domain')
xlabel('time (sec)')
ylabel('amplitude')



Y=fft(S_t_DSB);

title_="DSB Modulated signal in frequency domain";
x=f_sampling*linspace(-1,1,length(t));
plot_figures(3,x,abs(Y),title_,"f (Hz)","amplitude",212)


%--------- SSB ----------

S_t_SSB = speech_signal.*cos(2*pi*Fc*t)+hilbert(speech_signal).*sin(2*pi*Fc*t);

title_str = 'SSB Modulated signal in time domain';

plot_figures(4,t,S_t_SSB,title_str,'time (sec)','amplitude',211)


Y=fft(S_t_SSB);

title_str = 'SSB Modulated signal in frequency domain';
x=f_sampling*linspace(-1,1,length(t));
plot_figures(4,x,abs(Y),title_str,'f (Hz)','amplitude',212)




%**************************************************************************




%----------------------- Modulated Signal + Noise ------------------------


%-------Conventional AM + Noise ------

noisy_AM_20 = awgn(S_t_AM,20,'measured');
noisy_AM_10 = awgn(S_t_AM,-10,'measured');
Y1=fft(noisy_AM_20);
Y2=fft(noisy_AM_10);


title_str = "noisy AM modulated signal 20 DB in time domain";

plot_figures(5,t,noisy_AM_20,title_str,'time (sec)','amplitude',221)


title_str = "20 DB in freq domain";
x=f_sampling*linspace(-1,1,length(t));
plot_figures(5,x,abs(Y1),title_str,'freq (Hz)','amplitude',222)
ylim([-0.1 300])

title_str='noisy AM modulated signal -10 DB in time domain';
plot_figures(5,t,noisy_AM_10,title_str,'time (sec)','amplitude',223)
ylim([-0.1 50]);


title_str = "-10 DB in freq domain";
x=f_sampling*linspace(-1,1,length(t));
plot_figures(5,x,abs(Y2),title_str,'freq (Hz)','amplitude',224)
ylim([-0.1 2000]);




%-------DSB + Noise-------

noisy_DSB_20 = awgn(S_t_DSB,20,'measured');
noisy_DSB_10 = awgn(S_t_DSB,-10,'measured');
Y1=fft(noisy_DSB_20);
Y2=fft(noisy_DSB_10);

title_str='noisy DSB modulated signal 20 DB in time domain';
plot_figures(6,t,noisy_DSB_20,title_str,'time (sec)','amplitude',221)

title_str = "20 DB in freq domain";
x=f_sampling*linspace(-1,1,length(t));
plot_figures(6,x,abs(Y1),title_str,'freq (Hz)','amplitude',222)
ylim([-0.1 300])


title_str='noisy DSB modulated signal -10 DB in time domain';
plot_figures(6,t,noisy_DSB_10,title_str,'time (sec)','amplitude',223)

title_str = "-10 DB in freq domain";
x=f_sampling*linspace(-1,1,length(t));
plot_figures(6,x,abs(Y2),title_str,'freq (Hz)','amplitude',224)
ylim([-0.1 300])



%-------SSB + Noise -----



noisy_SSB_20 = awgn(S_t_SSB,20,'measured');
noisy_SSB_10 = awgn(S_t_SSB,-10,'measured');
Y1=fft(noisy_SSB_20);
Y2=fft(noisy_SSB_10);

title_str='noisy SSB modulated signal 20 DB in time domain';
plot_figures(7,t,noisy_SSB_20,title_str,'time (sec)','amplitude',221)

title_str = "20 DB in freq domain";
x=f_sampling*linspace(-1,1,length(t));
plot_figures(7,x,abs(Y1),title_str,'freq (Hz)','amplitude',222)
ylim([-0.1 200])

title_str='noisy SSB modulated signal -10 DB in time domain';
plot_figures(7,t,noisy_SSB_10,title_str,'time (sec)','amplitude',223)

title_str = "-10 DB in freq domain";
x=f_sampling*linspace(-1,1,length(t));
plot_figures(7,x,abs(Y2),title_str,'freq (Hz)','amplitude',224)
ylim([-0.1 200])




%**************************************************************************




%---------------------- Demodulated Signal + Noise ------------------------

%------------ Demodulated Conventional AM ---------


%s(t)=Ac*cos(2piFc*t)+Ac*Ka*m(t)*Cos(2piFc*t) --coherent demodulation--->
% s(t)-Ac*Cos(2piFc*t)= Ac*Ka*m(t)*Cos(2piFc*t)
%then we will need to multiply by cos(2piFc*t) same as DSB demodulation
% Ac*m(t)*cos(2piFc*t)*cos(2piFc*t)= 0.5*Ac*m(t)+0.5*Ac*cos(4piFc*t)
%then we will use low pass filter to extract the m(t)

demod_AM_20 = (noisy_AM_20-C).*cos(2*pi*Fc*t);

demod_AM_20 = lowpass(demod_AM_20,0.5); % Low Pass Filter


figure(8)
subplot(2,1,1)
plot(demod_AM_20*2/Ac)
title("demodulated AM in time domain at 20 DB")

subplot(2,1,2)
plot(abs(fft(demod_AM_20)))
title("freq domain at 20 DB")
ylim([-0.1 400])


%lets listen to the demodulated signal 

%sound(demod_AM_20)
% sounds good :) as the 20 DB SNR made the signal power more dominant than
% the noise



demod_AM_10 = (noisy_AM_10-C).*cos(2*pi*Fc*t);

demod_AM_10 = lowpass(demod_AM_10,0.5); % Low Pass Filter

figure(9)
subplot(2,1,1)
plot(demod_AM_10)
title("demodulated AM in time domain at 10 DB")


subplot(2,1,2)
plot(abs(fft(demod_AM_10)))
title("freq domain at 10 DB")
ylim([-0.1 1000])

%lets listen to the -10 DB demodulated signal 

%WARNING TURN DOWN THE VOLUME SO AS NOT TO HURT YOUR EAR
%sound(demod_AM_10)

% sounds aweful :( as the -10 DB SNR made the noise power more dominant than
% the signal


%------------ Demodulated DSB ---------


demod_DSB_20 = noisy_DSB_20.*cos(2*pi*Fc*t);

%s(t)=Ac*m(t)*Cos(2piFc*t) --coherent demodulation--->
%Ac*m(t)*cos(2piFc*t)*cos(2piFc*t)= 0.5*Ac*m(t)+0.5*Ac*cos(4piFc*t)
%then we will use low pass filter to extract the m(t)

demod_DSB_20 = lowpass(demod_DSB_20,0.5); % Low Pass Filter

figure(10)
subplot(2,1,1)
plot(demod_DSB_20)
title("demodulated DSB in time domain at 20 DB")

subplot(2,1,2)
plot(abs(fft(demod_DSB_20)))
title("freq domain at 20 DB")
ylim([-0.1 400])


%lets listen to the demodulated signal 

sound(demod_DSB_20)
% sounds good :) as the 20 DB SNR made the signal power more dominant than
% the noise



demod_DSB_10 = noisy_DSB_10.*cos(2*pi*Fc*t);

demod_DSB_10 = lowpass(demod_DSB_10,0.5); % Low Pass Filter

figure(11)
subplot(2,1,1)
plot(demod_DSB_10)
title("demodulated DSB in time domain at 10 DB")


subplot(2,1,2)
plot(abs(fft(demod_DSB_10)))
title("freq domain at 10 DB")
ylim([-0.1 400])

%lets listen to the -10 DB demodulated signal 

%WARNING TURN DOWN THE VOLUME SO AS NOT TO HURT YOUR EAR
%sound(demod_DSB_10)

% sounds aweful :( as the -10 DB SNR made the noise power more dominant than
% the signal


%------------ Demodulated SSB ---------

%--------------------------------------------------------------------------
