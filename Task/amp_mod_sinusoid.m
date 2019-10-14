% Amplitude modulated sinusoid

Ac=input('enter carrier signal amplitude');
Am=input('enter message signal amplitude');
fc=input('enter carrier frequency');
fm=input('enter message frequency');% fm<fc
m=input('enter modulation index');
t=input('enter time period');
t1=linspace(0,t,1000);
y1=sin(2*pi*fm*t1); % message signal
y2=sin(2*pi*fc*t1); % carrier signal
eq=(1+m.*y1).*(Ac.*y2);
subplot(311);
plot(t1,y1);
xlabel('Time');
ylabel('Amplitude');
title('Message signal')
subplot(312)
plot(t1,y2);
xlabel('Time');
ylabel('Amplitude');
title('Carrier signal');
subplot(313);
plot(t1,eq);
plot(t1,eq,'r');
xlabel('Time');
ylabel('Amplitude');
title('Modulated signal');