%time series analysis for chaotic time-series from AFM
clf;
clear all
AFM=load('myharmonicseries.dat');
tmin=0;
tmax=0.005;
mystep=1*10^-6;
Time=[mystep:mystep:tmax]';
figure(1)
scatter(Time,AFM,8)
xlabel('Time')
ylabel('AFM_measurements')

%%

%simple descriptive methods

figure(2)
hist(AFM)


%%
%Differenced time series first order (periodicity left,last point left)

n=length(AFM); %number of data
for i=1:n-1
    yAFM(i)= AFM(i+1) - AFM(i);
end
TimeNew=Time(1:n-1);
figure(3)
plot(Time(1:n-1),yAFM,'r')

%%
%Differenced time series second order (noise left,last two points left)

for i=1:n-2
ynoise(i)= AFM(i+2) -2*AFM(i+1)+ AFM(i);
end

figure(4)

plot(Time(1:n-2),ynoise,'r')

%%
%Correlation coefficient matrix

figure(5)

autocorr(AFM,50,[],[])
title('ACF for harmonic TS')

%%
%Fourier Analysis

Y=fft(AFM);
figure(6)
plot(Y,'ro')
title('Fourier Coefficients in the Complex Plane');
xlabel('Real Axis');
ylabel('Imaginary Axis');
%%
m=length(Y);
Y(1)=[];
power = abs(Y(1:floor(m/2))).^2;
nyquist = 1/2;
freq = (1:m/2)/(m/2)*nyquist;
figure(7)
loglog(freq,power,'g')
xlabel('Fourier freq')
ylabel('PSD')
title('PSD of a periodic regime')
figure (8)
plot(freq,power,'r')
xlabel('Fourier freq')
ylabel('PSD')
title('periodogram of a periodic regime')
%%
%LAG PLOTS
figure(9)
clear i
clear j
j=1;AFMnew=zeros(length(AFM)-1,1);
for i=2:(length(AFM)-1)
    j=j+1;
AFMnew(j,1)=AFM(i-1,1);
end
scatter(AFM(2:end-1),AFMnew(2:end))
