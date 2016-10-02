% Consecutive mass measuements of a hanging Silicon sphere
clf
clear all
case5=load('SI081208.txt');
Mass=case5(:,2);
Time=[0:20:20*(length(Mass)-1)];
figure(1)

set(gca,'FontSize',20,'FontName','Times')
scatter(Time,Mass,'k')
xlabel('Time (sec)','fontsize',25,'FontName','Times')
set(gca,'box','on')
ylabel('Mass (gr)','fontsize',25,'FontName','Times')
set(gca,'box','on')


%print -depsc figure1.eps
%%
%Methods for transform non-stationary to stationary ts
%Differenced time series first order (periodicity left,last point left)

n=length(Mass); %number of data
for i=1:n-1
    
    y(i) = Mass(i+1) - Mass(i);
   
end
TimeNew=Time(1:n-1);
figure(2)
plot(Time(1:n-1),y,'k')
title('First Differences');
xlabel('Time (sec)');
ylabel('First Diff');
%%

figure(3)
steps=round(n/4);
set(gca,'FontSize',20,'FontName','Times')

autocorr(y,steps,2,[])
set(findobj('Type','line'),'Color','k')
xlabel('Lag','fontsize',25,'FontName','Times')
grid off
ylabel('Sample ACF','fontsize',25,'FontName','Times')
grid off


%print -depsc -tiff figure3.eps
%%
%LAG PLOT for y
clf
figure(4)
clear i
clear j
ytrans=y';
j=1;ynew=zeros(length(y)-1,1);
for i=2:(length(y)-1)
    j=j+1;
    
ynew(j,1)=ytrans(i-1,1);
end
set(gca,'FontSize',20,'FontName','Times')
scatter(y(2:end-1),ynew(2:end),'k')

a=-0.7;b=4*10^-7;
ynew=a*y+b;
hold on
plot(y,ynew,'r')
set(gca,'XLim',[-1*10^-3 1*10^-3],'YLim',[-1*10^-3 1*10^-3])
xlabel('x_n','fontsize',25,'FontName','Times')
set(gca,'box','on')
ylabel('x_n_-_1','fontsize',25,'FontName','Times')
set(gca,'box','on')
hold off

%print -depsc -tiff figure4.eps
%%
%Fourier Analysis

Y1=fft(y);
figure(5)
plot(Y1,'ro')
title('Fourier Coefficients in the Complex Plane');
xlabel('Real Axis');
ylabel('Imaginary Axis');
%%
p=length(Y1);
Y1(1)=[];
power1 = abs(Y1(1:floor(p/2))).^2;

sampltime=20;
nyquist = 1/(2*sampltime);
freq1 = (1:p/2)./(p*sampltime);
%%

figure (6)
plot(freq1,power1,'k')
set(gca,'XLim',[0 0.025],'YLim',[0 1.5*10^-3])
set(gca,'FontSize',20,'FontName','Times')
xlabel('Frequency','fontsize',25,'FontName','Times')
ylabel('Power','fontsize',25,'FontName','Times')

%print -depsc -tiff figure6.eps
%%
figure (7)
loglog(freq1,power1,'k')
set(gca,'FontSize',20,'FontName','Times')
hold on
a1=3*10^-10;b1=1;power1=a1.*freq1.^(-b1);
plot(freq1(1:50),power1(1:50),'b')

a2=3*10^-7;b2=0;power1=a2.*freq1.^(-b2);
plot(freq1(47:150),power1(47:150),'r')

a3=10^-3;b3=-1.3;power1=a3.*freq1.^(-b3);
plot(freq1(150:end),power1(150:end),'g')

xlabel('Frequency','fontsize',25,'FontName','Times')
ylabel('Power','fontsize',25,'FontName','Times')
%print -depsc -tiff figure7.eps
hold off
%%
%Bretthorst Method
p=length(Y1);
Y1(1)=[];
power1 = abs(Y1(1:floor(p/2))).^2;

sampltime=20;
nyquist = 1/(2*sampltime);
freq1 = (1:p/2)./(p*sampltime);
period = 1./freq1;
dmean=mean(y(:).^2);
Comega=(1/p)*(real(Y1(1:floor(p/2))).^2+imag(Y1(1:floor(p/2))).^2);
Pomega=3*10^-21*(1-(2.*Comega/(p*dmean))).^((2-p)/2);
[mp,index2] = max(Pomega);
freq1(index2)

figure(8)
plot(freq1,Pomega,'r')
set(gca,'FontSize',20,'FontName','Times')
hold on
plot(freq1,power1,'k')
set(gca,'XLim',[0.01 0.025],'YLim',[0 1.5*10^-3])
xlabel('Frequency','fontsize',25,'FontName','Times')
ylabel('P(f|D,I)','fontsize',25,'FontName','Times')
hold off
print -depsc -tiff figure8.eps
%%
%Allan and Modified Allan
% Call Function Allan and define the 
% sampling period tau0

tau0=1;
[var2,varmod,time]=Allancalc(y,tau0);


%%

figure(9)
timenew=log2(time);
allannew=log2(varmod);
allannew2=log2(var2);

x=[15:55];x1=[90:1000];x2=[2:14];
loglog(time,varmod,'-sk')
%hold on
%loglog(time,var2,'-sk')
%line(x,10^-4.6.*x.^-0.51,'color','r','linewidth',2,'linestyle','--')
%line(x1,0.9*10^-5.7*x1.^-0,'color','g','linewidth',2,'linestyle','--')
line(x2,1.6*10^-4.1*x2.^-1.5,'color','b','linewidth',2,'linestyle','--')
set(gca,'FontSize',20,'FontName','Times')
xlabel('Time','fontsize',25,'FontName','Times')
ylabel('Mod Allan Dev','fontsize',25,'FontName','Times')
%hold off
print -depsc -tiff figure12.eps
