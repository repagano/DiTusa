close all
clear objFab
clc

%% Simulated data: dHvA oscilation, frequencies F1 and F2, xi given by
%Lifshitz-Kosevich formula
% F1 = 3753;
% F2 = 3697;
H = 10:.00001:55;
% xi is proportional to H here, looking Schomberg P. 56 (P. 90 PDF)
F1 = 3761;%1;
F2 = 3692;
c = 3e-4;%2.2;
mom = .90;%
x = 7.;
a = 14.7*mom;
T = .96;
xi = (F1^2*exp(-a*x./H).*cos(2*pi*F1./H)./sinh(a*T./H)+F2^2*exp(-a*x./H).*cos(2*pi*F2./H+pi/4)./sinh(a*T./H)).*(c./(H.^(5/2)));%%3687%
% xi = 2.2*(cos(2*pi*F2./H)+
objFab.xUp = H;
objFab.yUp = xi;
objFab.temp = 1; 

%%
plot(objFab.xUp,objFab.yUp)
% figure
% plot(1./H,xi)
%% Fourier transform simulated data
endFields = 15:.5:55;%[35];%
iendFields = 1./endFields;
iFFspan = abs((1/40-1/55));%.02;%
fRange = [4 4500];
for jj = 1:length(endFields)
    iendField = 1/endFields(jj);
    iFFrange = [iendField+iFFspan, iendField];
    FFrange = 1./iFFrange     ;                  
    FFT = FFTcalc(objFab,'dHvA',FFrange,'up',fRange,'filter');      
    FFTmaxFab(jj) = max(FFT.FFT);
    FFTrangeFab(jj) = mean(FFT.range);
    Bm_fab(jj) = (1/2*(1./FFT.range(1)+1/FFT.range(2)))^(-1);
%     
%     figure
%     plot(FFT.f,FFT.FFT)
end
% objFab.FFT = FFTcalc(objFab,endFields,iFFspan,[5 4500]);%R
FFTfab.max = FFTmaxFab;
FFTfab.peak = FFTrangeFab;

%% Plot FFT max of simulated data and actual data
% Aexper = 19000000*exp(-a*x./CFmax1)./sinh(a*T./CFmax1)./(CFmax1.^(5/2));
% figure
% plot(1./Bm,FFTmax)%FFTmax)%,'b')
% hold on
% plot(1./Bmmaxstar,FFTmaxstar,'*')
% plot(1./Bm1,FFTmax1,'Color','m','LineWidth',.7)%[0.4940, 0.1840, 0.5560])[0 .7 0]
% plot(1./Bmmaxstar1,FFTmaxstar1,'*')
% plot(1./Bm_fab,FFTfab.max,'r')%*.025/.73
% title('Beating between 3760 and 3690')
% legend('FFT','FFT select pnts','FFT with bp filter 3500 to 3900 T','bp FFT select pnts','simulated FFT')
% xlabel('1/B_m (1/T)')
% ylabel('Fourier Amp. (Arb. Units)')
% ylim([0 1700])




%% Use class to load data 
loc = 'C:\Users\repag\Documents\DiTusa\MagLab\031919_PdGa\031919\';
temps = [1.52 .62 5.71 3.96 8 10 12 14 18 3.27 2.24 .96];
files = dir(strcat(loc,'*.ASC'));
filesFFT = dir(strcat(loc,'FT*'));
[temps,I] = find(temps == .96);
files = files(I);
filesFFT = filesFFT(I);
%% low pass filter with cut off frequency at 40000
obj = dHvA(loc,temps,files,filesFFT);
endFields = 15:.5:55;%
iendFields = 1./endFields;
iFFspan = abs((1/40-1/55));

dHvA.FFTload(obj,endFields,iFFspan,[5 6000],'lowpass',40000);

%find max value
for ii = 1:length(obj.FFT.range)
    [~,Itop] = min(abs(obj.FFT.range(1).upTemp.f-3500));
    [FFTmax(ii), I(ii)] = max(obj.FFT.range(ii).upTemp.FFT(Itop:end));
    CFmax1(ii) = mean(obj.FFT.range(ii).upTemp.range);
    Bm(ii) = (1/2*(1/obj.FFT.range(ii).upTemp.range(1)+1/obj.FFT.range(ii).upTemp.range(2)))^(-1)
    I(ii) = I(ii)+Itop;
end


%% Plot FFT of data with low pass filter
FFTmaxstar = FFTmax(1:8:end);
Bmmaxstar = Bm(1:8:end);
CFmaxstar = CFmax1(1:8:length(FFTmax));
figure
leg = [];
plt_cnt = 1;
cnt = 1;
c = hot(length(1:8:length(FFTmax))+6);
for ii = 1:8:length(FFTmax)
    p(plt_cnt) = plot(obj.FFT.range(ii).upTemp.f,obj.FFT.range(ii).upTemp.FFT,...
        'Color',c(cnt,:),'LineWidth',1)%(Itop:end)
    hold on
    p(plt_cnt+1) = plot(obj.FFT.range(ii).upTemp.f(I(ii)),FFTmax(ii),...
        '*','Color',p(plt_cnt).Color,'LineWidth',1)
    leg = [leg,{sprintf('%0.0f mT^{-1}',(1/Bm(ii))*1e3)}];
%     cnt = cnt+1;
    plt_cnt = plt_cnt+2;
    cnt = cnt+1;
%     pause(.5)
end
ledg = legend(p(1:2:end),leg,'Location','East','FontSize',13)
title(ledg,'B_m^{-1} values')
xlabel('Frequency (T)')
ylabel('Fourier amp. (arb. units)')
xlim([0 5500])
ylim([0 500])
ax = gca;
ax.XTick=500:1000:4500;
ax.YTick=0:100:400;
ax.FontSize = 13;
str = 'PdGa, 0.96 K, Fourier Window 6.8 mT^{-1}'
% str = [str newline 'Mulitple B_m values'];
dim = [.13 .80 .55 .12];
dim2 = [.15 .75 .55 .12];
t = annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',13,'LineStyle','none');
% t2 = annotation('textbox',dim2,'String',str2,'FontSize',12,'Boxoff');%'FitBoxToText','on',
%% fourier transform data with band pass filter with frequencies at 3500 
%and 3900
obj1 = dHvA(loc,temps,files,filesFFT);
endFields = 15:.5:55;%
iendFields = 1./endFields;
iFFspan = abs((1/40-1/55));

dHvA.FFTload(obj1,endFields,iFFspan,[5 4500],'bandpass', [3500 3900]);
%find max value
for ii = 1:length(obj1.FFT.range)
    [FFTmax1(ii), I(ii)] = max(obj1.FFT.range(ii).upTemp.FFT);
    CFmax(ii) = mean(obj1.FFT.range(ii).upTemp.range);
    Bm1(ii) = (1/2*(1/obj1.FFT.range(ii).upTemp.range(1)+1/obj1.FFT.range(ii).upTemp.range(2)))^(-1)
end

%% Plot FFT of data with low pass filter
FFTmaxstar1 = FFTmax1(1:8:end);
Bmmaxstar1 = Bm1(1:8:end);
CFmaxstar = CFmax1(1:8:length(FFTmax1));
figure
cnt = 1;
for ii = 1:8:length(FFTmax1)
    plot(obj1.FFT.range(ii).upTemp.f,obj1.FFT.range(ii).upTemp.FFT)%(Itop:end)
    hold on
    plot(obj1.FFT.range(ii).upTemp.f(I(ii)),FFTmax1(ii),'*')
%     cnt = cnt+1;
end
%% Plot FFT max of simulated data and actual data
% Aexper = 19000000*exp(-a*x./CFmax1)./sinh(a*T./CFmax1)./(CFmax1.^(5/2));
figure
plot(1./Bm,FFTmax)%FFTmax)%,'b')
hold on
plot(1./Bmmaxstar,FFTmaxstar,'*')
plot(1./Bm1,FFTmax1,'Color','m','LineWidth',.7)%[0.4940, 0.1840, 0.5560])[0 .7 0]
% plot(1./Bmmaxstar1,FFTmaxstar1,'*')
plot(1./Bm_fab,FFTfab.max,'r')%
% title('Beating between 3760 and 3690')
legend('FFT','FFT select pnts','FFT with bp filter 3500 to 3900 T','simulated FFT','FontSize',12)
xlabel('1/B_m (1/T)')
ylabel('Fourier Amp. (Arb. Units)')
ylim([0 1900])
xlim([.02 .07])
ax = gca;
ax.XTick=.02:.01:.07;
ax.YTick=0:400:1600;
ax.FontSize = 12;
str = sprintf('PdGa, beeting between 3761 and 3692 T Fourier Amp.');
dim = [.15 .80 .55 .12];
t = annotation('textbox',dim,'String',str,'FitBoxToText','on','FontSize',12,'EdgeColor',[1 1 1]);
%%
% figure
% for ii = 1:2:length(obj.FFT.range)
% 
% plot(obj1.FFT.range(20).upTemp.f,obj1.FFT.range(20).upTemp.FFT)
% hold on
% pause(1)
% end



%%
% x = [-10:10];
% y = exp(1./x)./sinh(1./x);
% figure
% plot(x,y)
