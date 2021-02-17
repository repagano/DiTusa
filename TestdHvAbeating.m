% close all
% clear
clc
%% view interference on 
% F1 = 3753;
% F2 = 3697;
H = 10:.00001:55;
% xi is proportional to H here, looking Schomberg P. 56 (P. 90 PDF)
F1 = 3761%1;
F2 = 3692;
c = 6e-3;%2.2;
mom = .9;
x = 6.4;
a = 14.7*mom;
T = .96;
xi = (F1^2*exp(-a*x./H).*cos(2*pi*F1./H)./sinh(a*T./H)+F2^2*exp(-a*x./H).*cos(2*pi*F2./H+pi/4)./sinh(a*T./H)).*(c./(H.^(5/2)));%%3687%
% xi = 2.2*(cos(2*pi*F2./H)+
objFab.xUp = H;
objFab.yUp = xi;
objFab.temp = 1; 

% figure
% plot(1./H,xi)
%%
endFields = 15:.5:55;%[35];%
icenterFields = 1./endFields;
iFFspan = abs((1/40-1/55));%.02;%
fRange = [4 4500];
for jj = 1:length(endFields)
    iendField = 1/endFields(jj);
    iFFrange = [iendField+iFFspan, iendField];
    FFrange = 1./iFFrange     ;                  
    FFT = FFTcalc(objFab,'dHvA',FFrange,'up',fRange,'filter');      
    FFTmaxFab(jj) = max(FFT.FFT);
    FFTrangeFab(jj) = mean(FFT.range);
end
% objFab.FFT = FFTcalc(objFab,endFields,iFFspan,[5 4500]);%R
FFTfab.max = FFTmaxFab;
FFTfab.peak = FFTrangeFab;
%%
% figure
% plot(1./objFab.xUp,objFab.yUp)
%%
% figure
% plot(1./FFTrangeFab,FFTmaxFab)

%% Load data and Plot max of hig pass filtered data
% close all
% clear
clc
%% Use class to load data 
loc = 'C:\Users\repag\Documents\DiTusa\MagLab\031919_PdGa\031919\';
temps = [1.52 .62 5.71 3.96 8 10 12 14 18 3.27 2.24 .96];
files = dir(strcat(loc,'*.ASC'));
filesFFT = dir(strcat(loc,'FT*'));
[temps,I] = find(temps == .96);
files = files(I);
filesFFT = filesFFT(I);
%%
obj1 = dHvA(loc,temps,files,filesFFT);
%%
endFields = 15:.5:55;%
icenterFields = 1./endFields;
iFFspan = abs((1/35-1/55));
% %% Check span
% icenterField = 1./endFields;
% for ii = 1:length(endFields) 
%     chii = ii
%     s = [icenterFields(ii) icenterField(ii)+iFFspan]
%     span = 1./s
% end

%% Find maxs and plot
dHvA.FFTload(obj1,endFields,iFFspan,[5 4500]);
%%
for ii = 1:length(obj1.FFT.range)
    [FFTmax1(ii), I] = max(obj1.FFT.range(ii).upTemp.FFT);
    CFmax1(ii) = mean(obj1.FFT.range(ii).upTemp.range);
end

%%
figure
plot(1./CFmax,FFTmax)%,'b')
hold on
plot(1./CFmax1,FFTmax1,'Color',[0.4940, 0.1840, 0.5560])
plot(1./FFTfab.peak,FFTfab.max*.025/.73,'r')
title('Beating between 3760 and 3690')
legend('FFT','FFT with bp filter 3500 to 3900 T','simulated FFT')
xlabel('1/B_m (1/T)')
ylabel('Fourier Amp. (Arb. Units)')
ylim([0 1700])
%%
% figure
% for ii = 1:2:length(obj.FFT.range)
% 
% plot(obj.FFT.range(ii).upTemp.f,obj.FFT.range(ii).upTemp.FFT)
% hold on
% pause(1)
% end