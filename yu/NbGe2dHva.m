close all
clear
clc

%% Here the files to be loaded are chosen.
addpath C:\Users\repag\Documents\MATLAB\DiTusa
loc = 'C:\Users\repag\Documents\DiTusa\yu NbGe2\data2\';
temps = [3.49 2.67 2.18 2.38 2.03 1.77 1.54 1.55 1.11 3.98 3.06];%Temperatures
% corresponding to files to be loaded
filesH = dir(strcat(loc,'*P001.dat'));
varTemps = temps;
I = ismember(temps,varTemps);
temps = temps(I);%
files = filesH(I);

%%
obj = dHvA(loc,temps,files);

%% Here the 1/H windows that are selected and fourier transformed are defined
clc
close all
field1 = 20;
endFields = 60;%the endFields are the maximum field values of each window

iFFspan = abs((1/field1-1/endFields));%Here the width of the 1/H window is defined 
% startFields = 15;
% endFields = 1/(-iFFspan*2+1/startFields)%% Here the Fourier transform is performed  

dHvA.FFTload(obj,endFields,iFFspan,[10 5e3],'extFFT');%7e4],'extFFT')%
% %%
% for ii = 1%1:length(obj.FFT.range.upTemp)
%     figure
%     plot(1./obj.FFT.range.upTemp(ii).xspl,obj.FFT.range.upTemp(ii).yspl)
%     xlabel('Frequency (T)')
%     ylabel('Frequency after bs')
%     title('p007_012221P001.dat')
% end

%%
peakRangeDownH =  [350, 450; 550, 670; 770, 865; 2140, 2220; 3470, 3560; 3800, 3930; 4340, 4440]; %HFT, 20 to 60 T, falling field% [4340, 4440];%
peakRangeUpH = [160, 240; 350, 460; 560, 620; 770, 875]; %HFT, 20 to 60 T, rising field
peakRangeDownP = [500, 620];

% 15to55% [2666, 2858; 3194, 3400; 3558, 3875];% 15to35%
% peakRange
% dHvA.massLoad runs the method that calculates the effective masses.
dHvA.massLoad(obj,peakRangeDownP);

%%
figure
plot(obj.mass.range.upPeak.temp,obj.mass.range.upPeak.AoT,'*',...
    obj.mass.range.upPeak.temp,obj.mass.range.upPeak.AoTfitV)
%%
tv = obj.FFT.range.upTemp(1).f';
leg = []
figure
for ii = 1:length(obj.FFT.range(end).upTemp)%4%
%     figure
%     for jj = 1:length(obj.mass.range.upPeak)
%         peakA(jj) = obj.mass.range.upPeak(jj).A(ii)
%         peakF(jj) = obj.mass.range.upPeak(jj).maxFreq(ii)
%         chjj = jj
%     end
    c = parula(length(obj.FFT.range.upTemp));
    plot(obj.FFT.range.upTemp(ii).f,obj.FFT.range.upTemp(ii).FFT,'LineWidth',1)%,'Color',c(ii,:)
    hold on
%     plot(peakF,peakA,'*')
    txt = {strcat(num2str(obj.FFT.range.upTemp(ii).temp),' K')};
    txtpeak = {strcat(num2str(obj.FFT.range.upTemp(ii).temp),' K peak')};
    leg = [leg,txt,txtpeak];
end
xlabel('Frequency (T)')
ylabel('Fourier amp. (arb units)')
legend(leg)
title('P sample, range of 40 to 60 T, falling field')
%%
for ii = 1:length(obj.mass.range.upPeak)
    
    figure
    subplot(1,2,1)
    plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).A,'*')
    subplot(1,2,2)
    plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).AoT,'*')
    hold on
    plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).AoTfitV)
    tv = mean(obj.mass.range.upPeak(ii).maxFreq);
    suptitle(sprintf('max peak,  %0.1f T', tv))
end
%%
% for ii = 2%1:length(obj.FFT.range.upTemp)
%     figure
% %     plot(1./obj1.FFT.range.upTemp(ii).xspl,obj1.FFT.range.upTemp(ii).yspl,'r','LineWidth',1.1)
% %     hold on
% %     s = 
%     plot(1./obj.FFT.range.upTemp(ii).xspl,obj.FFT.range.upTemp(ii).yspl,'r','LineWidth',1.2)
%     hold on
%     plot(1./obj1.FFT.range.upTemp(ii).xspl,obj1.FFT.range.upTemp(ii).yspl,'b','LineWidth',.8)
%     legend('Rising Field','Falling Field')
% %     s = plot(obj.raw(4).xUp,obj.raw(4).yUp,'b')
%     s.Color(4) = 0.45
%     xlabel('Magnetic field 1/H (T)')
%     ylabel('dM/dH (arb. units)')
%     title_string = sprintf('sample H, fourier window 20 to 60 T, %.2f K temp',obj.FFT.range.upTemp(ii).temp)
%     title(title_string)%, lowpass filter')
% 
% end
%%

% ylim([0 1.6e6])

%label for P up
% text(1410,6e5,'1419 T','FontSize',12,'Rotation',90);
% text(3470,10e5,'3470 T','FontSize',12,'Rotation',90);

%labels for H up
% text(195,10.4e5,'196 T','FontSize',12,'Rotation',90);
% text(398,8.6e5,'398 T','FontSize',12,'Rotation',90);
% text(618,7e5,'618 T','FontSize',12,'Rotation',90);
% text(822,5e5,'822 T','FontSize',12,'Rotation',90);

%labels for H down
% text(135,12.5e5,'135 T','FontSize',12,'Rotation',90);
% text(354,14.1e5,'354 T','FontSize',12,'Rotation',90);
% text(611,6e5,'611 T','FontSize',12,'Rotation',90);
% text(838,5e5,'838 T','FontSize',12,'Rotation',90);
% text(2171,6.3e5,'2171 T','FontSize',12,'Rotation',90);
% text(2500,6.5e5,'2500 T','FontSize',12,'Rotation',90);
% text(3180,5e5,'3180 T','FontSize',12,'Rotation',90);
% text(3870,4e5,'3870 T','FontSize',12,'Rotation',90);
% text(4382,4e5,'4388 T','FontSize',12,'Rotation',90);

% set(p1Lab,)
% ylim([0, 350])

