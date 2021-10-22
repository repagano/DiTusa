 close all
clear
clc

%% Here the files to be loaded are chosen.
addpath C:\Users\repag\Documents\MATLAB\DiTusa
loc = 'C:\Users\repag\Documents\DiTusa\yu NbGe2\data2\';
temps = [3.49 2.67 2.18 2.38 2.03 1.77 1.54 1.55 1.11 3.98 3.06];%Temperatures
% corresponding to files to be loaded
filesH = dir(strcat(loc,'*H001.dat'));
varTemps = temps;
I = ismember(temps,varTemps);
temps = temps(I);%
files = filesH(I);

%%
obj = dHvA(loc,temps,files);

%% Here the 1/H windows that are selected and fourier transformed are defined
clc
close all
endFields = 60;%the endFields are the maximum field values of each window

iFFspan = abs((1/20-1/60));%Here the width of the 1/H window is defined 
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
for ii = 2%1:length(obj.FFT.range.upTemp)
    figure
%     plot(1./obj1.FFT.range.upTemp(ii).xspl,obj1.FFT.range.upTemp(ii).yspl,'r','LineWidth',1.1)
%     hold on
    s = plot(1./obj.FFT.range.upTemp(ii).xspl,obj.FFT.range.upTemp(ii).yspl,'b')
%     s = plot(obj.raw(4).xUp,obj.raw(4).yUp,'b')
    s.Color(4) = 0.45
    xlabel('Magnetic field 1/H (T)')
    ylabel('dM/dH (arb. units)')
%     legend('Raw','Low-pass filter')
    title_string = sprintf('fourier window 20 to 60 T, rising field, %.2f K temp',obj.FFT.range.upTemp(ii).temp)
    title(title_string)%, lowpass filter')

end
%%
tv = obj.FFT.range.upTemp(1).f';
leg = []
figure
for ii = 1:length(obj.FFT.range(end).upTemp)%4%
%     figure
    c = parula(length(obj.FFT.range.upTemp));
    plot(obj.FFT.range.upTemp(ii).f,obj.FFT.range.upTemp(ii).FFT,'LineWidth',1)%,'Color',c(ii,:)
    hold on
    txt = {strcat(num2str(obj.FFT.range.upTemp(ii).temp),' K')};
    leg = [leg,txt];
    hold on
%     pause(1)
    
%     tv(:,ii+1) = obj.FFT.range.upTemp(ii).FFT';
end
xlabel('Frequency (T)')
ylabel('Fourier amp. (arb units)')
legend(leg)
title('range of 20 to 60 T, rising field')
% ylim([0, 350])

%%
endV = length(obj.raw)
figure
for ii = 1:endV
    plot(obj.raw(ii).xDown,obj.raw(ii).yDown)
    hold on
end



%%
fH = strcat(loc,filesH(7).name);
dH = load(fH);
fP = strcat(loc,filesH(7).name)
dP = load(f2);

figure
plot(dH(:,1),dH(:,2))
title('H')
% hold on 
figure
plot(dP(:,1),dP(:,2))%./d2(:,2))
title('P')