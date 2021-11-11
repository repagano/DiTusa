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
% startFields = 15; endFields = 1/(-iFFspan*2+1/startFields)%% Here the
% Fourier transform is performed

dHvA.FFTload(obj,endFields,iFFspan,[10 5e3],'extFFT');%7e4],'extFFT')%
% %% for ii = 1%1:length(obj.FFT.range.upTemp)
%     figure
%     plot(1./obj.FFT.range.upTemp(ii).xspl,obj.FFT.range.upTemp(ii).yspl)
%     xlabel('Frequency (T)') ylabel('Frequency after bs')
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
