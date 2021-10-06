close all
clear
clc
%% Here the files to be loaded are chosen.
addpath C:\Users\repag\Documents\MATLAB\DiTusa
loc = 'C:\Users\repag\Documents\DiTusa\yu NbGe2\data\';
temps = [3.49 32.67 2.18 2.38 2.03 1.77 1.54 1.55 1.11 3.98 3.06];%Temperatures
% corresponding to files to be loaded
filesH = dir(strcat(loc,'*H001.dat'));
filesP = dir(strcat(loc,'*P001.dat'));
% filesFFT = dir(strcat(loc,'FT*'));
% here set the variable varTemps equal to the temperatures you wish to
% load. If you set varTemps = temps all files are loaded.
% [~,I] = max(temps)
% varTemps = temps;
% I = ismember(temps,varTemps);
% temps = temps(I);
% files = files(I);
% filesFFT = filesFFT(I);
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
%%
obj = dHvA(loc,temps,filesH,filesFFT);

%% Here the 1/H windows that are selected and fourier transformed are defined
endFields = 35;%the endFields are the maximum field values of each window

iFFspan = abs((1/20-1/60));%Here the width of the 1/H window is defined 
% startFields = 15;
% endFields = 1/(-iFFspan*2+1/startFields)%% Here the Fourier transform is performed  

dHvA.FFTload(obj,endFields,iFFspan,[200 5e3],'extFFT');%7e4],'extFFT')%