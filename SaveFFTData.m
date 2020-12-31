close all
clc
%% Code used for saving fourier transform
f = obj.FFT.range.upTemp(2).f';
FFT = obj.FFT.range.upTemp(2).FFT';
data1 = [f,FFT];
save('FT35to55_temp0.96','data1','-ascii')

f = obj.FFT.range.upTemp(3).f';
FFT = obj.FFT.range.upTemp(3).FFT';
data1 = [f,FFT];
save('FT35to55_temp1.52','data1','-ascii')

f = obj.FFT.range.upTemp(7).f';
FFT = obj.FFT.range.upTemp(7).FFT';
data1 = [f,FFT];
save('FT35to55_temp5.71.txt','data1','-ascii')

%%
T = table(['M';'F';'M'],[45;41;36],...
    {'New York, NY';'San Diego, CA';'Boston, MA'},[true;false;false]) 

%%
L1 = load('FT35to55_temp0.96')%'FT25.67to35_temp0.96')%
L2 = load('FT35to55_temp1.52')%'FT25.67to35_temp1.52')%
L3 = load('FT35to55_temp5.71')%'FT25.67to35_temp5.71')%
% L1 = ans;

figure
plot(L1(:,1),L1(:,2),L2(:,1),L2(:,2),L3(:,1),L3(:,2))

%%
L1 = load('FT25.67to35_temp0.96')%
L2 = load('FT25.67to35_temp1.52')%
L3 = load('FT25.67to35_temp5.71')%
% L1 = ans;

figure
plot(L1(:,1),L1(:,2),L2(:,1),L2(:,2),L3(:,1),L3(:,2))


%% Code used for saving and loading AoT and Temp data over multiple ranges 
%% daving data
clear AoT1 AoT2 AoT3 AoT4 AoT5 cf2 cf
cf = [];
for ii = 1:length(obj.FFT.range)
   
     AoT1(:,ii+1) = obj.mass.range(ii).upPeak(1).AoT;
     AoT2(:,ii+1) = obj.mass.range(ii).upPeak(2).AoT;
     AoT3(:,ii+1) = obj.mass.range(ii).upPeak(3).AoT;
     AoT4(:,ii+1) = obj.mass.range(ii).upPeak(4).AoT;
     AoT5(:,ii+1) = obj.mass.range(ii).upPeak(5).AoT;
     cfV = strcat('CF',num2str(mean(obj.FFT.range(ii).upTemp(1).range)));
     cf = [cf,{cfV}];
     cf2(:,ii) = obj.FFT.range(ii).upTemp(1).range'

    
end
cf ={'1' '2' '3' '4' '5' '6'}
AoT1(:,1) = [obj.mass.range(ii).upPeak(5).temp];
AoT2(:,1) = [obj.mass.range(ii).upPeak(5).temp];
AoT3(:,1) = [obj.mass.range(ii).upPeak(5).temp];
AoT4(:,1) = [obj.mass.range(ii).upPeak(5).temp];
AoT5(:,1) = [obj.mass.range(ii).upPeak(5).temp];

%%
AoT1t = array2table(AoT1)
AoT1t.Properties.VariableNames(1) = {'temp'}% 'a1.1' 'a2' 'a3' 'a4' 'a5' 'a6'}
writetable(AoT1t,'AoT1.csv','Delimiter',',','QuoteStrings',true)

AoT2t = array2table(AoT2)
AoT2t.Properties.VariableNames(1) = [{'temp'}]%,cf]
writetable(AoT2t,'AoT2.csv','Delimiter',',','QuoteStrings',true)

AoT3t = array2table(AoT3)
AoT3t.Properties.VariableNames(1) = [{'temp'}]%,cf]
writetable(AoT3t,'AoT3.csv','Delimiter',',','QuoteStrings',true)

AoT4t = array2table(AoT4)
AoT4t.Properties.VariableNames(1) = [{'temp'}]%,cf]
writetable(AoT4t,'AoT4.csv','Delimiter',',','QuoteStrings',true)

AoT5t = array2table(AoT5)
AoT5t.Properties.VariableNames(1) = [{'temp'}]%,cf]
writetable(AoT5t,'AoT5.csv','Delimiter',',','QuoteStrings',true)

cft = array2table(cf2)
% AoT5.Properties.VariableNames(1:ii+1) = ['1' '2']
writetable(cft,'cf.csv','Delimiter',',','QuoteStrings',true)

%%
% csvwrite('filename.csv',AoT1)
save('AoT_at_610T.txt','AoT1','-ascii')
save('AoT_at_775T.txt','AoT2','-ascii')
save('AoT_at_906T.txt','AoT3','-ascii')
save('AoT_at_2765T.txt','AoT4','-ascii')
save('AoT_at_3700T.txt','AoT5','-ascii')
% txt1 = 'AoT peak at  %5.2 T';
% title = sprintf('AoT peak at  %d T',mean(obj.mass.range(ii).upPeak(1).maxFreq))

%% loading data
clear 
close all
clc

addpath('C:\Users\repag\Documents\MATLAB\AoTvals')
fileID = fopen('C:\Users\repag\Documents\MATLAB\AoTvals\AoT1.csv')
dataArray = textscan(fileID,'%f%f%f%f%f%f%[^\n\r]','Delimiter',',','TextType',...
     'string','HeaderLines' ,1,'ReturnOnError', false,'EndOfLine', '\r\n')
fclose(fileID)
temp = dataArray{1}
figure
for ii = 2:length(dataArray)
    AoT610(:,ii-1) = dataArray{ii}
    subplot(2,3,ii-1)
    plot(temp,AoT610(:,ii-1),'*')
    suptitle('AoT610')
end

fileID = fopen('C:\Users\repag\Documents\MATLAB\AoTvals\AoT2.csv')
dataArray = textscan(fileID,'%f%f%f%f%f%f%[^\n\r]','Delimiter',',','TextType',...
     'string','HeaderLines' ,1,'ReturnOnError', false,'EndOfLine', '\r\n')
fclose(fileID)
temp = dataArray{1}
figure
for ii = 2:length(dataArray)
    AoT795(:,ii-1) = dataArray{ii}
    subplot(2,3,ii-1)
    plot(temp,AoT795(:,ii-1),'*')
    suptitle('AoT795')
end

%% 
fileID = fopen('C:\Users\repag\Documents\MATLAB\AoTvals\AoT3.csv')
dataArray = textscan(fileID,'%f%f%f%f%f%f%[^\n\r]','Delimiter',',','TextType', 'string',...
    'HeaderLines' ,1,'ReturnOnError', false,'EndOfLine', '\r\n')
fclose(fileID)
temp = dataArray{1}
figure
for ii = 2:length(dataArray)
    AoT850(:,ii-1) = dataArray{ii}
    subplot(2,3,ii-1)
    plot(temp,AoT850(:,ii-1),'*')
    suptitle('AoT850')
end

fileID = fopen('C:\Users\repag\Documents\MATLAB\AoTvals\AoT4.csv')
dataArray = textscan(fileID,'%f%f%f%f%f%f%[^\n\r]','Delimiter',',','TextType', 'string',...
    'HeaderLines' ,1,'ReturnOnError', false,'EndOfLine', '\r\n')
fclose(fileID)
temp = dataArray{1}
figure
for ii = 2:length(dataArray)
    AoT2750(:,ii-1) = dataArray{ii}
    subplot(2,3,ii-1)
    plot(temp,AoT2750(:,ii-1),'*')
    suptitle('AoT 2750')
end

fileID = fopen('C:\Users\repag\Documents\MATLAB\AoTvals\AoT5.csv')
dataArray = textscan(fileID,'%f%f%f%f%f%f%[^\n\r]','Delimiter',',','TextType', 'string',...
    'HeaderLines' ,1,'ReturnOnError', false,'EndOfLine', '\r\n')
fclose(fileID)
temp = dataArray{1}
figure
for ii = 2:length(dataArray)
    AoT3700(:,ii-1) = dataArray{ii}
    subplot(2,3,ii-1)
    plot(temp,AoT3700(:,ii-1),'*')
    suptitle('AoT3700')
end



%% Code used for saving full AoT (15 to 55) data, then ploting python mass calc
%% Save data
AoTFull = obj.mass.range.upPeak(1).temp
for ii = 1:length(obj.mass.range.upPeak)
    AoTFull(:,ii+1) = obj.mass.range.upPeak(ii).AoT
end
 
AoTFullt = array2table(AoTFull)
AoTFullt.Properties.VariableNames(1) = {'temp'}%, {'a1'}, {'a2'}, {'a3'}, {'a4'}, {'a5'}, {'a6'}
writetable(AoTFullt,'AoTFull.csv','Delimiter',',','QuoteStrings',true)

%% plot data

maxFreq = [ 606.44  ,  771.0279,  888.5672, 2763.2   , 3724.6   ];
massFull = [0.32288387, 0.51571273, 0.37455657, 0.65685059, 1.49424404];
massFull_error = [0.04712786, 0.07486916, 0.32655534, 0.10476229, 0.17324713];

figure
errorbar(maxFreq,massFull,massFull_error,'*')
xlabel('Frequency (T)')
ylabel('Effective mass')

fullAoTPfit = [1.36781632e+02, 8.80190220e+01, 5.50746089e+01, 3.67020108e+01,...
        2.42157795e+01, 1.93564116e+01, 1.20568834e+01, 7.12599692e+00,...
        4.65779221e+00, 3.08354750e+00, 2.05276283e+00, 9.15230904e-01;
       7.28715730e+01, 4.66315616e+01, 2.87632101e+01, 1.86502111e+01,...
        1.16432899e+01, 8.88434235e+00, 4.78725590e+00, 2.24795063e+00,...
        1.17701106e+00, 6.18171822e-01, 3.24941867e-01, 8.98317544e-02;
       4.50878082e+01, 2.89778299e+01, 1.80737351e+01, 1.19703146e+01,...
        7.79893041e+00, 6.16728163e+00, 3.71188517e+00, 2.07315542e+00,...
        1.28079256e+00, 7.98427507e-01, 4.99467661e-01, 1.96095708e-01;
       6.44518110e+01, 4.10097105e+01, 2.49324592e+01, 1.57322600e+01,...
        9.31112991e+00, 6.80315122e+00, 3.22365493e+00, 1.25250887e+00,...
        5.51717485e-01, 2.43251886e-01, 1.07268967e-01, 2.08613939e-02;
       1.89694913e+02, 1.13674000e+02, 5.97117969e+01, 2.91898023e+01,...
        1.10383550e+01, 5.79610381e+00, 1.13533646e+00, 1.34577685e-01,...
        2.08985398e-02, 3.24533079e-03, 5.03966887e-04, 1.21531279e-05];
    
figure
plot(AoTFull(:,1),AoTFull(:,6))
hold on
plot(AoTFull(:,1),fullAoTPfit(5,:),'*')
xlabel('Temp (K)')
ylabel('Amplitude/T (arb. units)')


%%
clear
clc

load('AoT_at_610T.txt')
load('AoT_at_775T.txt')
load('AoT_at_906T.txt')
load('AoT_at_2765T.txt')
load('AoT_at_3700T.txt')

AoT1t = array2table(AoT_at_610T)
AoT1t.Properties.VariableNames(1) = {'temp'}% 'a1.1' 'a2' 'a3' 'a4' 'a5' 'a6'}
writetable(AoT1t,'AoT1.csv','Delimiter',',','QuoteStrings',true)

AoT1t = array2table(AoT_at_775T)
AoT1t.Properties.VariableNames(1) = {'temp'}% 'a1.1' 'a2' 'a3' 'a4' 'a5' 'a6'}
writetable(AoT1t,'AoT2.csv','Delimiter',',','QuoteStrings',true)

AoT1t = array2table(AoT_at_906T)
AoT1t.Properties.VariableNames(1) = {'temp'}% 'a1.1' 'a2' 'a3' 'a4' 'a5' 'a6'}
writetable(AoT1t,'AoT3.csv','Delimiter',',','QuoteStrings',true)

AoT1t = array2table(AoT_at_2765T)
AoT1t.Properties.VariableNames(1) = {'temp'}% 'a1.1' 'a2' 'a3' 'a4' 'a5' 'a6'}
writetable(AoT1t,'AoT4.csv','Delimiter',',','QuoteStrings',true)

AoT1t = array2table(AoT_at_3700T)
AoT1t.Properties.VariableNames(1) = {'temp'}% 'a1.1' 'a2' 'a3' 'a4' 'a5' 'a6'}
writetable(AoT1t,'AoT5.csv','Delimiter',',','QuoteStrings',true)
