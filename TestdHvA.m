clear
close all
clc
 %% Here the files to be loaded are chosen.
addpath C:\Users\repag\Documents\MATLAB\DiTusa
loc = 'C:\Users\repag\Documents\DiTusa\MagLab\031919_PdGa\031919\';
temps = [1.52 .62 5.71 3.96 8 10 12 14 18 3.27 2.24 .96];%Temperatures
% corresponding to files to be loaded
files = dir(strcat(loc,'*.ASC'));
filesFFT = dir(strcat(loc,'FT*'));
% here set the variable varTemps equal to the temperatures you wish to
% load. If you set varTemps = temps all files are loaded.
% [~,I] = max(temps)
varTemps = temps;
I = ismember(temps,varTemps);
temps = temps(I);
files = files(I);
filesFFT = filesFFT(I);

%% Here the raw data is loaded into dHvA class
obj = dHvA(loc,temps,files,filesFFT);

%% Here the 1/H windows that are selected and fourier transformed are defined
endFields = 35;%the endFields are the maximum field values of each window

iFFspan = abs((1/15-1/35));%Here the width of the 1/H window is defined 
% startFields = 15;
% endFields = 1/(-iFFspan*2+1/startFields)%% Here the Fourier transform is performed  

dHvA.FFTload(obj,endFields,iFFspan,[200 5e3],'extFFT');%7e4],'extFFT')%


%%
% % T = csvread('C:\Users\repag\Documents\MATLAB\Manuscript 21 data\dMdH.csv')%array2table([1./obj.FFT.range.upTemp(4).xspl',obj.FFT.range.upTemp(4).yspl']);
% tv = [table2array(T),obj.FFT.range.upTemp(4).yspl'];
% figure
% plot(tv(:,1),tv(:,2),tv(:,1),tv(:,3))
% T = array2table(tv)
% T.Properties.VariableNames = [{'H'},{'dM/dH no filter'},{'dM/dH 10000 filter'}]
% writetable(T,'C:\Users\repag\Documents\MATLAB\Manuscript 21 data\dMdH.csv','Delimiter',',','QuoteStrings',true)
%% Here the effective masses of the different peaks are calculated 


% peakRange is an matrix whose rows define the range around each peak 
% identified on the Fourier spectrum. The columns correspond to different
% peaks. The range should be large enough to allow for variation in the position 
% of the peak as the different ranges and temperatures are scanned over, but 
% not be so large that another peak is selected. 

peakRange =  [2666, 2858; 3232, 3316; 3317.5, 3392.5; 3615, 3823];% 15to55% [2666, 2858; 3194, 3400; 3558, 3875];% 15to35%
% peakRange
% dHvA.massLoad runs the method that calculates the effective masses.
dHvA.massLoad(obj,peakRange);

%%

for ii = 1:length(obj.mass.range.upPeak)
    figure
    plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).AoT)
end

%%
for ii = 1:length(obj.raw)
    figure 
    plot(obj.raw(ii).xUp,obj.raw(ii).yUp,'r')
    hold on
    plot(obj.raw(ii).xDown,obj.raw(ii).yDown,'b')
end

%%
tv = obj.FFT.range.upTemp(1).f';
leg = []
figure
for ii = 1:length(obj.FFT.range(end).upTemp)%4%
%     figure
    c = parula(length(obj.FFT.range.upTemp));
    plot(obj.FFT.range.upTemp(ii).f,obj.FFT.range.upTemp(ii).FFT,'LineWidth',1,'Color',c(ii,:))
    hold on
    txt = {strcat(num2str(obj.FFT.range.upTemp(ii).temp),' K')};
    leg = [leg,txt];
    hold on
    pause(1)
    
    tv(:,ii+1) = obj.FFT.range.upTemp(ii).FFT';
end
xlabel('Frequency (T)')
ylabel('Fourier amp. (arb units)')
legend(leg)
title('range of 15 to 35 T, rising field')
%%
T = array2table(tv)
T.Properties.VariableNames(1) = {'Frequency (T)'};
T.Properties.VariableNames(2:end) = leg
writetable(T,'C:\Users\repag\Documents\MATLAB\Manuscript 21 data\FFT_15to35.csv','Delimiter',',','QuoteStrings',true)
%%
for ii = 1:length(obj.FFT.range.upTemp)
    figure
%     plot(1./obj1.FFT.range.upTemp(ii).xspl,obj1.FFT.range.upTemp(ii).yspl,'r','LineWidth',1.1)
%     hold on
    s = plot(obj.FFT.range.upTemp(ii).xspl,obj.FFT.range.upTemp(ii).yspl,'b')
%     s = plot(obj.raw(4).xUp,obj.raw(4).yUp,'b')
    s.Color(4) = 0.45
    xlabel('Magnetic field 1/H (T)')
    ylabel('dM/dH (arb. units)')
%     legend('Raw','Low-pass filter')
    title('15 to 55 T, 2.24 K')%, lowpass filter')
    pause(1)
end
%% generate data
% 
% x = 0:.1:10;
% y = x.*x + randn(size(x));
% w = linspace(.5, .7,length(x));
% x = x(:);
% y = y(:);
% w = w(:);
% %plot data
% plot(x,y,'.');
% %fit
% ft = fittype('poly2');
% cf = fit(x,y,ft,'Weight',w);
% % Plot fit
% hold on
% plot(cf,'fit',0.95);
% 
% %%
% for ii = 1:length(obj.mass.range.upPeak)
% %     figure
% %     plot(obj.mass.range.upPeak(ii).temp)
%     
%     figure
%     plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).A)
%     title(strcat('max freq',num2str(mean(obj.mass.range.upPeak(ii).maxFreq))))
%     xlabel('temp')
%     ylabel('A')
%     
%     
%     figure
%     plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).AoT,'*')
%     hold on
%     plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).AoTfitV)
%     title(strcat('max freq',num2str(mean(obj.mass.range.upPeak(ii).maxFreq))))
%     xlabel('temp')
%     ylabel('AoT')
%     
%     figure
%     mf = mean(obj.mass.range.upPeak(ii).maxFreq,obj.mass.range.upPeak(ii))
%     plot(obj.mass.range.upPeak(ii).m,'*' )
% end
%%
figure
for ii = [1,2,4,5]
    plot(mean(obj.mass.range.upPeak(ii).maxFreq),obj.mass.range.upPeak(ii).phase(3),'*')
    hold on
end
xlabel('Frequency (T)')
ylabel('phase')
title('phase calc, range 0.049 1/T')



%% phase plot with fixed frequency
figure 
leg = [];
for ii = 1:length(obj.mass.range.upPeak)
    chioi = ii
    mF(ii) = obj.mass.range.upPeak(ii).maxFreq(3)
    
    for ii = 1:length(obj.FFT.range.upTemp)
        chjj = ii
        dif = abs(obj.FFT.range.upTemp(ii).f - mF(ii));
%         figure
%         plot(dif)
        [~,phI] = min(dif);
        chphI = obj.FFT.range.upTemp(ii).f(phI)
        tempV(ii) = obj.FFT.range.upTemp(ii).temp;
        phaseV(ii,ii) = obj.FFT.range.upTemp(ii).phase(phI)
        
%         dif = diff(phaseV(ii,:))
        
        for kk = 1:length(phaseV(1,:))-1
            kk
            dif = diff([phaseV(ii,kk) phaseV(ii,kk+1)])
            sin = sign(phaseV(ii,kk)-phaseV(ii,kk+1))
            cnt = 0;
            if abs(dif) >= 2*pi*.65
                cnt = cnt+1
                phaseV(ii,kk+1) = phaseV(ii,kk+1)+sin*pi*2;
            end 
            
        end

    end
    
    
    plot(tempV,phaseV(ii,:))
    hold on
    txt = {strcat('Peak ', num2str(ii))};
    leg = [leg,txt];
                       
end
% sf = num2str(
legend(leg)
xlabel('Temp')
ylabel('phase')
title('phase at fixed frequency, field range 15 to 55 T')

%% Fourier transform from LANL and MATLAB. Fig 2 in document created using 
% this section with endfield  == 40 and iFFspan=abs((1/10-1/40))=0.075.
tempFFT = 1:11;
pos = 1;
% figure
leg = []
for ii = 1:12
    figure
%     subplot(2,2,pos)
%     hold on
%     plot(obj.raw(ii).f,obj.raw(ii).FFT,'LineWidth',1.5)
%     hold on
    plot(obj.FFT.range.upTemp(ii).f,obj.FFT.range.upTemp(ii).FFT)
    pos = pos+1;
    temp(ii) = obj.FFT.range.upTemp(ii).temp;
%     for jj = 1:4
%         plot(obj.mass.raw(jj).maxFreq(ii),obj.mass.raw(jj).A(ii),'*')
%     end
    txt = {strcat('T ', num2str(obj.FFT.range.upTemp(ii).temp))};
    leg = [leg,txt];
%     title(strcat(num2str(temp(ii)), ' K, fourier window 15 to 50 T'))
    % legend('LANL fourier transform','MATLAB fourier transform')
    xlabel('dHvA frequency, F (T)')
    ylabel('Amplitude (arb. units)')
    pause(1)
    %used to make Fig.2.1
end
% title(strcat(num2str(temp 'K, fourier window 15 to 35 T')
% legend('LANL fourier transform','MATLAB fourier transform')
set(gca,'Ytick',0:2:10)
set(gca,'Xtick',0:1000:4500)
xlabel('dHvA frequency, F (T)')
ylabel('Amplitude (arb. units)')
legend(leg)
%used to make Fig.2.1

%% Plot masses vs frequency
figure

for ii = 1:length(obj.FFT.range)
%     jj
    for ii = 1:length(obj.mass.range(1).upPeak)
        massPlot(ii) = obj.mass.range(ii).upPeak(ii).m;
        freqPlot(ii) = mean(obj.mass.range(ii).upPeak(ii).maxFreq);
        errorPlot(ii) = mean(obj.mass.range(ii).upPeak(ii).AoTrms);
        rawMassPlot(ii) = obj.mass.raw(ii).m;
        rawFreqPlot(ii) = mean(obj.mass.raw(ii).maxFreq);
        rawErrorPlot(ii) = mean(obj.mass.raw(ii).AoTrms);
    end
    leg{ii} = num2str(obj.FFT.range(ii).upTemp(1).range);
    errorbar(freqPlot,massPlot,errorPlot,'o')
    hold on
    errorbar(rawFreqPlot,rawMassPlot,rawErrorPlot,'*')%freqPlot,massPlot,'*',
    hold on

end
% errorbar(rawFreqPlot,rawMassPlot,rawErrorPlot,'*')
title('Fourier Window 15 to 35')
ylabel('Effective mass, m*')
xlabel('F(T)')
% legend('m*','m* from LANL ft')
%used to make fig.4

%% Plot AoT vs. temp from my MATLAB fft

for ii = 1:length(obj.FFT.range)
%     jj
figure
    for ii = 1:length(obj.mass.range(1).upPeak)
        aotfitV = obj.mass.range.upPeak(ii).AoTfitV;
        ait = obj.mass.range(1).upPeak(ii).AoT; 
        figure
%         subplot(2,2,ii)
%         figure(1)
        plot(obj.mass.range(ii).upPeak(ii).temp,aotfitV,obj.mass.range(ii).upPeak(ii).temp,ait,'*')
        hold on
%         (obj.mass.range(jj).upPeak(ii).temp,aotfitV1,obj.mass.range(jj).upPeak(ii).temp,ait1,'*')
%         plot(obj.mass.range(jj).upPeak(ii).temp,(obj.mass.range(jj).upPeak(ii).AoT),'*',...
%             obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).AoTfitV)
%         maxFreq = obj.mass.range.upPeak(kk).maxFreq;
%         title(num2str(mean(maxFreq)))
        maxFreq = obj.mass.range(ii).upPeak(ii).maxFreq;
        title(strcat('Peak Frequency (T) =',num2str(round(mean(maxFreq)))));
        xlabel('Temperature (K)')
        ylabel('Amplitude/T (arb. units)')
        
    end
    suptitle('MATLAB ft, Fourier Window 15 to 35 T')
end
%used to make fig. 3.a
legend('1','1*')




%% plot dM/dH vs 1/H. Fig 1 in document created using this section with
% endfield  == 40 and iFFspan=abs((1/10-1/40))=0.075.
temp = 1:12;%[1];
for ii = length(obj.FFT.range)
    for ii = temp
        
        figure(ii)
        subplot(2,1,1)
%         plot(obj.FFT.range(jj).upTemp(ii).x,obj.FFT.range(jj).upTemp(ii).y,'Linewidth',1.2)
%         hold on
        plot(obj.FFT.range(ii).upTemp(ii).xspl,obj.FFT.range(ii).upTemp(ii).yspl,'Linewidth',.9)
        xlabel('1/H')
        ylabel('dM/dH')
        legend('dM/dH vs 1/H','Spline fit')
        subplot(2,1,2)
        plot(1./obj.FFT.range(ii).upTemp(ii).x,obj.FFT.range(ii).upTemp(ii).y,'Linewidth',1)
        hold on
%         subplot(2,1,2)
        plot(1./obj.FFT.range(ii).upTemp(ii).xspl,obj.FFT.range(ii).upTemp(ii).yspl,'Linewidth',.4)
        xlabel('H')
        ylabel('dM/dH')
        legend('dM/dH vs 1/H','Spline fit')
    end      
end
suptitl1 = strcat(num2str(obj.raw(ii).temp),'K');
suptitle(strcat('Fig 1: T=',suptitl1))
%used to make Fig.1
%%

for ii = 1:length(obj.FFT.range.upTemp)
    ii
    fdif = length(obj.FFT.range.upTemp(ii).f)
%     ch2 = length(obj.FFT.range.upTemp(ii+1).f)
    
%     figure
%     plot(fdif)    
end

%%
figure
for ii = 1:5
    figure
%     subplot(2,1,1)
%     plot(obj.mass.range(1).upPeak(ii).temp,obj.mass.range(1).upPeak(ii).phase)
%     subplot(2,1,2)
    plot(obj.mass.range(1).upPeak(ii).temp,obj.mass.range(1).upPeak(ii).maxFreq)
end


%%
figure
for ii = 1:5%length(obj.mass.range.upTemp
    
    subplot(3,2,ii)
    plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).maxFreq)
    xlabel('temp (K)')
    ylabel('frequency')
    title(strcat('peak ',num2str(ii)))
end
suptitle('max frequency with fixed frequency range')



%%
leg = [];
figure
for ii = 1:5
    plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).phase)
    hold on
    leg = [leg,{strcat('Peak ',num2str(ii))}]
end
legend(leg)
title('Peaks at max frequency')


    
%%
for ii = 1:12
%     subplot(1,2,1)
    figure
    plot(obj.raw(ii).f,obj.raw(ii).FFT,'LineWidth',1.5)
    hold on
    plot(obj.FFT.range.upTemp(ii).f,obj.FFT.range.upTemp(ii).FFT)
%     pos = pos+1;
%     for jj = 1:4
%         plot(obj.mass.raw(jj).maxFreq(ii),obj.mass.raw(jj).A(ii),'*')
%     end
end
title('Fig 2.1U: 0.96 K, fourier window 15 to 35 T')
legend('LANL fourier transform','MATLAB fourier transform')
xlabel('dHvA frequency, F (T)')
ylabel('Amplitude (arb. units)')

% subplot(1,2,2)
for ii =1:4
    plot(obj.mass.raw(ii).temp,obj.mass.raw(ii).A,'LineWidth',2)
    hold on
    plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).A)
    
end
title('Fig 2.2U')
xlabel('Temp (K)')
ylabel('Fourier Amplitude (arb. units)')
legend('Frequency ~600 matlab','Frequency ~600 Lanl','Frequency ~775 matlab','Frequency ~775 Lanl',...
    'Frequency ~2770 matlab','Frequency ~2770 Lanl','Frequency ~3700 matlab','Frequency ~3700 Lanl')
%% Compleat phase plot
% this section with endfield  == 40 and iFFspan=abs((1/10-1/40))=0.075.
% tempFFT = 1:11;
pos = 1;
    
for ii = 1:length(obj.FFT.range(1).upTemp)
%     obj.FFT.range.upTemp(ii).phase(1:2:end) = obj.FFT.range.upTemp(ii).phase(1:2:end)+pi
%     subplot(2,2,pos)
    figure
    plot(obj.FFT.range.upTemp(ii).f,obj.FFT.range.upTemp(ii).phase)
    pos = pos+1;
   pause(1)

%     for jj = 1:4
%         plot(obj.mass.raw(jj).maxFreq(ii),obj.mass.raw(jj).A(ii),'*')
%     end
end
% title('Fig 2: 0.96 K, fourier window 15 to 35 T')
% legend('LANL fourier transform','MATLAB fourier transform')
% xlabel('dHvA frequency, F (T)')
% ylabel('Amplitude (arb. units)')



%% AoT vs temp plot from LANL fft
figure
for kk = 1:5
    subplot(3,2,kk)
    plot(obj.mass.raw(kk).temp,obj.mass.raw(kk).AoT,'*',...
        obj.mass.raw(kk).temp,obj.mass.raw(kk).AoTfitV)
    maxFreq = obj.mass.raw(kk).maxFreq;
    title(strcat('Peak Frequency (T) =',num2str(round(mean(maxFreq)))));
    xlabel('Temperature (K)')
    ylabel('Amplitude/T (arb. units)')
    
end

suptitle('LANL ft, Fouriet Window 15 to 35 T')
%used to make fig. 3.b

%% plot Fourier transform for all temps 
cnt = [1 2];% 3 4 5 6 7 10];
sp = 1;
figure
for ii = 1:length(obj.FFT.range)
%     subplot(2,3,sp)
    for ii = 1:length(obj.raw)
%         figure(jj)
        leg{ii} = [strcat(num2str(obj.FFT.range(ii).upTemp(ii).temp), 'K')];
%         subplot(2,4,sp)
        plot(obj.FFT.range(ii).upTemp(ii).f,obj.FFT.range(ii).upTemp(ii).FFT)
        hold on
        title(num2str(obj.FFT.range(ii).upTemp(ii).range))
    
    end

sp = sp+1;
titl1 = strcat('CF = ',num2str(round(mean(obj.FFT.range(ii).upTemp(1).range))));
titl = strcat(titl1,'T');
title(titl)
xlabel('frequency (T)')
ylabel('Amplitude (arb. units)')

end
legend(leg)

%% plot fourier transform for temps individually 
cnt = [1 2];% 2 3 4 5 6 7 10];
sp = 1;
figure
for ii = 1:length(obj.FFT.range)%cnt
    ii
    
    subplot(2,3,sp)    
    plot(obj.FFT.range(ii).upTemp(2).f,obj.FFT.range(ii).upTemp(2).FFT,'r')
    hold on
    plot(obj.FFT.range(ii).upTemp(3).f,obj.FFT.range(ii).upTemp(3).FFT,'b')
    plot(obj.FFT.range(ii).upTemp(7).f,obj.FFT.range(ii).upTemp(7).FFT,'Color',[0 .7 .3])
    titl1 = strcat('CF = ',num2str(round(mean(obj.FFT.range(ii).upTemp(1).range))))
    titl = strcat(titl1,'T')
    title(titl)
    xlabel('frequency (T)')
    ylabel('Amplitude (arb. units)')
    sp = sp+1
end
legend('0.96 K','1.52 K','5.71 K','Location','northeast')



%% Plot mass vs center field and frequency vs center field
% figure
for ii = 1:length(obj.FFT.range)
    CF(ii) = mean(obj.FFT.range(ii).upTemp(1).range);
    for ii = [1:4]%:length(obj.mass.range(ii).upPeak)
        chjj = ii
        massCF(ii,ii) = obj.mass.range(ii).upPeak(ii).m;
        errorCF(ii,ii) = obj.mass.range(ii).upPeak(ii).AoTrms;
        maxPeakSTD(ii,ii) = std(obj.mass.range(ii).upPeak(ii).maxFreq);
        maxPeakCF(ii,ii) = mean(obj.mass.range(ii).upPeak(ii).maxFreq);
    end
end

figure(5)
figure(6)
for kk = [1 4] %:length(maxPeakCF(:,1))
    figure(5)
    
    errorbar(CF,maxPeakCF(kk,:),maxPeakSTD(kk,:),'*')
    hold on
    
    figure(6)
    title('center field vs mass')
%     plot(CF,massCF(kk,:),'*')
    errorbar(CF,massCF(kk,:),errorCF(kk,:),'-*')
%     legend(
    hold on
end
figure(5)
title('center field vs max peaks')
legend('Peek 1','Peek 2','Peek 4','Peek 5')
xlabel('Center Field (T)')
ylabel('Peek Frequency (T)')

figure(6)
title('center field vs mass')
legend('Peek 1','Peek 2','Peek 4','Peek 5')
xlabel('Center Field (T)')
ylabel('Effective Mass')


%% Plot max peak amplitude vs center peaks
peaks = [1 2 4 5];
figure(1)
figure(2)
for kk = 1:4
    for ii = 1:length(obj.mass.range)
        CF(ii) = mean(obj.FFT.range(ii).upTemp(1).range);
        A(:,ii) = obj.mass.range(ii).upPeak(kk).A;
        phase(:,ii) = obj.mass.range(ii).upPeak(kk).phase;
%         if isempty(Aval) == 0
%             A(:,ii) = Aval;
%         end
    end
   
    
    for ii = 1:length(A(:,1))
%         figure(1)
%         subplot(2,2,kk)
%         leg{jj} = strcat(num2str(obj.FFT.range(kk).upTemp(jj).temp), 'K');
%         plot(CF,A(jj,:))
%         xlabel('Center Field (T)')
%         ylabel('Fourier amp. (arb. units)')
%         hold on
        
        figure(2)
%         subplot(2,2,kk)
        leg{ii} = strcat(num2str(obj.FFT.range(kk).upTemp(ii).temp), 'K');
        plot(CF,phase(ii,:))
        xlabel('Center Field (T)')
        ylabel('phase (rad)')
        hold on
    end

    title(strcat('Peak',num2str(peaks(kk))))
%     leg
    for kk = 4
        figure(1)
        subplot(2,2,kk)
        legend(leg)
        
        figure(2)
        subplot(2,2,kk)
        legend(leg)
    end
%     pause(1)
end
%% calcluate and plot masses multipl times to test fit
% figure
% for jj = 1:15
%     chjj = jj
%     dHvA.massLoad(obj,peakRange)
%     for ii = 1:length(obj.mass.range(1).upPeak)
%         plot(mean(obj.mass.range(1).upPeak(ii).maxFreq),obj.mass.range(1).upPeak(ii).m,'*')%,'o',obj.mass.raw(ii).maxFreq,obj.mass.raw(ii).m,'*')%,obj.mass.down(ii).maxFreq,obj.mass.down(ii).m,'x',
%         hold on
%     end
% pause(.5)
% end
% xlabel('Frequency (T)')
% ylabel('mass')
% legend('my FT','LANL FT')


%% Phase and amplitude subplot
% for ii = 1:length
Fmax = [obj.mass.range.upPeak(1).maxFreq(1),obj.mass.range.upPeak(2).maxFreq(1),....
    obj.mass.range.upPeak(3).maxFreq(1),obj.mass.range.upPeak(4).maxFreq(1)];

Pmax = [obj.mass.range.upPeak(1).phase(1),obj.mass.range.upPeak(2).phase(1),....
    obj.mass.range.upPeak(3).phase(1),obj.mass.range.upPeak(4).phase(1)];

figure
subplot(2,1,1)
plot(obj.FFT.range.upTemp(1).f,obj.FFT.range.upTemp(1).phase,Fmax,Pmax,'*r')
title('phase from FFT')
subplot(2,1,2)
plot(obj.FFT.range.upTemp(1).f,obj.FFT.range.upTemp(1).FFT)
title('amplitude from FFT')

% figure
% plot(obj.FFT.range.upTemp(1).xspl,obj.FFT.range.upTemp(1).yspl)
%% Plot phase for 4 different peaks, all temps 
% 
% figure
% for ii = 1:4
%     ii
%     p=obj.mass.range.upPeak(ii).phase
%     I = p<pi
%     p(I) = p(I)+pi
%     phase(:,ii) = p
%     plot(obj.mass.range.upPeak(ii).temp,phase(:,ii))
%     hold on
% end

%% Plot phase for 4 different peaks, all temps 

figure
for ii = 1:4
    obj.mass.range.upPeak(ii).phase([2 9 10 11 12])...
        = obj.mass.range.upPeak(ii).phase([2 9 10 11 12])+pi
    plot(obj.mass.range.upPeak(ii).temp,obj.mass.range.upPeak(ii).phase)
    hold on
%     pause(1)
end
legend('peak 1','peak 1','peak 3','peak 4')
xlabel('Temp')
ylabel('phase') 
