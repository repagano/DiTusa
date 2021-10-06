clear
close all
clc
%% Here the files to be loaded are chosen.
loc = 'C:\Users\repag\Documents\DiTusa\MagLab\031919_PdGa\031919\';
temps = [1.52 .62 5.71 3.96 8 10 12 14 18 3.27 2.24 .96];%Temperatures
% corresponding to files to be loaded
files = dir(strcat(loc,'*.ASC'));
filesFFT = dir(strcat(loc,'FT*'));
% here set the variable temps equal to the data temperatures corresponding 
% to the files you wish to load in the expression find(temps == temps) by 
% setting temps == temps, all files are loaded.
[~,I] = find(temps == temps);
temps = temps(I);
files = files(I);
filesFFT = filesFFT(I);

%% Here the raw data is loaded into dHvA class
obj = dHvA(loc,temps,files,filesFFT);

%% Here the 1/H windows that are selected and fourier transformed are defined
% endFields = [15:8:55];%the endFields are the maximum field values of each window
% clear
% num = 5
iFFspan = abs((1/35-1/55))%Here the width of the 1/H window is defined 
endFields = [15:4:55]%flip([55:-6:15])%
% for jj = 1:length(endFields)
%     iendField = 1/endFields(jj)
%     iFFrange(jj,:) = [iendField+iFFspan iendField]
%     FFrange = 1./iFFrange
% %     pause(1)
% end
% cf = mean(FFrange,2)
% startFields = 15
% % 
% for ii = 1:num
%     endFields(ii) = 1/(-iFFspan*ii+1/startFields);
% %         ,1/(-iFFspan*2+1/startFields),...
% %     1/(-iFFspan*3+1/startFields),1/(-iFFspan*4+1/startFields)]
% end

%% Here the Fourier transform is performed  
dHvA.FFTload(obj,endFields,iFFspan,[200 4500],'lowpass',10000);

%% Plot FFT
figure
leg = [];
range = []
c = parula(length(obj.FFT.range));%
c = [0 0 0;c(1:end-1,:)]
c = pink(length(obj.FFT.range)+3)
for ii = 1:length(obj.FFT.range)
    range = [range; obj.FFT.range(ii).upTemp(1).range];
    cf(ii) = mean(obj.FFT.range(ii).upTemp(1).range);
    leg = [leg,{strcat(num2str(1/cf(ii)),' 1/T')}];
    plot(obj.FFT.range(ii).upTemp(2).f,obj.FFT.range(ii).upTemp(2).FFT,...
        'LineWidth',1.2,'Color',c(ii,:))
    hold on
end
leg = legend(leg)

title(leg,'1/B_m')
xlabel('Frequency (T)')
ylabel('Fourier amp. (arb. units)')
title('range 0.002 1/T range, interval 15:8:55')


%% Here the effective masses of the different peaks are calculated 
% peakRange is an matrix whose rows define the range around each peak 
% identified on the Fourier spectrum. The columns correspond to different
% peaks. The range should be large enough to allow for variation in the position 
% of the peak as the different ranges and temperatures are scanned over, but 
% not be so large that another peak is selected. 
peakRange = [490 708; 708 850; 2600 2900; 3546 3880]; %850 980;
% dHvA.massLoad runs the method that calculates the effective masses.
dHvA.massLoad(obj,peakRange);

%%
clc  
saveAoT1(obj,'C:\Users\repag\Documents\MATLAB\AoTvals\AoTrng15_1_55_25to55\','AoT15_1_55rng')
%%
close all
for ii = 1:length(obj.FFT.range)
    leg = [];
    figure
   for jj = 1:length(obj.FFT.range(1).upTemp)  
       
        plot(obj.FFT.range(ii).upTemp(jj).f,obj.FFT.range(ii).upTemp(jj).FFT)
        hold on
        leg = [leg,{strcat(num2str(obj.FFT.range(ii).upTemp(jj).temp),' K')}];
        xlabel('Frequency (T)')
        ylabel('Fourier amp (arb. units)')
        cf = (mean(obj.FFT.range(ii).upTemp(jj).range));
   end
   legend(leg)
   title(sprintf('Center field %0.1f, 0.015 1/T fourier range', cf))
end
%%
close all
for ii = 1:length(obj.FFT.range)
%     figure 
%     plot(obj.FFT.range(ii).upTemp(3).f,obj.FFT.range(ii).upTemp(3).FFT)
    cf(ii) = mean(obj.FFT.range(ii).upTemp(3).range);
    rngV(ii,:) = obj.FFT.range(ii).upTemp(3).range;
    for jj = 1:length(obj.mass.range(1).upPeak)
        jj
%         figure(1)
        
        A(jj,ii) = obj.mass.range(ii).upPeak(jj).A(3)
        mf(jj,ii) = mean(obj.mass.range(ii).upPeak(jj).maxFreq)
    end
 
end
%%
figure
txt = []
for ii = 1:length(A(:,1))
    plot(cf,A(ii,:),'-*')
    hold on
    txt = [txt,{num2str(mean(mf(ii,:)))}];
    pause(1)
    
end
legend(txt)
%% Plot max peaks vs. temps
clearvars leg
figure
for jj = 1:length(obj.mass.range)
    %%
%     subplot(2,3,jj) 
%     figure
%     jj = 8;
    for ii = 4%:length(obj.mass.range(jj).upPeak)
        plot(obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).A,'-*')
        hold on
        leg{ii} = strcat('f = ', num2str(round(mean(obj.mass.range(jj).upPeak(ii).maxFreq))))
    end
    xlabel('T (K)')
    ylabel('Amplitude (arb. units)')
    title(num2str(obj.FFT.range(jj).upTemp(1).range))
    
end
legend(leg)

%%
figure 
leg = [];
rng = 4
for ii = 1:length(obj.mass.range(rng).upPeak)
    
    mF(ii) = obj.mass.range(rng).upPeak(ii).maxFreq(3)
    
    for jj = 1:length(obj.FFT.range(rng).upTemp)
        dif = abs(obj.FFT.range.upTemp(jj).f - mF(ii));
%         figure
%         plot(dif)
        [~,phI] = min(dif);
        chphI = obj.FFT.range.upTemp(jj).f(phI)
        tempV(jj) = obj.FFT.range.upTemp(jj).temp;
        phaseV(jj) = obj.FFT.range.upTemp(jj).phase(phI);
    end
    
    
    plot(tempV,phaseV)
    hold on
    txt = {strcat('Peak ', num2str(ii))}
    leg = [leg,txt];
                       
end
legend(leg)
xlabel('Temp')
ylabel('phase')
tlte('phase at fixed frequency')

%% plot dM/dH vs 1/H. Fig 1 in document created using this section with
% endfield  == 40 and iFFspan=abs((1/10-1/40))=0.075.
temp = 1:length(obj.raw);
for jj = length(obj.FFT.range)
    for ii = temp
        
        figure(ii)
        subplot(2,1,1)
        plot(obj.FFT.range(jj).upTemp(ii).x,obj.FFT.range(jj).upTemp(ii).y,'Linewidth',1.2)
        hold on
        plot(obj.FFT.range(jj).upTemp(ii).xspl,obj.FFT.range(jj).upTemp(ii).yspl,'Linewidth',.4)
        xlabel('1/H')
        ylabel('dM/dH')
        legend('dM/dH vs 1/H','Spline fit')
        subplot(2,1,2)
        plot(1./obj.FFT.range(jj).upTemp(ii).x,obj.FFT.range(jj).upTemp(ii).y,'Linewidth',1)
        hold on
        plot(1./obj.FFT.range(jj).upTemp(ii).xspl,obj.FFT.range(jj).upTemp(ii).yspl,'Linewidth',.4)
        xlabel('H')
        ylabel('dM/dH')
        legend('dM/dH vs 1/H','Spline fit')
    end      
end
suptitl1 = strcat(num2str(obj.raw(ii).temp),'K');
suptitle(strcat('Fig 1: T=',suptitl1))

%% plot Fourier transform for all temps 
close all
temp = 1:12;%[2 3 4 5 6 7 10];
sp = 1;

leg = [];
figure 
for jj = 1:length(obj.FFT.range)
%     subplot(2,3,sp)

      
%     subplot(2,3,jj)
    figure
    for ii = temp%1:length(obj.raw)
%         figure(jj)
%         leg{jj} = strcat(CFnum2str(round(mean(obj.FFT.range(jj).upTemp(ii).range))),'T');
        obj.FFT.range(jj).upTemp(ii).range
%         figure(ii)
%         txtT = {sprintf('CF = %.0f T',round(mean(obj.FFT.range(jj).upTemp(ii).range)))};
        txt = {strcat('Temp ', num2str(obj.FFT.range(jj).upTemp(ii).temp))};
        leg = [leg,txt];
%         subplot(2,4,sp)
        plot(obj.FFT.range(jj).upTemp(ii).f,obj.FFT.range(jj).upTemp(ii).FFT) 
        title(strcat('CF = ',num2str(mean(obj.FFT.range(jj).upTemp(ii).range))))%txtT)
        hold on
%         pause(1)
        
    end

sp = sp+1;
% titl1 = strcat('CF = ',num2str(round(mean(obj.FFT.range(jj).upTemp(1).range))));
% titl = strcat(titl1,'T');
% txt = sprintf('1/H window = %0.1d (1/T)',iFFspan)
% title(txt)  
% title(titl)
xlabel('frequency (T)')
ylabel('Amplitude (arb. units)')
legend(leg)

end

% Used to creat fig.5a-c (fig 5c has a different lable but is the third figure) 
% allthought two different ranges was scanned over. 
%% Plot max peaks vs. temps
clearvars leg
figure
for jj = 1:length(obj.mass.range)
    %%
    subplot(2,3,jj) 
%     jj = 8;
    for ii = 1:length(obj.mass.range(jj).upPeak)
        plot(obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).A,'-*')
        hold on
        leg{ii} = strcat('f = ', num2str(round(mean(obj.mass.range(jj).upPeak(ii).maxFreq))));
    end
    xlabel('T (K)')
    ylabel('Amplitude (arb. units)')
    title(num2str(obj.FFT.range(jj).upTemp(1).range))
    
end
legend(leg)

%% Plot AoT vs. temp from my MATLAB fft
figure
for jj = 1:length(obj.FFT.range)
    %%
%     figure
%     jj = 8;
    subplot(2,3,jj)
    for ii = 1:length(obj.mass.range(1).upPeak)
        
%         subplot(2,2,ii)
%         figure
        plot(obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).AoT,'*',...
            obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).AoTfitV)
%         maxFreq = obj.mass.range.upPeak(kk).maxFreq;
%         title(num2str(mean(maxFreq)))
        xlabel('Temperature (K)')
        ylabel('Amplitude/T (arb. units)')
        
    end
    titl = num2str(mean(obj.FFT.range(jj).upTemp(1).range));
    suptitle(titl)
end
endCnt = 23
%% Plot mass vs center field and frequency vs center field
% figure
if ~exist('endCnt')
    endCnt = 1;
end
% chendCnt = endCnt
for ii = 1:length(obj.FFT.range)
    CF(ii) = mean(obj.FFT.range(ii).upTemp(2).range);
    for jj = 1:length(obj.mass.range(ii).upPeak)
%         chjj = jj 
        massCF(jj,ii) = obj.mass.range(ii).upPeak(jj).m;
        errorCF(jj,ii) = obj.mass.range(ii).upPeak(jj).AoTrms;
        MyErrorCF(jj,ii) = sum((obj.mass.range(ii).upPeak(jj).AoT-obj.mass.range(ii).upPeak(jj).AoTfitV).^2);
        maxPeakSTD(jj,ii) = std(obj.mass.range(ii).upPeak(jj).maxFreq);
        maxPeakCF(jj,ii) = mean(obj.mass.range(ii).upPeak(jj).maxFreq);
        maxPhaseCF(jj,ii) = mean(obj.mass.range(ii).upPeak(jj).phase);
        maxA1(:,ii) = obj.mass.range(ii).upPeak(1).A;
        maxA2(:,ii) = obj.mass.range(ii).upPeak(2).A;
        maxA3(:,ii) = obj.mass.range(ii).upPeak(3).A;
        maxA4(:,ii) = obj.mass.range(ii).upPeak(4).A;
        
    end
end
% 
% figure(endCnt)
% figure(endCnt+1)
% figure(endCnt+2)

for kk = 1:2%:length(maxPeakCF(:,1))
%     figure;%(endCnt)
%     errorbar(CF,maxPeakCF(kk,:),maxPeakSTD(kk,:),'-*')
%     title(num2str(kk))
%     hold on
%     
    figure(1)%(endCnt+1)
    plot(CF,massCF(kk,:),'-*')
%     subplot(2,1,1)
%     errorbar(CF,massCF(kk,:),errorCF(kk,:),'-*')
    xlabel('Center Field (T)')
    ylabel('effective mass')
%     subplot(2,1,2)
%     plot(CF,massCF(kk,:),'-*')
%     xlabel('Center Field (T)')
%     ylabel('effective mass')
    suptitle(strcat('Peak ',num2str(kk)))
    hold on
%     
%     figure(endCnt+2)
%     plot(CF,errorCF(kk,:),'-*',CF,MyErrorCF(kk,:),'o')
%     hold on
%     
%     figure(endCnt+3)
%     subplot(2,2,kk)
%     plot(CF,maxPhaseCF(kk,:),'*-')
%     hold on
% % 
%     figure
%     subplot(2,1,1)
%     plot(CF,errorCF(kk,:),'-*')
%     title('Sum of Squares Due to Error (SSE) vs CF')
%     xlabel('CF')
%     ylabel('SSE')
%     subplot(2,1,2)
%     plot(CF,massCF(kk,:),'-*')
%     title('effective mass vs cf')
%     xlabel('CF')
%     ylabel('m*')
%     suptitle(strcat(num2str(kk),' peak, range .0068 1/T'))
    
end
%%
figure(1)%(endCnt)
title('center field vs max peaks')
legend('Peak 1','Peak 2','Peak 4','Peak 5')
xlabel('Center Field (T)')
ylabel('Peak Frequency (T)')

figure(endCnt+1)
title('Fig. 6.b m* vs center field')
legend('Peak 1','Peak 2','Peak 4','Peak 5')
xlabel('Center Field (T)')
ylabel('effective mass')

figure(endCnt+2)
title('mass error vs center field')
legend('Peak 1','Peak 2','Peak 4','Peak 5')
xlabel('Center field (T)')
ylabel('mass sse')

figure(endCnt+3)
suptitle('phase vs center field')
legend('Peak 1','Peak 2','Peak 4','Peak 5')
xlabel('Center field (T)')
ylabel('phase')

%%
for jj = 1:4
    mp(jj,:) = maxPhaseCF(jj,:);
    for ii = 1:length(maxPhaseCF(1,:))-1

        dif = mp(jj,ii)-mp(jj,ii+1);
        cnt = 0;
        if abs(dif) >= pi/2
            cnt = cnt+1
            mp(jj,ii+1) = mp(jj,ii+1)+sign(dif)*pi/2;

        end
    end

    figure
    plot(CF,mp(jj,:),'-*r',CF,maxPhaseCF(jj,:),'-o')
    title(strcat('peak ',num2str(jj)))
    legend('corected phase', 'origional phase')
    xlabel('Center field (T)')
    ylabel('phase')
end
%%
figure
clear leg
% leg = []
for ii = 1:12
%     figure
    leg{ii} = strcat(num2str(obj.mass.range(1).upPeak(4).temp(ii)),' K')
    plot(CF,maxA4(ii,:))
    hold on
    xlabel('cf')
    ylabel('max amplitude')
%     title(num2str(ii))
    
end
legend(leg,'Location','northeastoutside')

% Used to creat Fig 6 allthought two different ranges was scanned over.