% Git stop being a little bitch and work!! 
clear obj endField
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
num = 1
endFields = 55%20:5:55;
iFFspan = abs((1/25-1/55))/num;%Here the width of the 1/H window is defined 
% 
% startFields = 15
% for ii = 1:num
%     ii
%     endFields(ii) = 1/(-iFFspan*ii+1/startFields)
%     
% end

%% Here the Fourier transform is performed  
dHvA.FFTload(obj,endFields,iFFspan,[300 4500]);

%% Here the effective masses of the different peaks are calculated 
% peakRange is an matrix whose rows define the range around each peak 
% identified on the Fourier spectrum. The columns correspond to different
% peaks. The range should be large enough to allow for variation in the position 
% of the peak as the different ranges and temperatures are scanned over, but 
% not be so large that another peak is selected. 
peakRange = [490 708; 708 850; 2600 2900;2901 3545; 3546 3880]; %850 980;
% dHvA.massLoad runs the method that calculates the effective masses.
dHvA.massLoad(obj,peakRange);


%% Fourier transform as a function of cftemp for mult cfs
close all
tempv = obj.mass.range(1).upPeak(1).temp;
for ii = 1:length(obj.FFT.range)
    CF = mean(obj.FFT.range(ii).upTemp(1).range);
    ttext = sprintf('CF = %.2g T',CF);
%     figure
    c = hsv(length(tempv));
    leg = [];
    for jj = 1:length(tempv)
        
%         plot(obj.FFT.range(ii).upTemp(jj).f,obj.FFT.range(ii).upTemp(jj).FFT,'LineWidth',1.25,'Color',c(jj,:))
%         hold on
%         pause(.5)
        ltext = strcat(num2str(obj.FFT.range(ii).upTemp(jj).temp),' K');
        leg = [leg,{ltext}];%{num2str(obj.FFT.range(ii).upTemp(1).range)}];
        legend(leg)
        xlabel('Frequency (T)')
        ylabel('Fourier amp. (arb. units)')

    %     p = obj.mass.range(ii).upPeak( 
    end
    title(ttext)
end
%%
% figure(1)
% leg2 = [];
% titv = [{'610 T Peak'},{'795 T Peal'},{'2767 T Peak'},{'3292 T Peak'},{'3700 T Peak'}];
for ii = 1:5
%     figure
    leg1 = []
    for jj = 1:length(obj.mass.range)
        cf(jj) = mean(obj.FFT.range(jj).upTemp(1).range)
        A(jj) = (obj.mass.range(jj).upPeak(ii).A(3))
        
%         plot(obj.mass.range(jj).upPeak(ii).temp,obj.mass.range(jj).upPeak(ii).A,'LineWidth',1)
%         hold on
%         ltext = sprintf('CF = %.2g T',cf(jj))
%         leg1 = [leg1,{ltext}]
%         title(titv(ii))
%         xlabel('Temp (K)')
%         ylabel('Fourier amp. (arb. units)')
    end
%     legend(leg1)
%         
%     figure(1)
%     plot(cf,A','LineWidth',1)
%     hold on
%     leg2 = [leg2,{num2str(ii)}] 
end
% legend(titv)
% xlabel('CF (T)')
% ylabel('Fourier amp. (arb. units)')
% title('1.52 T')
%% Fourier transform as a function of cf for mult temps
% tempv = obj.mass.range(1).upPeak(1).temp;
% for jj = 1:length(tempv)
%     
%     figure
%     leg = []
%     for ii = 1:length(obj.FFT.range)
%         CF = mean(obj.FFT.range(ii).upTemp(1).range);        
%         plot(obj.FFT.range(ii).upTemp(jj).f,obj.FFT.range(ii).upTemp(jj).FFT,'LineWidth',1.25)
%         hold on
%         pause(.5)
%         ltext = sprintf('CF = %.2g T',CF)
%         leg = [leg,{ltext}]%{num2str(obj.FFT.range(ii).upTemp(1).range)}];
%         legend(leg)
%         tv = strcat(num2str(obj.FFT.range(ii).upTemp(jj).temp),' K')
%         title(tv)
%         xlabel('Frequency (T)')
%         ylabel('Fourier amp. (arb. units)')
% 
%     %     p = obj.mass.range(ii).upPeak( 
%     end
% end


%%
% rangeI = 1;
% for ii = 1:length(obj.FFT.range(1).upTemp)
%    
%     ff(ii) = obj.FFT.range(rangeI).upTemp(ii).xspl(1);
%     ffspl(ii) = obj.FFT.range(rangeI).upTemp(ii).x(1);
%     
% end
% 
% run = [1:length(obj.FFT.range(1).upTemp)];
% figure
% plot(run,ff,'*r',run,ffspl,'*r')
% hold on
%%
clc
clear phaseV difP phase signV2 r
% close all
% figure
for rI = 1:length(obj.FFT.range)
%     close all
    for tI = 1:length(obj.FFT.range(rI).upTemp)
        
%         This Line !!!
        phaseV = (obj.FFT.range(rI).upTemp(tI).phase);%unwrap

%         phaseV(1) = phaseV(1)+pi;
%         phaseV1 = phaseV;
        for pii = 1:length(obj.mass.range(rI).upPeak)
            chpii = pii
            maxPhase(tI,pii) = obj.mass.range(rI).upPeak(pii).phase(tI);
            maxf(tI,pii) = obj.mass.range(rI).upPeak(pii).maxFreq(tI);
        end
        
%         plot(obj.FFT.range(rI).upTemp(tI).f,phaseV,maxf(tI,:),maxPhase(tI,:),'*')%obj.FFT.range(rI).upTemp(tI).phase
%         xlabel('frequency (T)')
%         ylabel('phase')
%         title('hanning window')
%         hold on
%         pause(1)
        
        r(rI).temp(tI).phase = phaseV;
        r(rI).temp(tI).f = obj.FFT.range(rI).upTemp(tI).f;
        clear phaseV %maxPhase maxf
%         for fI = 1:length(obj.FFT.range(rI).upTemp(1).phase)-1
% 
%             difP(fI) = (diff([phaseV(fI) phaseV(fI+1)]));
%             signV2(fI) = sign(difP(fI));
%             cnt = 0;
%             if abs(difP(fI)) > pi/2 && abs(difP(fI))< 3*pi/2
%                 cnt = cnt+1;
%                 phaseV(fI+1) = phaseV(fI+1)-signV2(fI)*pi;
%                 elseif abs(difP(fI)) >3*pi/2
%                     phaseV(fI+1) = phaseV(fI+1)-signV2(fI)*2*pi;
%             else 
%                 phaseV(fI+1) = phaseV(fI+1);
%             end
%             
% %             figure(1)
% %             plot(obj.FFT.range(rI).upTemp(tI).f,phaseV)
% %             pause(.0001)
% 
%         end

%         NEXT TWO LINES!!!  
%         r(rI).temp(tI).phase = phaseV;
%         r(rI).temp(tI).f = obj.FFT.range(rI).upTemp(tI).f;
        
%         r(rI).temp(tI).maxPhase = [obj.mass.range(rI).upPeak(1).phase(tI),...
%             obj.mass.range(rI).upPeak(2).phase(tI), obj.mass.range(rI).upPeak(3).phase(tI),...
%             obj.mass.range(rI).upPeak(4).phase(tI)];
%         r(rI).temp(tI).maxf = [obj.mass.range(rI).upPeak(1).maxFreq(tI),...
%             obj.mass.range(rI).upPeak(2).maxFreq(tI), obj.mass.range(rI).upPeak(3).maxFreq(tI),...
%             obj.mass.range(rI).upPeak(4).maxFreq(tI)];
        

%         figure(1)
%         plot(obj.FFT.range(rI).upTemp(tI).f,r(rI).temp(tI).phase)
%         hold on
%         plot(r(rI).temp(tI).maxf,r(rI).temp(tI).maxPhase)
    
    end
end

%%
clear phaseV phaseAvg dif1 dif2 signV CF fI rI shiftshi
clc
% close all
rangnum = length(obj.mass.range);%1;%
phaseV1 = [];
phaseV2 = [];
phaseV3 = [];
phaseV4 = [];
cfLeg = [];

for rI = 1:rangnum
    chhh = rI;
%     figure 
    legp2 = [];
    for pI = 1:length(obj.mass.range(rI).upPeak)
        chii = pI;
        mF(pI) = obj.mass.range(rI).upPeak(pI).maxFreq(3);
                
        for tI = 1:7%length(obj.FFT.range(hh).upTemp)
            
            chjj = tI;
            dif = abs(obj.FFT.range(rI).upTemp(tI).f - mF(pI));
    %         figure
    %         plot(dif)
            [~,phI] = min(dif);
%             chphI = obj.FFT.range(hh).upTemp(jj).f(phI)
            tempV(tI) = obj.FFT.range(rI).upTemp(tI).temp;
            phaseV(pI,tI) = r(rI).temp(tI).phase(phI);
            fV(pI,tI) = obj.FFT.range(rI).upTemp(tI).f(phI);
%             phaseV(fI,jj) = obj.FFT.range(hh).upTemp(jj).phase(phI);
            
%             if ii == 2
%                 phaseV(ii,jj) = phaseV(ii,jj)-2*pi
%             end

%             chhij1 = [rI,pI,tI];
            for kk = 1:length(phaseV(1,:))-1
            if tI>1
                
                dif1(pI,tI-1) = diff([phaseV(pI,tI-1) phaseV(pI,tI)]);
                signV2(pI,tI-1) = sign(dif1(pI,tI-1));
                cnt = 0;
                if abs(dif1(pI,tI-1)) >= 2*pi*.70
%                     cnt = cnt+1
                    phaseV(pI,tI) = phaseV(pI,tI)-signV2(pI,tI-1)*pi*2;
                end
            end
            end

        end
        
        chphasVr = phaseV(pI,:);
        mv = mean(phaseV(pI,:));
%         chphaseV = phaseV
        phaseAvg(pI,rI) = mean(phaseV(pI,:));
       
        
%         
        
        if rI > 1  
%             hich = [ii hh]
            difm = [phaseAvg(pI,rI-1) phaseAvg(pI,rI)];
            difV = diff([phaseAvg(pI,rI-1) phaseAvg(pI,rI)]);
            difAvgV(pI,rI-1) = difV;
            dif2(pI,rI-1) = difV;%diff([phaseAvg(ii,hh-1) phaseAvg(ii,hh)]);
            signV2(pI,rI-1) = sign(difV);
            cnt = 0;
            
            shift(pI,rI-1) = 0;
            if abs(dif2(pI,rI-1)) >= 2*pi*.65
                shift(pI,rI) = 1;
                phaseAvg(pI,rI) = phaseAvg(pI,rI)-signV2(pI,rI-1)*pi*2;
                
            elseif abs(dif2(pI,rI-1)) >= pi*.65
                shift(pI,rI) = 1;
                phaseAvg(pI,rI) = phaseAvg(pI,rI)-signV2(pI,rI-1)*pi;
                
            end
        end
%         pause(1)
        CF(rI)= mean(obj.FFT.range(rI).upTemp(1).range);
        
%         figure(rI)
% %         subplot(ceil(rangnum/2),2,hh)
%         xlabel('Temp')
%         ylabel('phase')
% %         subplot(2,2,hh)
%         plot(tempV,phaseV(pI,:),'*-')
%         hold on
        txt = {strcat('Peak ', num2str(pI))};
        legp2 = [legp2,txt];
        
%         title(strcat('1/CF =', num2str(1/CF(rI))))
%         legend(legp2)
        
%         figure(rI+rangnum)
%         chhi2 = [rI pI]
        switch pI
            case 1
                phaseV1(rI,:) = phaseV(pI,:);
                cfLeg = [cfLeg,{strcat('CF = ',num2str(CF(rI)))}];
            case 2
                phaseV2(rI,:) = phaseV(pI,:);
            case 3 
                phaseV3(rI,:) = phaseV(pI,:);
            case 4
                phaseV4(rI,:) = phaseV(pI,:);

        end
        

    end
    

    
%     figure(1)
%     plot(fV(:,1),phaseV(:,1),'*',r(rI).temp(tI).f,r(rI).temp(tI).phase)
%     hold on
    
%     clear phasV
end

phaseAvg(1,:) = phaseAvg(1,:)+pi;

% figure(1)
% plot(CF,phaseAvg,'-*')
% xlabel('CF')
% ylabel('Phase avraged over temp')
% legend(legp2)
% title('30 to 55 to 25')


%%
clear rng ii
rng = 6
rv(rng).cf = CF;

for ii = 1:length(phaseAvg(:,1))
    ii
    rv(rng).peak(ii).pa = phaseAvg(ii,:);
   
end
% legend('Test 1')


%%
% titlV = [{'Peak 610 T, Fig 4.a'},{'Peak 795 T, Fig 4.b'},{'Peak 2750 T, Fig 4.c'},{'Peak 3700 T, Fig 4.d'}]
% figure
% for jj = 1:rng
%     for ii = 1:4
%         subplot(2,2,ii)
%         plot(rv(jj).cf(1:end),rv(jj).peak(ii).pa(1:end),'*-')
%         hold on
%         xlabel('Center field (T)')
%         ylabel('Avg. Phase')
%         title(titlV(ii))
% %         pause(1)
%     end
% end
% 
% legend('rng .049 1/T','rng .032 1/T','rng .022 1/T','rng .015 1/T','rng .0104 1/T','rng .0062 1/T')

   


%%
% % temp = obj.mass.range(1).upPeak(1).temp
% % maxf_Hann = maxf
% % maxPhase_range = maxPhase
% % 
% % figure
% % plot(temp,maxPhase')
% % % plot(temp,maxPhase_HannAzero(:,1),'*-',temp,maxPhase_Hann(:,1),'*-',temp,maxPhase_raw(:,1),'*-')
% % xlabel('temp')
% % ylabel('phase')
% % legend('zero padding and hanning','hanning','raw data')
% % title('peak 610 T')
% % % figure
% % % plot(obj.FFT.range(1).upTemp(1).f, obj.FFT.range(1).upTemp(1).phase)
% % % 
% % % figure
% % % plot(obj.FFT.range(1).upTemp(1).f,phaseV1)
% % % 
% % % figure
% % % plot(obj.FFT.range(1).upTemp(1).f(1:end-1),difP)
% % 
% % 
% % %%
% % for ii = 1:length(obj.FFT.range)
% %     tng = obj.FFT.range(ii).upTemp(1).range
% % end
