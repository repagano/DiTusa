%% This is the code from TestdHvA
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

dHvA.FFTload(obj,endFields,iFFspan,[200 5e3],'extFFT');

%%
figure
plot(obj.FFT.range.upTemp(1).f(1:end),obj.FFT.range.upTemp(1).FFT)%(4450:end)
hold on
% x1 = [3919 4083 4141 4320 4394 4562 4607 4697];
% y1 = [43.81 41.39 35.86 33.78 19.98 18.92 7.738 20.09];
% plot(x1,y1,'*')
plot(obj.FFT.range.upTemp(2).f(1:end),obj.FFT.range.upTemp(2).FFT(1:end))
plot(obj.FFT.range.upTemp(3).f(1:end),obj.FFT.range.upTemp(3).FFT(1:end))
plot(obj.FFT.range.upTemp(4).f(1:end),obj.FFT.range.upTemp(4).FFT(1:end))
plot(obj.FFT.range.upTemp(5).f(1:end),obj.FFT.range.upTemp(5).FFT(1:end))
legend('0.62 K','0.96 K', '1.52 K','2.24 K','3.27 K')
% ylim([0 50])
%%
for ii = 2:length(x1)
%     chy1 = x1(ii-1)
%     chy2 = x1(ii)
    dif = (x1(ii-1)-x1(ii))
end

%%
clear ii pp0 ii L FFT f pointNum0
close all
pointNum0 = 0;
for hh = 1%:length(obj.FFT.range)
    figure(hh)
    leg = [];
    for ii = 1:length(obj.FFT.range(hh).upTemp)
%         chii1 = ii
        L = length(obj.FFT.range(hh).upTemp);
        FFT = obj.FFT.range(hh).upTemp(ii).FFT;
        f = obj.FFT.range(hh).upTemp(ii).f;
        minprom = std(FFT);
        [pp0(ii).range(hh).pks,pp0(ii).range(hh).locs,pp0(ii).range(hh).widths,pp0(ii).range(hh).proms] =...
            findpeaks(FFT,f,'MinPeakProminence',minprom);
        pp0(ii).range(hh).ind = true(1,length(pp0(ii).range(hh).pks));
        pp0(ii).range(hh).temp = obj.FFT.range(hh).upTemp(ii).temp;
        beginLocs0(ii) = (pp0(ii).range(hh).locs(1));
        hit = 0;
        cnt = 0;                                  
        
%         subplot(L,1,ii)
        plton = 1;
        if plton == 1
            p1 = plot(f,FFT);%,pp0(ii).locs,pp0(ii).pks,'*');
            hold on
%             set(get(get(p1(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            tval = strcat(num2str(obj.FFT.range(hh).upTemp(ii).temp),' K');
            leg = [leg,{tval}];%,strcat('p0 ',{num2str(ii)})];
            
        end
        
        pointNum0 = pointNum0+length(pp0(ii).range(hh).locs);
    end    
end

%%
clc
clear chdif beginLocs nnReove pp1 m pointNum1
for jj = 1:length(pp0(1).range)
    beginLocs = beginLocs0(:,hh);
    [minV,I] = min(beginLocs);
    nnV = 1:length(pp0);
    nnRemove = nnV~=I;
    nnV = nnV(nnRemove);
    mm = 1;
    pointNum1 = 0;

    while ~all(isnan(beginLocs))
    %     chmm(mm) = mm;
        V0 = find(pp0(I).range(jj).ind, 1, 'first');
        pp1(mm).range(jj).pks(1) = pp0(I).range(jj).pks(V0);
        pp1(mm).range(jj).locs(1) = pp0(I).range(jj).locs(V0);
        pp1(mm).range(jj).widths(1) = pp0(I).range(jj).widths(V0);
        pp1(mm).range(jj).temp(1) = pp0(I).range(jj).temp;
        pp0(I).range(jj).ind(V0) = false;
%         chdif(mm,I) = NaN;
%         chwidth(mm,I) = NaN;
        clear beginLocs
        if length(pp0(I).range(jj).widths)>V0
            beginLocs(I) = pp0(I).range(jj).locs(V0+1);
        else
            beginLocs(I) = NaN;
        end
        cnt = 2;
        for nn = nnV      
    %         chnn = nn
            V1 = find(pp0(nn).range(jj).ind, 1, 'first');
            chind = pp0(nn).range(jj).ind;
            if isempty(V1)
                beginLocs(nn) = NaN;
            else            
                dif = abs(pp0(nn).range(jj).locs(V1)-pp0(I).range(jj).locs(V0));
                halfWidth = pp0(I).range(jj).widths(V0)/2;
        %         chdif(mm,nn) = dif;
        %         chwidth(mm,nn) = pp0(I).widths(V0);
                if dif <= halfWidth
        %             disp('bigPP')
        %             chpp = pp0(I).pks(mm)
                    pp1(mm).range(jj).pks(cnt) = pp0(nn).range(jj).pks(V1);
                    pp1(mm).range(jj).locs(cnt) = pp0(nn).range(jj).locs(V1);
                    pp1(mm).range(jj).widths(cnt) = pp0(nn).range(jj).widths(V1);
                    pp1(mm).range(jj).temp(cnt) = pp0(nn).range(jj).temp;
%                     leg = [leg,{num2str(pp1(mm).range(jj).temp(cnt))}];
                    if length(pp0(nn).range(jj).widths)>V1
                        beginLocs(nn) = pp0(nn).range(jj).locs(V1+1);
                    else
                        beginLocs(nn) = NaN;
                    end
                    pp0(nn).range(jj).ind(V1) = false;
                    cnt = cnt+1;
                else %if pp0(nn).locs(V1)<pp0(I).locs(V0)-pp0(I).widths(V0)
                    beginLocs(nn) = pp0(nn).range(jj).locs(V1);
                end
            end
        end
        nlen = length(beginLocs);
        endlen = length(pp0);
%         chbeginLocs = beginLocs
        if nlen<endlen
            diflen = endlen-nlen;
            beginLocs(nlen+1:endlen) = ones(1,diflen)*nan;%zeros(1,diflen);%
        end
        nnRemove2 = beginLocs == 0;%isnan(beginLocs);
        beginLocs(nnRemove2) =  NaN;
        chbeginLocs = beginLocs
        [minV,I] = min(beginLocs) 
        nnV = 1:length(pp0);
        nnRemove1 = nnV~=I; %&& ~isnan(BeginLocs);
        nnRemove = nnRemove1 == true | nnRemove2 == true
        nnV = nnV(nnRemove)
        pointNum1 = pointNum1+length(pp1(mm).range(hh).locs);

        mm = mm+1;
    %     chpp1p = pp1(V0).pks
    %     chpp1w = pp1(V0).locs    

    end
end
    

%%
leg1 = [];
Ipks = []
cnt = 1;
for kk = 1:length(pp1(1).range)
    peakRange = [];
    for oo = 1:length(pp1)
        chpp1 = pp1(oo).range(kk).pks;
        figure(kk)
        plot(pp1(oo).range(kk).locs,pp1(oo).range(kk).pks,'*')
        hold on
%         leg1 = [leg1,{num2str(oo)}];

        if length(pp1(oo).range(kk).pks)>=4
            Ipks = [Ipks,oo]
            pp(cnt).range(kk).pks = pp1(oo).range(kk).pks;
            pp(cnt).range(kk).locs = pp1(oo).range(kk).locs;
            pp(cnt).range(kk).widths = pp1(oo).range(kk).locs;
            pp(cnt).range(kk).temp = pp1(oo).range(kk).temp;
            widthMean = mean(pp1(oo).range(kk).widths); 
            locMean = mean(pp1(oo).range(kk).locs);
            peakRange = [peakRange;locMean-widthMean/2 locMean+widthMean/2];            

            figure(1)
            plot(pp(cnt).range(kk).locs,pp(cnt).range(kk).pks,'o')
%             xline(locMean-widthMean/2)
%             xline(locMean+widthMean/2)
            cnt = cnt+1;
        end
    end
    rng(kk).peakRange = peakRange
end
% legend(leg)
xline(517,'r','LineWidth',1.5)
xlabel('Frequency (T)')
ylabel('Fourier amp. (arb. units)')
peakRange1 = [peakRange(1:5,:);3066 3165; 3252 3321; 3322 3381;peakRange(6:end,:)]
legend(leg)
%%
dHvA.massLoad(obj,peakRange(1:end,:));

%%
clc  
saveAoT(obj,'C:\Users\repag\Documents\MATLAB\AoTvals\AoTrng15to35\','full15to35')

%% Make table 
clear T maf
maf = [];
txt = [];
% figure
for ii = 1:length(obj.mass.range(1).upPeak)%length(obj.mass.range(1).upPeak)-
    objpeak = obj.mass.range.upPeak(ii);
%     for jj = 1:length(pp)
%         p.Tem(jj) = pp(jj).range.pks
    figure
    plot(objpeak.temp,objpeak.A,'*-')
%     hold on
%     plot(pp(ii).range.temp,pp(ii).range.pks,'o')
%     plot(pp(ii).
    txt = [txt,{sprintf('Peak frequency %.1f T',mean(objpeak.maxFreq))}];
    txt = [{sprintf('dHvA %.1f T',mean(objpeak.maxFreq))},{sprintf('peakFinder %.1f T',mean(pp(ii).range.locs))}];
    legend(txt)
    maf = [maf; [obj.mass.range(1).upPeak(ii).m, mean(obj.mass.range(1).upPeak(ii).maxFreq)]] 
end
legend(txt)
T = array2table(maf)
T.Properties.VariableNames = [{'mass'},{'peakFrequency'}]
%%
writetable(T,'C:\Users\repag\Documents\MATLAB\AoTvals\AoTrng15to55\massAndMaxfreq.csv','Delimiter',',','QuoteStrings',true)
%%

% for ii = 1:length(pp)
%     figure
%     plot(sort(temps),pp(ii).range.pks)
% end

%% This is the code from TestdHvArange it should be exactly the same
% clear ii pp0 ii L FFT f pointNum0
% close all
% pointNum0 = 0;
% for hh = 1:length(obj.FFT.range)
%     figure(hh)
%     leg = [];
%     for ii = 1:length(obj.FFT.range(hh).upTemp)
% %         chii1 = ii
%         L = length(obj.FFT.range(hh).upTemp);
%         FFT = obj.FFT.range(hh).upTemp(ii).FFT;
%         f = obj.FFT.range(hh).upTemp(ii).f;
%         minprom = std(FFT);
%         [pp0.range(hh).p(ii).pks,pp0.range(hh).p(ii).locs,pp0.range(hh).p(ii).widths,pp0.range(hh).p(ii).proms] =...
%             findpeaks(FFT,f,'MinPeakProminence',minprom);
%         pp0.range(hh).p(ii).ind = true(1,length(pp0.range(hh).p(ii).pks));
%         pp0.range(hh).p(ii).temp = obj.FFT.range(hh).upTemp(ii).temp;
%         beginLocs0(ii) = (pp0.range(hh).p(ii).locs(1));
%         hit = 0;
%         cnt = 0;                                  
%         
% %         subplot(L,1,ii)
%         plton = 1;
%         if plton == 1
%             p1 = plot(f,FFT);%,pp0(ii).locs,pp0(ii).pks,'*');
%             hold on
% %             set(get(get(p1(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%             leg = [leg,{num2str(ii)}];%,strcat('p0 ',{num2str(ii)})];
%             legend(leg)
%         end
%         
%         pointNum0 = pointNum0+length(pp0.range(hh).p(ii).locs);
%     end    
% end
% 
% %%
% 
% clc
% clear chdif beginLocs nnReove pp1 m pointNum1
% for jj = 1:length(pp0.range)
%     beginLocs = beginLocs0(:,jj);
%     [minV,I] = min(beginLocs);
%     nnV = 1:length(pp0.range(jj).p);
%     nnRemove = nnV~=I;
%     nnV = nnV(nnRemove);
%     mm = 1;
%     pointNum1(jj) = 0;
%     while ~all(isnan(beginLocs))
%     %     chmm(mm) = mm;
%         V0 = find(pp0.range(jj).p(I).ind, 1, 'first');
%         pp1.range(jj).p(mm).pks(1) = pp0.range(jj).p(I).pks(V0);
%         pp1.range(jj).p(mm).locs(1) = pp0.range(jj).p(I).locs(V0);
%         pp1.range(jj).p(mm).widths(1) = pp0.range(jj).p(I).widths(V0);
%         pp1.range(jj).p(mm).temp(1) = pp0.range(jj).p(I).temp;
%         pp0.range(jj).p(I).ind(V0) = false;
% %         chdif(mm,I) = NaN;
% %         chwidth(mm,I) = NaN;
%         clear beginLocs
%         if length(pp0.range(jj).p(I).widths)>V0
%             beginLocs(I) = pp0.range(jj).p(I).locs(V0+1);
%         else
%             beginLocs(I) = NaN;
%         end
%         cnt = 2;
%         for nn = nnV      
%     %         chnn = nn
%             V1 = find(pp0.range(jj).p(nn).ind, 1, 'first');
%             chind = pp0.range(jj).p(nn).ind;
%             if isempty(V1)
%                 beginLocs(nn) = NaN;
%             else            
%                 dif = abs(pp0.range(jj).p(nn).locs(V1)-pp0.range(jj).p(I).locs(V0));
%                 halfWidth = pp0.range(jj).p(I).widths(V0)/2;
%         %         chdif(mm,nn) = dif;
%         %         chwidth(mm,nn) = pp0(I).widths(V0);
%                 if dif <= halfWidth
%         %             disp('bigPP')
%         %             chpp = pp0(I).pks(mm)
%                     pp1.range(jj).p(mm).pks(cnt) = pp0.range(jj).p(nn).pks(V1);
%                     pp1.range(jj).p(mm).locs(cnt) = pp0.range(jj).p(nn).locs(V1);
%                     pp1.range(jj).p(mm).widths(cnt) = pp0.range(jj).p(nn).widths(V1);
%                     pp1.range(jj).p(mm).temp(cnt) = pp0.range(jj).p(nn).temp;
%                     if length(pp0.range(jj).p(nn).widths)>V1
%                         beginLocs(nn) = pp0.range(jj).p(nn).locs(V1+1);
%                     else
%                         beginLocs(nn) = NaN;
%                     end
%                     pp0.range(jj).p(nn).ind(V1) = false;
%                     cnt = cnt+1;
%                 else %if pp0(nn).locs(V1)<pp0(I).locs(V0)-pp0(I).widths(V0)
%                     beginLocs(nn) = pp0.range(jj).p(nn).locs(V1);
%                 end
%             end
%         end
%         nlen = length(beginLocs);
%         endlen = length(pp0.range(jj).p);
% %         chbeginLocs = beginLocs
%         if nlen<endlen
%             diflen = endlen-nlen;
%             beginLocs(nlen+1:endlen) = ones(1,diflen)*nan;%zeros(1,diflen);%
%         end
%         nnRemove2 = beginLocs == 0;%isnan(beginLocs);
%         beginLocs(nnRemove2) =  NaN;
%         chbeginLocs = beginLocs
%         [minV,I] = min(beginLocs) 
%         nnV = 1:length(pp0.range(jj).p);
%         nnRemove1 = nnV~=I; %&& ~isnan(BeginLocs);
%         nnRemove = nnRemove1 == true | nnRemove2 == true
%         nnV = nnV(nnRemove)
%         pointNum1(jj) = pointNum1(jj)+length(pp1.range(jj).p(mm).locs);
% 
%         mm = mm+1;
%     %     chpp1p = pp1(V0).pks
%     %     chpp1w = pp1(V0).locs    
% 
%     end
% end
%     
% 
% %%
% clear pp widthMean locMean 
% leg1 = [];
% cnt = 1;
% for kk = 1:length(pp1.range)
%     peakRange = [];
%     for oo = 1:length(pp1.range(kk).p)
%         chpp1 = pp1.range(kk).p(oo).pks;
%         figure(kk)
%         plot(pp1.range(kk).p(oo).locs,pp1.range(kk).p(oo).pks,'*')
%         hold on
%         leg1 = [leg1,{num2str(oo)}];
% 
%         if length(pp1.range(kk).p(oo).pks)>5
%             pp.range(kk).p(cnt).pks = pp1.range(kk).p(oo).pks;
%             pp.range(kk).p(cnt).locs = pp1.range(kk).p(oo).locs;
%             widthMean = mean(pp1.range(kk).p(oo).widths); 
%             locMean = mean(pp1.range(kk).p(oo).locs);
%             peakRange = [peakRange;locMean-widthMean/2 locMean+widthMean/2];
% 
%             figure(kk)
%             plot(pp.range(kk).p(cnt).locs,pp.range(kk).p(cnt).pks,'o')
%         end
%     end
%     rng(kk).peakRange = peakRange
%     clear oo
% end
% legend(leg1)