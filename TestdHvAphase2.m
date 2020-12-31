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
                
        for tI = 1:length(obj.FFT.range(rI).upTemp)
            
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
        phaseAvg(pI,rI) = (phaseV(pI,3));
       
        
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
            if abs(dif2(pI,rI-1)) >= 2*pi*.60
                shift(pI,rI) = 1;
                phaseAvg(pI,rI) = phaseAvg(pI,rI)-signV2(pI,rI-1)*pi*2;
                
            elseif abs(dif2(pI,rI-1)) >= pi*.60
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
        
        title(strcat('1/CF =', num2str(1/CF(rI))))
        legend(legp2)
        
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



% figure
% plot(CF,phaseAvg,'-*')
% xlabel('CF')
% ylabel('Phase avraged over temp')
% legend(legp2)
% title('30 to 55 to 25')

%%
figure
plot([fV(1,3),fV(2,3),fV(4,3),fV(5,3)],[phaseAvg(1),phaseAvg(2),phaseAvg(4),phaseAvg(5)],'b*','MarkerSize',8)
hold on
plot([fV(1,3),fV(2,3),fV(4,3),fV(5,3)],[phaseAvg_1(1),phaseAvg_1(2),phaseAvg_1(4),phaseAvg_1(5)],'r*','MarkerSize',8)
% plot(fV(2,3),phaseAvg_1(2),'b*','MarkerSize',10)
% plot(fV(4,3),phaseAvg_1(4),'b*','MarkerSize',10)
% plot(fV(5,3),phaseAvg_1(5),'b*','MarkerSize',10)

% plot(fV(1,3),phaseAvg(1),'r*','MarkerSize',7)
% plot(fV(2,3),phaseAvg(2),'r*','MarkerSize',7)
% plot(fV(4,3),phaseAvg(4),'r*','MarkerSize',7)
% plot(fV(5,3),phaseAvg(5),'r*','MarkerSize',7)
ylabel('Phase')
xlabel('Frequency (T)')
legend('0.049 1/T range','0.038 1/T range')%
title('1.52 K temp data set')
%%
close all
figure
subplot(4,1,1)
plot(CF,phaseAvg(1,:),'-o')
title('610 T peak')
xlabel('temp')
ylabel('Phase')
subplot(4,1,2)
plot(CF,phaseAvg(2,:),'-o')
xlabel('temp')
ylabel('Phase')
title('795 T peak')
subplot(4,1,3)
plot(CF,phaseAvg(4,:),'-o')
xlabel('temp')
ylabel('Phase')
title('2750 T peak')
subplot(4,1,4)
plot(CF,phaseAvg(5,:),'-o')
xlabel('temp')
ylabel('Phase')
title('3700 T peak')
suptitle('1.52 K temp phase calc')

%%
figure
plot(CF,phaseV1(:,3),'-o')
xlabel('temp')
ylabel('Phase')
title('peak 610T')
% legend(cfLeg,'Location','northeastoutside')



figure
plot(tempV,phaseV2,'-o')
xlabel('temp')
ylabel('Phase')
title('peak 795T')
legend(cfLeg)

figure
plot(tempV,phaseV3,'-o')
xlabel('temp')
ylabel('Phase')
title('peak 2750T')
legend(cfLeg)

figure
plot(tempV,phaseV4,'-o')
xlabel('temp')
ylabel('Phase')
title('peak 3700T')
legend(cfLeg)

%%
% 
% % pAvg(
% figure
% leg =[];
% rII = 1%tI
% phaseV2 = phaseV;
% fV2 = fV;
% for tII = 1:tI
% %      plot(r(rII).temp(tII).f,r(rII).temp(tII).phase,fV(:,tII),phaseV(:,tII),'*')
% %     plot(r(rI).temp(tI).phase,r(rI).temp(tI).f,fV(:,tII),phaseV(:,tII),'*')%,'o',fV2(:,tII),phaseV2(:,tII),'*')
%     
%     hold on
% %     pause(4)
% %     leg = [leg,{strcat(num2str(tempV(tII)),' K')},{strcat(num2str(tempV(tII)),' K max peak val')}]
% end
% % legend(leg,'Location','eastoutside')
% 
% %     figure
% %     legend(leg)
% %     suptitle('phase at fixed frequency')
% % phaseAvg(5,:) = phaseAvg(5,:)-2*pi;
% % phaseAvg1 = phaseAvg; 
% % CF1 = CF;
% % figure
% % for ii = 1:length(phaseAvg(:,1))    
%     
% %     subplot(2,2,3)
% %     subplot(3,2,ii)
% % %     figure(ii)
% %     plot(1./CF,phaseAvg(ii,:),'*-')
% %     hold on
% %     plot(1./CF1,phaseAvg1(ii,:),'d-')
% %     xlabel('1/CF')
% %     ylabel('avg pahse')
% % end
% % legend('15 to 55 range split into fifths','15 to 55 range split into fourths')
% 
% 
% % figure
% % plot(CF)
% % 
% % figure
% % plot(CF1)
% % title('phase avraged over all temps vs center field')
% % xlabel('1/CF (1/T)')
% % ylabel('avg phase')
% % legend('peak 1','peak 2','peak 3','peak 4','peak 5')