function [mass] = massCalc(FFTobj,peakRange,pCh)
% massCalc takes the peak Fourier amplitudes from multiple temperture runs
% and calculated the effective mass corisponding to each peak

% FFTobj(FFTobj,peakRange,'off') inputs multiple structs contained in 
% FFTobj with frequency data(FFTobj.f), Fourier amplitude data
% (FFTobj.FFT), temperature data (FFTobj.temp), and range data 
% (FFTobj.range). peakRange is an N by 2 array where N is the number of 
% peaks The first and second columb identify the begining and end of the
% range of interest in the data corisponding to each peak.

% FFTobj(FFTobj,peakRange,'on') must include the phase (FFTobj.phase)
% produced by a Fourier transform as well as the data indicated above. The
% function will then extract the phase corisponding to the peak Fourier
% amplitude as well as calculating the effective mass. 
%     figure(1)
%     figure(2)
    for jj = 1:length(FFTobj)
%         chjj = jj
        f_og = FFTobj(jj).f;   
        fI = f_og>=peakRange(1) & f_og<=peakRange(2);
        f = f_og(fI);
        FFT_og = FFTobj(jj).FFT;
        FFTverror = mean(FFT_og);       
        
        FFTvs = FFTobj(jj).FFT(fI); 
        [A(jj,1),maxI] = max(FFTvs);
        
%         figure(1)
%         plot(f_og,FFTobj(jj).FFT)
%         hold on
%         plot(f,FFTvs)
%         yline(FFTverror)
%         plot(f(maxI),A(jj),'b*')
        
        weight(jj,1) = FFTverror/A(jj);
        temp(jj,1) = FFTobj(jj).temp;
        Brange = FFTobj(jj).range;
        MF = f(maxI);
        maxFreq(jj) = MF;
        
%         figure(2)
%         errorbar(MF,A(jj),weight(jj),'*')
%         hold on
        
        if strcmp(pCh,'on')
            phaseV = FFTobj(jj).phase(fI);
            phase = phaseV(maxI);
            if phase < 0 
                mass.phase(jj,1) = phase;
            else
                mass.phase(jj,1) = phase;
            end
        end

    end
    %AoT: amplitude over temperature 
    AoT = A./temp;
%     Weights = A.*.1./temp;
    
    %Fit AoT vs T data to dHvA mass equation
    Bm = 1/(1/2*(1/Brange(1)+1/Brange(end)));
    AoTfittype = fittype(@(a,b,x) a./sinh(b*x));%fittype(@(a,b,c,x) a./sinh(b*x)+c);%fittype(@(a,b,x) (a.*14.69*b/Bm)./sinh(14.69*b*x./Bm));
    AoTfitopt = fitoptions(AoTfittype);
    AoTfitopt.startpoint = [2,10];%[3.3 .29 0];
%     Weights = 
%     AoTfitopt.Weights = 1./(weight);%[.62 .96 1.52 2.24 3.27 3.96 5.71 8 10 12 14 18]./18;%ones(1,length(temp));%
%     AoTfittype = fittype(a*asinh(x*b),'coefficients',{'a','b'})
    bUpper = 3*14.69/Bm;
    [AoTfit, ~] = fit(temp,AoT,AoTfittype,AoTfitopt);%,'StartPoint',[0 .07]
    % define mass struct 
%     chAoTfit = AoTfit
    bErrorV = confint(AoTfit);
    bError =  diff([bErrorV(1,2) bErrorV(2,2)]);
    mass.AoTfit = AoTfit;
    AoTfitV = AoTfit(temp);%(1):.1:temp(end)
    m = abs(AoTfit.b)*Bm/14.69;
    mError = bError*Bm/14.69;

    mass.m = m;
    mass.A = A;
    mass.AoT = AoT;
    mass.AoTfitV = AoTfitV;
    mass.AoTrms = mError;%gof.sse;%*Bm/14.69;
%     mass.phase = phase;
    mass.temp = temp;
    mass.maxFreq = maxFreq;
        
end