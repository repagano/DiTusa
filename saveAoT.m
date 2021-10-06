function [sFile] = saveAoT(obj,path,name)
    for ii = 1:length(obj.FFT.range)
        ii
        for jj = 1:length(obj.mass.range(ii).upPeak)
            AoTar(:,jj+1) = obj.mass.range(ii).upPeak(jj).AoT;
            mass(jj,1) = obj.mass.range(ii).upPeak(jj).m;
            maxFreq(jj,1) = mean(obj.mass.range(ii).upPeak(jj).maxFreq);
            
        end
        
        AoTar(:,1) = [obj.mass.range(ii).upPeak(1).temp];
        AoTtab = array2table(AoTar);
        AoTtab.Properties.VariableNames(1) = {'temp'};
        
        sFile = sprintf('%sAoT%s.csv',path,name)
        writetable(AoTtab,sFile,'Delimiter',',','QuoteStrings',true)%'AoTfull15to35.csv'
        
        massarray = [mass,maxFreq];
        masstab = array2table(massarray)
        masstab.Properties.VariableNames(2) = {'maxFreq'};
        masstab.Properties.VariableNames(1) = {'massMatlabCalc'};
        masstab
        
        sFile2 = sprintf('%sMaxFreq_MATLABmass%s.csv',path,name)
        writetable(masstab,sFile2,'Delimiter',',','QuoteStrings',true)%'AoTfull15to35.csv'
%         cf = obj.FFT.range(ii).upTemp(1).range
        
%         obj.range(ii).range
%         mp = mean(obj.mass.range(ii).upPeak(jj).maxFreq)

        
    end
end