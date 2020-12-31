function [] = saveAoT1(obj,path,name)
    for ii = 1:length(obj.mass.range(1).upPeak)
        ii
        for jj = 1:length(obj.mass.range)
            AoTar(:,jj+1) = obj.mass.range(jj).upPeak(ii).AoT;
            mass(jj,1) = obj.mass.range(jj).upPeak(ii).m;
            maxFreq(jj,1) = mean(obj.mass.range(jj).upPeak(ii).maxFreq);
            cfar(:,jj) = obj.FFT.range(jj).upTemp(1).range';
            
        end
        
        AoTar(:,1) = [obj.mass.range(ii).upPeak(1).temp];
        AoTtab = array2table(AoTar);
        AoTtab.Properties.VariableNames(1) = {'temp'};
        
        sFile = sprintf('%sAoT%s%0.0f.csv',path,name,ii);
        writetable(AoTtab,sFile,'Delimiter',',','QuoteStrings',true)%'AoTfull15to35.csv'        
        

    end
    
    cftab = array2table(cfar);
    sFile2 = sprintf('%scf.csv',path,name);
    writetable(cftab,sFile2,'Delimiter',',','QuoteStrings',true)
    
end