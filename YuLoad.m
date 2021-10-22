function raw = YuLoad(sortedTemp,filename)
    
    dataArray = load(filename);
    FF = dataArray(:,1);
    xi = dataArray(:,2);
    
    [~,Ttop] = max(FF);
    Tbot = 1.5;               
    FFup = [FF(1:Ttop);zeros(length(FF(Ttop:end)),1)];
    
    FFupInd = FFup >= Tbot & FFup <= Ttop;
    raw.xUp = FFup(FFupInd);
    raw.yUp = xi(FFupInd);

    FFdown = [zeros(Ttop,1);FF(Ttop:end)];
    FFdownInd = find(FFdown >= Tbot & FFdown <= Ttop);
    raw.xDown = FF(FFdownInd(end):-1:FFdownInd(1));
    raw.yDown = xi(FFdownInd(end):-1:FFdownInd(1));                
    raw.temp = sortedTemp;
    raw.dataType = 'Yu';
    
    
end