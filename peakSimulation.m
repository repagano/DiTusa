close all
clear 
clc
%%
for ii = 1:length(rngV(:,1))
    H = [rngV(ii,1):.0001:rngV(ii,2)];
    xi = 2.2*exp(-17./H).*cos(2*pi*3768./H+pi*3786/2)
    objFab.xUp = H;
    objFab.yUp = xi;
    objFab.temp = 1;
    
end