clear all;
close all;

Ini_Bms = [12:2:50];
equ_dB = 1;
dB = 10;
% Ini_Bms = [20:2:40];
% equ_dB = 2;
% dB = 1/100;

range = [4300 4450]; % find peak range
drange = 10;

mapomT = [];
filename_Fdep = ['C:\Users\17132\Google Drive\NbGe2\oscillation data\data_analysis\tmp\dataFdep_long_axis_',num2str(Ini_Bms(1)),'_',num2str((Ini_Bms(end)-Ini_Bms(1))/length(Ini_Bms)),'_',num2str(Ini_Bms(end)),'.mat'];
if(isfile(filename_Fdep))
    load(filename_Fdep);
else
    dataFdep.nominalH = [];
    dataFdep.H = [];
    dataFdep.mass = [];
    dataFdep.masserr = [];
    dataFdep.pha = [];
    dataFdep.phaerr = [];
    dataFdep.loc = [];
    dataFdep.locerr = [];
    dataFdep.mT = {};
    for i = 1:length(Ini_Bms)
        Ini_Bm = Ini_Bms(i);
        analysis_ParaToLongAxis_field_dependence_section;
        dataFdep.nominalH(end+1) = Ini_Bm;
        dataFdep.H(end+1) = Bm;
        dataFdep.mass(end+1) = pm.ms(pks_i);
        dataFdep.masserr(end+1) = pm.ms_err(pks_i);
        dataFdep.pha(end+1) = pm.pha(pks_i);
        dataFdep.phaerr(end+1) = pm.pha_err(pks_i);
        dataFdep.loc(end+1) = pm.loc(pks_i);
        dataFdep.locerr(end+1) = pm.loc_err(pks_i);
%         mapomT(:,end+1) = pm.mT;
    end
end
%%
figure();
subplot(1,3,1);
errorbar(dataFdep.H,abs(dataFdep.mass),dataFdep.masserr,'or');
ylim([0 3]);
ylabel('Effective Mass');
subplot(1,3,2);
errorbar(dataFdep.H,dataFdep.pha,dataFdep.phaerr,'or');
ylim([-pi pi]);
ylabel('Phase (pi)');
subplot(1,3,3);
errorbar(dataFdep.H,dataFdep.loc,dataFdep.locerr,'or');
% ylim([-0 1]);
ylabel('Fre (T)')
save(filename_Fdep,'dataFdep');

% %%
% figure();
% [xx,yy] = meshgrid(T_arr,dataFdep.H);
% fp = surf(xx,yy,mapomT');
% set(fp,'EdgeColor','none');
% % shading interp;
% caxis([0  2]);
% % ylim([0 3]);
% % xlim([20 100]);
% colorbar();
% xlabel('T (K)');
% ylabel('H (T)');
% title('map of mT');
% 
