cb = get_colorbar;
Dir = 'C:\Users\17132\Google Drive\NbGe2\oscillation data\data_analysis\TempDep_data\';
tmpDir = 'C:\Users\17132\Google Drive\NbGe2\oscillation data\data_analysis\tmp\';
titlename = 'H // long - axis';
H_direction = 'a_long';
H_raw_d = 1;
H_raw_u = 60;
% Bm = 25;
% dB = 1/100;
if(equ_dB ~= 1)
    Bm = Ini_Bm;
    H_d = 1/(1/Bm+dB);
    H_u = 1/(1/Bm-dB);
else
    H_d = Ini_Bm - dB;
    H_u = Ini_Bm + dB;
    Bm = 1/(0.5*(1/H_d + 1/H_u));
end
disp(['H_d = ',num2str(H_d),'; H_u = ',num2str(H_u)]);
confi_level = 0.68;  % 1 signma,, 2 sigma 0.95

% for ID = 1:11, T = 
T_arr_init = [3.49, 2.67, 2.18, 2.38, 2.03, 1.77, 1.54, 1.55, 1.11, 3.98, 3.06 ];
start_of_fileID = 1;
end_of_fileID = 11;
T_arr = T_arr_init(start_of_fileID:end_of_fileID);

i = 0;
for fileID = start_of_fileID:end_of_fileID
    if(fileID < 10)
        filename = ['p00',num2str(fileID),'_012221H001.dat'];
    else
        filename = ['p0',num2str(fileID),'_012121H001.dat'];
    end
    i = i + 1;
    % load_single_file(filename,flat of decreasing field, field range, if_plot)
    % set field range [8,60].  when min< 8, there are bad points
    dat(i) = load_single_file([Dir,filename],1,[H_raw_d,H_raw_u],[H_d,H_u],0);   % decreasing field
%     dat(i) = load_single_file([Dir,filename],2,[H_raw_d,H_raw_u],[H_d,H_u],0);   % Increasing field
%     dat(i) = load_single_file([Dir,filename],3,[H_raw_d,H_raw_u],[H_d,H_u],0);  % both

end
disp('Reading datafiles completed ...');

T_arr_tmp = T_arr;
dat_tmp = dat;
[T_arr,tmp_index] = sort(T_arr_tmp,'descend');
dat = dat_tmp(tmp_index);

%% make plot of amplitude vs F with hann
FF = [];
num = length(T_arr);

% f_a = figure();
% f_p = figure();
cols_num = 3;
rows_num = ceil(num/3);
for i = 1:num
    [x, Ind]= unique(dat(i).invH);
    y = dat(i).F(Ind);

    L_signal = 1e6;
    x1 = x(1):(x(end)-x(1))/(L_signal-1):x(end);
    y1 = interp1(x,y,x1);
    tmp = isnan(y1);
    y1(tmp) = [];
    x1(tmp) = [];

    nfft = 2^(nextpow2(L_signal)+2);
    yFFT = fft(y1'.*hann(L_signal),nfft);
    yFFT = yFFT(1:nfft/2);
    pFFT = angle(yFFT);

    T = (x1(end)-x1(1))/length(x1);
    Fs = 1/T;
%     figure(f_a);
%     subplot(rows_num,cols_num,i);
    xFFT = 0:Fs/nfft:(Fs/2-Fs/nfft);
    yFFT = 2*abs(yFFT/L_signal);
%     plot(xFFT,yFFT,'-b');
%     hold on
%     xlim([0 5000]);
%     title(['FFT_Amp,T = ',num2str(T_arr(i)),'K']);
%     figure(f_p);
%     subplot(rows_num,cols_num,i);
%     plot(xFFT,pFFT,'-b');
%     hold on
%     xlim([0 5000]);
%     title(['FFT_Amp,T = ',num2str(T_arr(i)),'K']);
    
    dat_fft(i).x = xFFT;
    dat_fft(i).a = yFFT;
    dat_fft(i).p = pFFT;
end
%%
f_tmp = figure();
for i = 1:num
   plot(dat_fft(i).x,dat_fft(i).a,'-','color',cb(i).c);
   hold on
end
% xlim([0 5000]);
str_legend = 'legend(';
for i = 1:num
   str_legend = [str_legend,'''',num2str(T_arr(i)),''',']; 
end
str_legend(end) = [];
str_legend = [str_legend,');'];
eval(str_legend);
title(['FFT Amplitude at B_m = ',num2str(Bm),'T']);
xlabel('Frequency (T)');
ylabel('Amplitude (a.u.)');
xlim(range);

%% find peaks in fft data
% if(isfile([tmpDir,'T_dep_fft_',H_direction,'.m']) > 0)
%    load([tmpDir,'T_dep_fft_',H_direction,'.m']);
% else
%     disp('No fft data files found!');
% end

f_tmp2 = figure();
for i = 1:num
    tmp_index = dat_fft(i).x > (range(1) - i*drange) & dat_fft(i).x < (range(2) + i*drange);
   [peaks(i).pks,peaks(i).index] = findpeaks(dat_fft(i).a(tmp_index),'SortStr','descend','Npeaks',1);
   [peaks(i).pks,peaks(i).locs,peaks(i).w,peaks(i).pro] = findpeaks(dat_fft(i).a(tmp_index),dat_fft(i).x(tmp_index),'SortStr','descend','Npeaks',1);
   figure(f_tmp);
   hold on
   plot(peaks(i).locs,peaks(i).pks,'o','color',cb(i).c);
   figure(f_tmp2);
   tmp_p = dat_fft(i).p(tmp_index);
   peaks(i).pha = tmp_p(peaks(i).index);
   subplot(2,2,1);
   plot(peaks(i).locs,peaks(i).pks,'o','color',cb(i).c);
   hold on
   subplot(2,2,2);
   plot(peaks(i).locs,peaks(i).w,'o','color',cb(i).c);
   hold on
   subplot(2,2,3);
   plot(peaks(i).locs,peaks(i).locs,'o','color',cb(i).c);
   hold on
   subplot(2,2,4);
   plot(peaks(i).locs,peaks(i).pha,'o','color',cb(i).c);
   hold on
end

subplot(2,2,1);
title(['Peak values, B_m = ',num2str(Bm),'T']);
xlim(range);
subplot(2,2,2);
title('Peak Width');
xlim(range);
subplot(2,2,3);
title('Peak Location');
xlim(range);
subplot(2,2,4);
title('Peak phase');
xlim(range);

%% calculat effective mass of the first frequency
Bm = 1./((1/H_d + 1/H_u)/2);
pks_i = 1;
for i = 1:num
   if(length(peaks(i).pks)>0)
        A(i) = peaks(i).pks(pks_i);
        AoT(i) = peaks(i).pks(pks_i)/T_arr(i); 
        Loc(i) = peaks(i).locs(pks_i);
        Pha(i) = peaks(i).pha(pks_i);
   else
       A(i) = nan;
       AoT(i) = nan;
       Loc(i) = nan;
       Pha(i) = nan;
   end
end
tmp_index = isnan(A);
A(tmp_index) = [];
AoT(tmp_index) = [];
Loc(tmp_index) = [];
T_arr(tmp_index) = [];
num = num - sum(tmp_index);

figure();
subplot(1,2,1);
plot(T_arr,A,'or');
xlabel('T (K)');
ylabel('Amplitude');

subplot(1,2,2);
plot(T_arr,AoT,'or','MarkerFaceColor','r');
hold on
fpeak = fittype(['pre*(14.96*m/',num2str(Bm),')/sinh(14.69*m*T/',num2str(Bm),')'],'independent','T');
% fitp = fit(T_arr(T_arr<3.5)',AoT(T_arr<3.5)',fpeak,'start',[0.1,100]);
fitp = fit(T_arr',AoT',fpeak,'start',[1,10]);
errs=confint(fitp,confi_level);
pm.ms(pks_i) = fitp.m;
pm.ms_err(pks_i) = (errs(2,1)-errs(1,1))/2;
pm.loc(pks_i) = mean(Loc);
pm.loc_err(pks_i) = std(Loc);
pm.pha(pks_i) = mean(Pha);
pm.pha_err(pks_i) = std(Pha);
pm.T{pks_i} = T_arr;
% pm.mT{pks_i} = log(14.69*fitp.m*fitp.pre/Bm./AoT + sqrt((14.69*fitp.m*fitp.pre/Bm./AoT).^2+1))*Bm/14.96;
fitlinex=(1:0.1:4);
fitlinep=feval(fitp,fitlinex);
plot(fitlinex,fitlinep,'-k');
xlabel('T (K)');
ylabel('Amp/T');
disp(['Effective mass = ',num2str(pm.ms(pks_i)),' ',char(177),' ',num2str(pm.ms_err(pks_i))]);
disp(['Frequency is ',num2str(pm.loc(pks_i)),' ',char(177),' ',num2str(pm.loc_err(pks_i))]);
disp(['Phase is ',num2str(pm.pha(pks_i)),' ',char(177),' ',num2str(pm.pha_err(pks_i))]);