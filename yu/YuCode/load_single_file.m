function res = load_single_file(filename,flag_dec,field_bkg_range,field_range,if_plot)
    % if_plot == 1; make plots
    % flag_dec == 1; choose field down scan
    % field_range = [min,max];
    % frequency  = frequency/ 1e5
    fid = fopen(filename);
    H = []; Fre = []; Amp = [];
    while 1
        S = fgetl(fid);
        if (S==-1) break; end
        s = textscan(S,'%s','Delimiter','\t');
        s = s{:};
        H = [H; str2num(s{1})];
        Fre = [Fre; str2num(s{2})];
        Amp = [Amp; str2num(s{3})];
    end
    fclose(fid);

    index = find(H == max(H),1);
    if(flag_dec == 1)
        % decreasing field
        H_tmp = H(index+1:end);
        Fre_tmp = Fre(index+1:end);
        Amp_tmp = Amp(index+1:end);
    elseif(flag_dec == 2)
        H_tmp = H(1:index);
        Fre_tmp = Fre(1:index);
        Amp_tmp = Amp(1:index); 
    else
        H_tmp_d = H(index+1:end);
        Fre_tmp_d = Fre(index+1:end);
        Amp_tmp_d = Amp(index+1:end);
        H_tmp_u = H(1:index);
        Fre_tmp_u = Fre(1:index);
        Amp_tmp_u = Amp(1:index); 
    end

    % fit polynomial
    if(flag_dec <= 2)
        [H,tmp] = sort(H_tmp,'ascend');
        Fre = Fre_tmp(tmp);
        Amp = Amp_tmp(tmp);
        
        fit_range = H > field_bkg_range(1) & H < field_bkg_range(2);
    
        dat_H = H(fit_range);
        dat_F = Fre(fit_range);
        dat_F = dat_F/1e5;

        f_res = polyfit(dat_H,dat_F,3);
        bkg_H = dat_H;
        bkg_F = polyval(f_res,bkg_H);

        index = dat_H > field_range(1) & dat_H < field_range(2);
        res.H = dat_H(index);
        res.invH = 1./dat_H(index);
        res.F = (dat_F(index) - bkg_F(index))*1e5;
    else    %% combine both inc & dec
        % field down
        fit_range = H_tmp_d > field_range(1) & H_tmp_d < field_range(2);
    
        dat_H = H_tmp_d(fit_range);
        dat_F = Fre_tmp_d(fit_range);
        dat_F = dat_F/1e5;

        f_res = polyfit(dat_H,dat_F,3);
        bkg_H = dat_H;
        bkg_F = polyval(f_res,bkg_H);

        tmpres.H = dat_H;
        tmpres.invH = 1./dat_H;
        tmpres.F = (dat_F - bkg_F)*1e5;
        %  field up
        fit_range = H_tmp_u > field_range(1) & H_tmp_u < field_range(2);
    
        dat_H = H_tmp_u(fit_range);
        dat_F = Fre_tmp_u(fit_range);
        dat_F = dat_F/1e5;

        f_res = polyfit(dat_H,dat_F,3);
        bkg_H = dat_H;
        bkg_F = polyval(f_res,bkg_H);
        
        tmpres.H = [tmpres.H;dat_H];
        tmpres.invH = [tmpres.invH;1./dat_H];
        tmpres.F = [tmpres.F;(dat_F - bkg_F)*1e5];
        
        index = tmpres.H > field_range(1) & tmpres.H < field_range(2);
        tmpres.H = tmpres.H(index);
        tmpres.invH = tmpres.invH(index);
        tmpres.F = tmpres.F(index);
        
        [res.H,tmp_index] = unique(tmpres.H);
        res.invH = tmpres.invH(tmp_index);
        res.F = tmpres.F(tmp_index);       
    end
    
    if(if_plot == 1)
       figure();
       subplot(1,2,1);
       plot(dat_H,dat_F,'-r');
       hold on
       plot(bkg_H,bkg_F,'-k');
       subplot(1,2,2);
       plot(res.H,res.F,'-r');
    end

end