function [FFT] = FFTcalc(obj,dataType,FFrange,FFdir,varargin)
   % FFTcalc Fourier transforms dHvA data.
   % 
   % FFT = FFTcalc(raw, dataType, FFrange, FFdir, fRange)the raw data
   % (raw) is loaded into FFTcalc as an struct: raw.x and raw.y. The data
   % falling within the specified range (FFrange) is extracted,
   % and fourier transforms the data after subtracting a polinomial fit
   % (the degre varies with the type of data specified by variable dataType, 
   % and may require the user to manualy change the degree in this
   % function), adding zero pading, and a Hanning window. FFdir specifies
   % wether the data is taking from the rising pulse ('up') or falling
   % pulse ('down'). 
%     obj
    varargin
    fRange = varargin{1};
    FFdirNum = strcmp(FFdir,'up');
    Tbot = FFrange(1)-1;
    Ttop = FFrange(2)+1;
    chTbot = FFrange(1);
    chTtop = FFrange(2);

    switch FFdir
        case 'up'   
            FF = obj.xUp;%Down(end:-1:1);
            xi = obj.yUp;%Down(end:-1:1);
%             FF(high:end) = zeros(length(FF(high:end)),1);  
            FFInd = FF >= Tbot & FF <= Ttop;
            xinv = FF(FFInd);
%             int = .0000001;%change to the power of 8
            x = 1./xinv;%
            N = 2^20;%round(abs(x(1)-x(end))/int); %Total number of points of array fed into FFT 
            
            if rem(N,2) == 1
                N = N+1;
            end
            
        case 'down'
            FF = obj.xDown;
            xi = obj.yDown;
            FF(1:high) = zeros(high,1);
            FFInd = FF >= Tbot & FF <= Ttop;
            xinv = FF(FFInd);
            int = .0000001;
            x = 1./xinv;
            if rem(N,2) == 1
                N = N+1;
            end
%             chN = N
            xspl = linspace(x(end),x(1),N);
    end

    n = 2^15 ; %round(N*2/3); %number of points fit to the spline %;
    xspl = linspace(1/chTtop,1/chTbot,n); %
%     chx = []
%     chorig = [x(end),x(1)]
%     chch = [1/chTtop,1/chTbot]
    int = xspl(2)-xspl(1);
    y = xi(FFInd);
    yog = y;
    
    if strcmp(dataType,'dHvA')
        degree =2;
    elseif strcmp(dataType,'rfPD')
        degree = 2;
    elseif strcmp(dataType,'whatever')
        degree = 0;
    end
    yfit = polyfit(xinv,y,degree);
    yfitV = polyval(yfit,xinv);
    
%     figure
%     subplot(2,1,1)
%     plot(xinv,y,'r',xinv,yfitV,'b')
    
    y = y-yfitV;
    [x,Iuniq] = unique(x);
    y = y(Iuniq);
    
    % spline interpolation
    yspl = interp1(x,y,xspl);%CHANGED !!! from spline
    
    Length=abs(xspl(end)-xspl(1));
    Datapoints = length(yspl);
    fs=Datapoints/Length;
%     var3fftc = varargin

    figure
    plot(xspl,yspl)
    
    if length(varargin) == 3
%         varargin3 = varargin{3}
        switch varargin{2}
            case 'lowpass' 
                yspl = lowpass(yspl,varargin{3},fs);%10000
            case 'bandpass'
                yspl = bandpass(yspl,varargin{3},fs);%[3500 3900]
        end
    end

%     subplot(2,1,2)
%     plot(1./xspl,yspl,'r')%,1./x,y,'.b')

%     chyspl1 = yspl(1)
%     chysplend = yspl(end)
% 

     
%     figure
%     plot(1./xspl,yspl1,'r')
%     pause(.5)
%     hold on55
%     plot(1./xspl,yspl,'b')
%     temp = strcat(num2str(obj.temp),'K');
%     title(temp)
    
    %Hanning window
    Datapoints = length(yspl);
    hwind = hann(Datapoints)';
    ysplW = yspl.*hwind;
    
    %Zeri pading
    nZero = (N-n)/2;
    xfft1 = xspl(1)-floor(nZero)*int;
    xfftEnd = xspl(end)+ceil(nZero)*int;   
    zeros1 = zeros(1,floor(nZero));
    zerosEnd = zeros(1,ceil(nZero));
        
    xfft = linspace(xfft1,xfftEnd,N);%xspl;%
    yfft = [zeros1,ysplW,zerosEnd];%ysplW;%yspl;%
         
    figure
    plot(x,y,'b',xfft,yfft,'r','LineWidth',1.2)
    
    on = 0;
    if on == 1
        figure
        subplot(2,1,1)
        plot(xinv,yog,'b',xinv,yfitV,'r')
        subplot(2,1,2)
        plot(xinv,y)
        temp = strcat(num2str(obj.temp),'K');
        suptitle(strcat( 'T=', temp))
    end
    
    %Calculate fourier transform and define FFT struct 
    Datapoints = length(yfft);
    FFTc = fft(yfft);%./length(yfft);%fft(yv);%   
    disp(FFTc)
    FFTval = abs(FFTc(1:Datapoints/2+1));
    
%     figure
%     plot(yfft)
    
%     tollerence = max(FFTval)/100
%     FFTc(abs(FFTc)<tollerence) = 0;
    FFTphase = angle(FFTc);
    FFTphaseVal = FFTphase(1:Datapoints/2+1);
    FFT.x = x;
    FFT.y = y;
    FFT.xspl = xspl;
    FFT.yspl = yspl;%yv;%
    
    Length=abs(xfft(end)-xfft(1)); 
    fs=Datapoints/Length; 
    f = fs/2*linspace(0,1,Datapoints/2+1);
    fI = f>=fRange(1) & f<=fRange(2);
    
    FFT.f = f(fI);
    FFT.FFT = FFTval(fI);
    FFT.phase = FFTphaseVal(fI);
    FFT.temp = obj.temp;
    FFT.range = FFrange;

end