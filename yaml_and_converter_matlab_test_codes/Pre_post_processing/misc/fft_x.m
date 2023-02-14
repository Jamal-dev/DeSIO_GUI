function [ft_peak,f,P] = fft_x(fig,u,t,range)
    L  = length(u); % Length of signal
    T  = t(2)-t(1);      % Sampling period
    Fs = 1/T;            % Sampling frequency
    Y  = fft(u);
    
    P2 = abs(Y/L);
    P1 = P2(1:fix(L/2)+1);
    P1(2:end-1) = 2*P1(2:end-1);
    f = Fs*(0:(L/2))/L *2*pi;
%     f = f;
    [val,inzm] = max(P1);
%     P1 = P1/val;
    figure(fig); hold on; grid on; %axis([0 5 0 1.1*max(P1)]);
    title('Single-Sided Amplitude Spectrum');
    xlabel('\omega rad/s'); ylabel('|fft(u)|');
    rgb_color = RGB_Color;
%     plot(f,P1,'-','color',rgb_color(i,:),'linewidth',1.8); 
    [c,inzf,val] = find(f>=range(1) & f<=range(2));
    plot(f(inzf),P1(inzf),'-','linewidth',1.0); 
%     text(f(inzm),P1(inzm),[num2str(f(inzm)) ' Hz']);
    [val,inz] = findpeaks(P1(inzf));
    f = f(inzf);
    P = P1(inzf);
    ft_peak = f(inz);
return