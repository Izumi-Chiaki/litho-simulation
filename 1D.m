clear,clc
close all

%% parameters

NA=1.1;       %NA
lambda=193;   %wavelength
sigma1=0.85;  %source radius
sigma2=0.6;   
Fs=10;        %sampling frequence

%%

T=1/Fs;    
x=[-500:T:500];  % x range,symmetric to the Y-axis
L=length(x);

%% inconherent

fx=x*Fs/(max(x)-min(x));
fx_p_source=fx(find(abs(fx)<=sigma1*NA/lambda & abs(fx)>=sigma2*NA/lambda));
p_source_num=length(fx_p_source);
I_m=zeros(1,length(fx));

%% mask imf

mask=zeros(1,length(x));
width=70;
for k=1:length(x)
    if ( abs(x(k))<=width/2 ) 
        mask(k)=0;
    else
        mask(k)=1;
    end
end
figure(1);
plot(x,mask,'linewidth',1);
xlabel('position/nm');
ylabel('Mask transmission');
title('Mask transmission')
T_m=fft(mask);
T_m=fftshift(T_m);
%fx=fft(x);
%figure(2)
%plot(fx,abs(T_m));

%% cycle

for j=1:p_source_num

    %% pupil imf

    for k=1:L
        if abs(fx(k)-fx_p_source(j))>NA/lambda
            T_m(k)=0;
        end
    end

    %% wafer

    %E_m=ifftshift(T_m);
    E_m=ifft(T_m);
    I_m=I_m+abs(E_m).^2;
end
I_m=I_m/p_source_num;
figure(3);
plot(x,I_m,'linewidth',1);
xlabel('position/nm');
ylabel('Intensity/au');
title('Aerial Image')