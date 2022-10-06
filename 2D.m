clear,clc
close all

%% parameters

NA=1.1;       %NA
lambda=193;   %wavelength
sigma1=0.85;  %source radius
sigma2=0.6;   
Fs_x=5;
Fs_y=0.1;       %sampling frequence

%%

T_x=1/Fs_x;
T_y=1/Fs_y;   
x=[-500:T_x:500];  % x range,symmetric to the Y-axis
y=[-1000:T_y:1000];

%% inconherent

fx=x*Fs_x/(max(x)-min(x));
fy=y*Fs_y/(max(y)-min(y));
f_2D=zeros(length(y),length(x));
p_source_num=0;
for j=1:length(x)
    for k=1:length(y)
        f_2D(k,j)=(fy(k)^2+fx(j)^2).^0.5;
        if f_2D(k,j)<=sigma1*NA/lambda & f_2D(k,j)>=sigma2*NA/lambda
            p_source_num=p_source_num+1;
            f_p_source(p_source_num)=f_2D(k,j);
        end
    end
end
I_m=zeros(length(fy),length(fx));

%% mask imf

mask=zeros(length(y),length(x));
width=70;
length_line=1000;
for j=1:length(x)
    for k=1:length(y)
        if ( abs(x(j))<=width/2 | ( x(j)>=105 & x(j)<=175 ) | ( x(j)>=-175 & x(j)<=-105 )) & ( abs(y(k))<=length_line/2 )
            mask(k,j)=0;
        else
            mask(k,j)=1;
        end
    end
end
figure(1);
meshgrid(x,y);
colormap(othercolor('Pastel13'));
contourf(x,y,mask);
xlabel('x');
ylabel('y');
colorbar
title('Mask transmission');
T_m=fft2(mask);
T_m=fftshift(T_m);
%fx=fft(x);
%figure(2)
%plot(fx,abs(T_m));

%% cycle

for j=1:p_source_num

    %% pupil imf

    for k=1:length(x)
        for t=1:length(y)
            if abs(f_2D(t,k)-f_p_source(j))>NA/lambda
                T_m(t,k)=0;
            end
        end
    end

    %% wafer

    %E_m=ifftshift(T_m);
    E_m=ifft2(T_m);
    I_m=I_m+abs(E_m).^2;
end
I_m=I_m/p_source_num;
figure(3);
meshgrid(x,y);
colormap(othercolor('Pastel13'));
contourf(x,y,I_m);
xlabel('x')
ylabel('y')
colorbar
title('Aerial Image')