clear,clc
close all

%% parameters

NA=0.57;       %NA
lambda=365;   %wavelength
sigma1=0.7;  %source radius
sigma2=0;  
Fs=1;        %sampling frequence of x
% alpha=0.0005;
dill_a=0.00075;
dill_b=0.00005;
dill_c=0.0025; %the photoresist Dill's coefficient
lamp_power=30000; %lamp_power
dose=2000;        %dose
t_tot=dose/lamp_power;
nt=50;            %sampling number of t
nz=50;            %sampling number of z
thk=1000;         %photoresist thickness
r_max=50;         %maximum development rate
r_min=0.8;        %minimum development rate
ord_rea=2;        %reaction order
m_th=0.01;        %threshold concentration

%%

T=1/Fs;    
x=[-1000:T:1000];  % x range,symmetric to the Y-axis
L=length(x);

%% inconherent

fx=x*Fs/(max(x)-min(x));
fx_p_source=fx(find(abs(fx)<=sigma1*NA/lambda & abs(fx)>=sigma2*NA/lambda));
p_source_num=length(fx_p_source);
I_m=zeros(1,length(fx));

%% mask imf

mask=ones(1,length(x));
%width=70;
%pitch=140;
%mask=(square(2*pi/pitch*x)+1)/2;
%x=x-pitch/4;
%mask(find(x(x<-35)))=1;
for j=1:length(x)
    if x(j)<=-125|x(j)>=125
        mask(j)=0;
    end
end
figure(1);
plot(x,mask,'linewidth',1);
xlabel('position/nm');
ylabel('Mask transmission');
title('Mask transmission');
%xlim([-1000,1000]);
T_m=fft(mask);
T_m=fftshift(T_m);
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
title('Aerial Image');
%xlim([-2000,3000]);
%hold on

%% exposure
dz=thk/nz;
z=linspace(0,thk,dz);
dt=t_tot/nt;
m_xz=ones(length(z),length(x));
[x_grid,z_grid]=meshgrid(x,z);
[I_m_grid,z_grid]=meshgrid(I_m,z);
I_xz=ones(length(z),length(x));
I_xz=I_m_grid.*exp(-(dill_a*m_xz+dill_b).*z_grid);
% for i=1:length(z)
%     I_xz(i,:)=I_m*exp(-alpha*z(i));
% end
% figure(4);
% colormap(othercolor('Pastel13'));
% contourf(x_grid,z_grid,I_xz);
% set(gca,'YDir','reverse'); 
% xlabel('x/nm');
% ylabel('resist thickness/nm');
% colorbar
% title('Original exposition');
for t_e=dt:dt:t_tot
    m_xz=m_xz.*exp(-dill_c*I_xz*dt*lamp_power);
    I_xz=I_m_grid.*exp(-(dill_a*m_xz+dill_b).*z_grid);
end
figure(5);
colormap(othercolor('Pastel13'));
contourf(x_grid,z_grid,I_xz);
set(gca,'YDir','reverse'); 
xlabel('x/nm');
ylabel('resist thickness/nm');
colorbar
title('Exposition');

%% Development
a_R_xz=(ord_rea+1)/(ord_rea-1)*(1-m_th)^ord_rea;
R_xz=r_max*((a_R_xz+1)*(1-m_xz).^ord_rea)./(a_R_xz+(1-m_xz).^ord_rea)+r_min;
t_xz=cumtrapz(1./R_xz,1);
figure(6);
colormap(othercolor('Pastel13'));
contourf(x_grid,z_grid,t_xz);
set(gca,'YDir','reverse'); 
xlabel('x/nm');
ylabel('resist thickness/nm');
colorbar
title('Needed development time');