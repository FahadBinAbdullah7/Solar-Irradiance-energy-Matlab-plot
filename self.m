
clc
clear all
close all

l=23.750; % LATTITUDE
n=219; %DAY 219, AUG 7


a= 23.45*(sind((n+284)*(360/365))) %declination angle

tlt = l ;    % TILT ANGLE                

AA = 1160+75*sind((360/365)*(n-275));
kk = 0.174+0.035*sind((360/365)*(n-100));
cc = 0.095+(0.04*sind((360/365)*(n-100)));
alb = 0.2; %GROUND reflectance


Ws = acosd((-tand(l)*tand(a)))  %sunrise angle

% local time not equal to solar time. 
%12 = SOLAR NOON TIME ; local time = 11;12;13


Sr=12-((l/219)*(acosd(-tand(l)*tand(a))))%sunrise time
Ss=12+((l/219)*(acosd(-tand(l)*tand(a))))% sunset time

intsr=floor(Sr)+1;
intss=floor(Ss);    %% JUST FOR MATLAB TIME INTERVAL SMOOTHNESS


T=Ss-Sr;  %total day time

timle=[Sr,intsr:0.25:intss,Ss]                                                   % 60/15 =4 ; 1/4 = ?

p=length(timle);

for i=1:(p)
    
    ws=(-Ws+(((2*Ws)/T)*(timle(i)-Sr))) %hour angle

  TTIME = timle(i)  
    

A=asind((sind(a)*sind(l))+(cosd(a)*cosd(l)*cosd(ws))) %solar altitude angle

Za=90-A;%zenith angle

AM=(1/cosd(Za)) %air mass

AM2=(1/sind(A)) %air mass

 fys=asind((cosd(a)*sind(ws))/cosd(A)); %SOLAR AZIMUTH
 
 % fi2=asind(sind(ws)/(sind(bft))) ;  %if surface azimuthal exist
% kosh= (cosd(A)* cosd(fys-fi2)*sind(bft))+(sind(A)*cosd(bft))

 kosh= (cosd(A)* cosd(fys-0)*sind(tlt))+(sind(A)*cosd(tlt)) %incident angle



Ib = AA*(exp(-kk*AM)); %BEAM IRRADIANCE


if(Ib==inf)
    Ib=0;
else
    Ib=AA*(exp(-kk*AM));
end

Ib

% Io=1367*((0.7)^(AM^(0.678)));

Idrect(i) = Ib*kosh;

refactpf = ((1-cosd(tlt))/2);
difactmf = ((1+cosd(tlt))/2); %sky


 idt(i) = cc*Ib*difactmf; %DIFFUSE IRRADIANCE
 
 irt(i) = alb*Ib*(sind(A)+cc)*refactpf;   %REFLECTANCE IRRADIANCE

Ibeam = Idrect(i)% BEAM IRRADIANCE STRIKE TO THE COLLECTOR
Idiffuse = idt(i)
Ireflect = irt(i)


total(i)= irt(i)+idt(i)+Idrect(i);


end

total_irradiance=sum(total)
total_direct = sum(Idrect)
total_diffusion = sum(idt)
total_reflect = sum(irt)

direct_portion = (total_direct/total_irradiance)*100
diffusion_portion = (total_diffusion/total_irradiance)*100
reflect_portion = (total_reflect/total_irradiance)*100


Energy_S=total_irradiance*0.25
monotb=total;

figure(12)
plot (timle,total,'k',timle,Idrect,'-.',timle,idt,'-',timle,irt,':','LineWidth',2),legend('Total IRRADIANCE', 'DIRECT IRRADIANCE','DIFFUSE IRRADIANCE','REFLECTED IRRADIANCE');
xlabel('Time(Hours)')
ylabel('Irradiance(Watt/m^2)')
str=sprintf('Solar Irrediance');
title(str);