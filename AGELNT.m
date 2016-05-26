
close all
clear all 
%%
% Extraction of s-parameters using Nicolson-Ross weir
% zakari i. emurana
% MSC.Thesis

Y=xlsread('DEPAULA_MAGNITYDE&PHASE.xlsx');
f_data=Y(1:1001,1);


S11phase=(Y(1:1001,4));
S21phase=(Y(1:1001,5));


 S11mag=(Y(1:1001,2));
 S21mag=(Y(1:1001,3));


% determine complex parameters

S11 = S11mag .* exp(j * (3.141 * S11phase/180));
S21 = S21mag .* exp(j * (3.141 * S21phase/180));


L=11.45*10^-3;  % sample length
w=45.72*10^3; % cutoff wavelenght
f=6.5571*10^9; % cutoff frequency

C = 2.997925*10^8; % speed of light


C = 2.997925*10^8; % speed of light

X= (S11.^2-S21.^2+1)./(2.*S11);

RFC_posit = X + sqrt (X.^2 - 1);
RFC_neget = X - sqrt (X.^2 - 1);

TC= (S11+S21-RFC_neget)./1-(S11+S21).*RFC_neget;
 

p=log(1./TC);
p1=p*pi./180;


V= -1./((1./2.*pi*L).*p1);
 







lamda_not= 2*3.142*f_data./3*10^8; % free spcae wavelength

% permeability 

 UR= 1+RFC_neget./V.*(1-RFC_neget).*sqrt((1./(lamda_not).^2)-(1./(w).^2));
 
UR_1=real(UR);
plot(f_data,UR_1)
hold on
plot(f_data,imag(UR))


ER= ((lamda_not.^2./UR).*((1./w^2)-(1./2.*pi*L.*p1).^2));

figure
 plot(f_data,real(ER))
hold on
plot(f_data,imag(ER))



