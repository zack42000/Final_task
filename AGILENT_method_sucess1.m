%%
% Extraction of s-parameters using Nicolson-Ross weir
% zakari i. emurana
% MSC.Thesis

Y1=xlsread('real&image_depual.xlsx');
f_data=Y1(1:1001,1);

S11imag=Y1(1:1001,2);
S21imag=Y1(1:1001,3);
S11real=Y1(1:1001,4);
S21real=Y1(1:1001,5);

% COMPLEX
S11=complex(S11real,S11imag);
S21=complex(S21real,S21imag);


L=11.45*10^-3;  % sample length
%w=45.72*10^-3; % cutoff wavelenght
w=pi/22.86;
f=6.5571*10^9; % cutoff frequency

C = 2.997925*10^8; % speed of light

X= (S11.^2-S21.^2+1)./(2.*S11);

RFC_posit = X + sqrt (X.^2 - 1);
RFC_neget = X - sqrt (X.^2 - 1);

rfc_n=abs((RFC_neget));
rfc_p=abs(RFC_posit);

RFC=( abs ( RFC_posit ) <=1) .* RFC_posit + ( abs ( RFC_neget ) <=1) .* RFC_neget ;


TC= (S11+S21-RFC)./1-(S11+S21).*RFC;

p=log(1./TC);

%p1=p*pi./180;
V= -1./((1./2.*pi*L).*p);

lamda_not= 2*3.142*f_data./3*10^8; % free spcae wavelength

 
 UR= 1+RFC./V.*(1-RFC).*sqrt((1./(lamda_not).^2)-(1./(w).^2));
 
 
 UR_1=real(UR);
 figure
plot(f_data,UR_1)
hold on
plot(f_data,imag(UR))
xlabel('frequency')
ylabel('U* and U*')




ER= ((lamda_not.^2./UR).*((1./w.^2)-(1./2.*pi*L.*p).^2));

figure
 plot(f_data,real(ER))
hold on
plot(f_data,imag(ER))
xlabel('frequency')
ylabel('E* and E*')


