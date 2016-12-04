clear
L=0.01; imax=12; 
Ei=linspace(0,1,imax-1);
beta=1.2;
beta1=beta+1;
beta2=beta-1;
beta3=(beta1/beta2).^(2*Ei-1);
num=(beta1*beta3)-beta2;
den=2*(1+beta3);
x=L*num./den;
y=zeros(1,imax-1); 

xc(2:imax-1)=(x(2:imax-1)+x(1:imax-2))/2;
xc(1)=x(1);xc(imax)=x(imax-1);
yc=zeros(1,imax);
figure(1)
title('Non uniform Grid using Algebraic Method');
hold on
plot(x,y,'m-s')
plot(xc,yc,'k-o')