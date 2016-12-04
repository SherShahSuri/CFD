clear
L1=1;L2=1;   imax=12;jmax=12;Ei=linspace(0,1,imax-1);Fi=linspace(0,1,jmax-1);
Beta=1.2;
Beta_p1=Beta+1;Beta_m1=Beta-1;Beta_p1_div_m1=(Beta_p1/Beta_m1).^(2*Ei-1);
num=(Beta_p1*Beta_p1_div_m1)-Beta_m1;den=2*(1+Beta_p1_div_m1);
x=L1*num./den;
y=L2*num./den;

xc(2:imax-1)=(x(2:imax-1)+x(1:imax-2))/2;
xc(1)=x(1);xc(imax)=x(imax-1);

yc(2:imax-1)=(y(2:imax-1)+y(1:imax-2))/2;
yc(1)=y(1);yc(imax)=y(imax-1);

figure(1)
for i=1:imax-1
    for j=1:imax-1
plot(x(i),y(j),'m-s')
hold on
    end
end
for i=1:imax
    for j=1:jmax
        plot(xc(i),yc(j),'k-o')
    end
end
