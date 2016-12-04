clear all;
clc; close all;
epsilon_st=0.001;
type=input('Enter gridtype = 1,2 or 3 for O,C or H gridtype respectively: ');
if type == 1
 imax=9; jmax=4;
 dzi=1/(imax-1); dzeta=1/(jmax-1);
 xold=ones(jmax,imax); yold=ones(jmax,imax);
 x=ones(jmax,imax); y=ones(jmax,imax); z=ones(jmax,imax);
 % % O-type grid generation Boundary Conditions % %
 for i=1:imax-1
 x(1,i)=5+2*cos(pi*(i-1)/4);% Inner Circle Boundary Condition on X
 y(1,i)=5+2*sin(2*pi-pi*(i-1)/4);% Inner Circle Boundary Condition on Y
 end
 % Outer Boundary Conditions
 x(jmax,1)=10; y(jmax,1)=5;
 x(jmax,2)=10; y(jmax,2)=0;
 x(jmax,3)=5; y(jmax,3)=0;
 x(jmax,4)=0; y(jmax,4)=0;
 x(jmax,5)=0; y(jmax,5)=5;
 x(jmax,6)=0; y(jmax,6)=10;
 x(jmax,7)=5; y(jmax,7)=10;
 x(jmax,8)=10; y(jmax,8)=10;
 x(jmax,9)=10; y(jmax,9)=5;
 % Connecting line Boundary Condition
 x(1,1)=7; x(2,1)=8; x(3,1)=9; x(4,1)=10; y(1:4,1)=5;
 x(1,imax)=7; x(2,imax)=8; x(3,imax)=9; x(4,imax)=10; y(1:4,imax)=5;
end

if type== 2
 imax=16; jmax=5;
 dzi=1/(imax-1); dzeta=1/(jmax-1);
 xold=rand(jmax,imax); yold=rand(jmax,imax);
 x=rand(jmax,imax); y=rand(jmax,imax); z=ones(jmax,imax);
 % % C-type grid generation Boundary Conditions % % %
 %Inner Boundary Conditions
 for i=4:13
   x(1,i)=5+2*cos(2*pi*(i-4)/9);
   y(1,i)=5+2*sin(2*pi-2*pi*(i-4)/9);
 end
 %Outer Boundary Condition
 x(1:5,1)=10;
 for j=1:5
   y(j,1)=5-10*(j-1)/8;
 end
 x(1:5,imax)=10;
 for j=1:5
   y(j,imax)=5+10*(j-1)/8;
 end
 y(jmax,1:6)=0;
 for i=1:6
    x(jmax,i)=10-10*(i-1)/5;
 end
 x(jmax,6:11)=0;
 for i=6:11
    y(jmax,i)=10*(i-6)/5;
 end
 y(jmax,11:16)=10;
 for i=11:16
    x(jmax,i)=10*(i-11)/5;
 end
 %Connecting line Boundary Condition
 y(1,1:4)=5;
 for i=1:4
    x(1,i)=11-i;
 end
 y(1,13:16)=5;
 for i=13:16
    x(1,i)=7+(i-13);
 end
end

if type==3;
 imax=11; jmax=5;
 dzi=1/(imax-1); dzeta=1/(jmax-1);
 xold=rand(jmax,imax); yold=rand(jmax,imax);
 x=rand(jmax,imax); y=rand(jmax,imax); z=ones(jmax,imax); 
 % % H-type grid generation Boundary Conditions % % %
 %Inner Boundary Conditions
 for i=4:8
   x(1,i)=5+2*cos(pi*(8-i)/4); 
   y(1,i)=2*sin(pi*(8-i)/4); 
 end
 % outer Boundary Condition
 x(1:jmax,1)=0;
 for j=1:jmax
    y(j,1)=5*(j-1)/4;
 end
 y(jmax,1:imax)=5;
 for i=1:imax
    x(jmax,i)=i-1;
 end
 x(1:jmax,imax)=10;
 for j=1:jmax
    y(j,imax)=5*(j-1)/4;
 end 
 % Connecting line Boundary Condition
 y(1,1:4)=0;
 for i=1:4
     x(1,i)=i-1;
 end
 y(1,8:11)=0;
 for i=8:11
    x(1,i)=i-1;
 end
end

Error=1;
while Error>=epsilon_st
   xold=x;
   yold=y;
   for i=2:imax-1
       for j=2:jmax-1
           A(j,i)=((xold(j+1,i)-xold(j-1,i))/(2*dzeta))^2+((yold(j+1,i)-yold(j-1,i))/(2*dzeta))^2;
           B(j,i)=((xold(j,i+1)-xold(j,i-1))*(xold(j+1,i)-xold(j-1,i))+(yold(j,i+1)-yold(j,i-1))*(yold(j+1,i)-yold(j-1,i)))/(4*dzi*dzeta);
           C(j,i)=((xold(j,i+1)-xold(j,i-1))/(2*dzi))^2+((yold(j,i+1)-yold(j,i-1))/(2*dzi))^2;
           Dx(j,i)=-2*B(j,i)*(xold(j+1,i+1)+xold(j-1,i-1)-xold(j+1,i-1)-xold(j-1,i+1))/(4*dzi*dzeta);
           Dy(j,i)=-2*B(j,i)*(yold(j+1,i+1)+yold(j-1,i-1)-yold(j+1,i-1)-yold(j-1,i+1))/(4*dzi*dzeta);
           x(j,i)=(A(j,i)*(x(j,i+1)+x(j,i-1))*(dzeta*dzeta)+C(j,i)*(x(j+1,i)+x(j-1,i))*(dzi*dzi)+Dx(j,i)*dzi*dzi*dzeta*dzeta)/(2*A(j,i)*(dzeta*dzeta)+2*C(j,i)*(dzi*dzi));
           y(j,i)=(A(j,i)*(y(j,i+1)+y(j,i-1))*(dzeta*dzeta)+C(j,i)*(y(j+1,i)+y(j-1,i))*(dzi*dzi)+Dy(j,i)*dzi*dzi*dzeta*dzeta)/(2*A(j,i)*(dzeta*dzeta)+2*C(j,i)*(dzi*dzi));
       end
   end
   Error=max(max(abs(x(j,i)-xold(j,i)),abs(y(j,i)-yold(j,i))));
   disp(Error)
end

figure;
if (type==3)
    for k=2:jmax-1
      plot (x(k,:),y(k,:)+5,'r-o','linewidth',1.5,'markeredgecolor','k','markersize',9,'markerfacecolor','y'); hold on;
    end
    for k=2:imax-1
     plot (x(:,k),y(:,k)+5,'r-o','linewidth',1.5,'markeredgecolor','k','markersize',9,'markerfacecolor','y'); hold on;
    end
    plot (x(1,:),((-1)*y(1,:))+5,'r-','linewidth',1.5,'markeredgecolor','k','markersize',9,'markerfacecolor','y'); hold on;
    for k=2:jmax-1
     plot (x(k,:),((-1)*y(k,:))+5,'r-','linewidth',1.5,'markeredgecolor','k','markersize',9,'markerfacecolor','y'); hold on;
    end
    for k=2:imax-1
     plot (x(:,k),((-1)*y(:,k))+5,'r-','linewidth',1.5,'markeredgecolor','k','markersize',9,'markerfacecolor','y'); hold on;
    end
    plot (x(1,:),y(1,:)+5,'k-o','linewidth',2,'markeredgecolor','k','markersize',9,'markerfacecolor','b'); hold on;
    plot (x(jmax,:),y(jmax,:)+5,'k-o','linewidth',2,'markeredgecolor','k','markersize',9,'markerfacecolor','b'); hold on;
    plot (x(:,1),y(:,1)+5,'k-o','linewidth',2,'markeredgecolor','k','markersize',9,'markerfacecolor','b'); hold on;
    plot (x(:,imax),y(:,imax)+5,'k-o','linewidth',2,'markeredgecolor','k','markersize',9,'markerfacecolor','b'); hold on;
    hold off;
else
    for k=2:jmax-1
     plot (x(k,:),y(k,:),'r-o','linewidth',1.5,'markeredgecolor','k','markersize',9,'markerfacecolor','y'); hold on;
    end
    for k=2:imax-1
     plot (x(:,k),y(:,k),'r-o','linewidth',1.5,'markeredgecolor','k','markersize',9,'markerfacecolor','y'); hold on;
    end
    plot (x(1,:),y(1,:),'k-o','linewidth',2,'markeredgecolor','k','markersize',9,'markerfacecolor','b'); hold on;
    plot (x(jmax,:),y(jmax,:),'k-o','linewidth',2,'markeredgecolor','k','markersize',9,'markerfacecolor','b'); hold on;
    plot (x(:,1),y(:,1),'k-o','linewidth',2,'markeredgecolor','k','markersize',9,'markerfacecolor','b'); hold on;
    plot (x(:,imax),y(:,imax),'k-o','linewidth',2,'markeredgecolor','k','markersize',9,'markerfacecolor','b'); hold on;
    hold off; 
end