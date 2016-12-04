clear
% Grid Generation.
L=6; H=1;
imax=62; jmax=22;
xi=linspace(0,1,imax-1);
Beta=1.2;
Beta_p1=Beta+1;Beta_m1=Beta-1;Beta_p1_div_m1=(Beta_p1/Beta_m1).^(2*xi-1);
num=(Beta_p1*Beta_p1_div_m1)-Beta_m1;den=2*(1+Beta_p1_div_m1);
x=L*num./den;

yi=linspace(0,1,jmax-1);
Beta=1.2;
Beta_p1=Beta+1;Beta_m1=Beta-1;Beta_p1_div_m1=(Beta_p1/Beta_m1).^(2*yi-1);
num=(Beta_p1*Beta_p1_div_m1)-Beta_m1;den=2*(1+Beta_p1_div_m1);
y=H*num./den;

xc(2:imax-1)=(x(2:imax-1)+x(1:imax-2))/2;
xc(1)=x(1);xc(imax)=x(imax-1);

yc(2:jmax-1)=(y(2:jmax-1)+y(1:jmax-2))/2;
yc(1)=y(1);yc(jmax)=y(jmax-1);

% Ploting the grid points.
for i=1:imax-1
    plot(x(i),yc(2:jmax-1),'m-s'); hold on;
end
for j=1:jmax-1
    plot(xc(2:imax-1),y(j),'m-s'); hold on;
end
for i=1:imax
    plot(xc(i),yc,'g-o'); hold on;
end
hold off;

% Geometric Properties
Dx(2:imax-1)=x(2:imax-1)-x(1:imax-2);
dx(1:imax-1)=xc(2:imax)-xc(1:imax-1);
Dy(2:jmax-1)=y(2:jmax-1)-y(1:jmax-2);
dy(1:jmax-1)=yc(2:jmax)-yc(1:jmax-1);

 for i=2:imax-2           % Weights of Advection scheme
         
         if (i==2)
             w1q_x_plus(i)=((Dx(i))*((2*Dx(i))))/((Dx(i)+Dx(i+1))*(Dx(i+1)+2*Dx(i)));
             w2q_x_plus(i)=((Dx(i))*((2*Dx(i))))/((Dx(i)+Dx(i+1))*(Dx(i))); 
             w3q_x_plus(i)=((Dx(i))*((Dx(i+1))))/((Dx(i))*(Dx(i)+2*Dx(i+1)));
         else 
              w1q_x_plus(i)=((Dx(i))*((2*Dx(i))+Dx(i-1)))/((Dx(i)+Dx(i+1))*(Dx(i+1)+Dx(i-1)+2*Dx(i))); 
              w3q_x_plus(i)=((Dx(i))*((Dx(i+1))))/((Dx(i)+Dx(i-1))*(Dx(i+1)+Dx(i-1)+2*Dx(i)));
              w2q_x_plus(i)=((Dx(i))*((2*Dx(i))+Dx(i-1)))/((Dx(i)+Dx(i+1))*(Dx(i-1)+Dx(i)));
         end
     end
     
     for j=2:jmax-2 
        
         if (j==2)
                   w1q_y_plus(j)=((Dy(j))*((2*Dy(j))))/((Dy(j)+Dy(j+1))*(Dy(j+1)+2*Dy(j)));
                   w2q_y_plus(j)=((Dy(j))*((2*Dy(j))))/((Dy(j)+Dy(j+1))*(Dy(j)));
                   w3q_y_plus(j)=((Dy(j))*((Dy(j+1))))/((Dy(j))*(Dy(j)+2*Dy(j+1)));
             else 
                   w1q_y_plus(j)=((Dy(j))*((2*Dy(j))+Dy(j-1)))/((Dy(j)+Dy(j+1))*(Dy(j+1)+Dy(j-1)+2*Dy(j)));
                   w2q_y_plus(j)=((Dy(j))*((2*Dy(j))+Dy(j-1)))/((Dy(j)+Dy(j+1))*(Dy(j-1)+Dy(j)));
                   w3q_y_plus(j)=((Dy(j))*((Dy(j+1))))/((Dy(j)+Dy(j-1))*(Dy(j+1)+Dy(j-1)+2*Dy(j)));
         end
     end


% User Input.
rho = 10.0; cp = 1.0; k = 1; mu=1;
u=1; v=0; alpha=k/(rho*cp); Re=(rho*H*u)/mu; Pr=(mu*cp)/k;
epsilon_st=0.000001; epsilon=0.0001; Dt=1;
mx=rho*u; my=rho*v;


% Coefficient of implicit LAEs
 for i=2:imax-1
    for j=2:jmax-1 
     aE(j,i)=(((k)/dx(i)))*Dy(j);
     aW(j,i)=(((mx)* cp)+((k)/dx(i-1)))*Dy(j);
     aN(j,i)=((((k)/dy(j)))*Dx(i));
     aS(j,i)=((((my)* cp)+((k)/dy(j-1)))*Dx(i));
    end
 end
 
 for i=2:imax-1
    for j=2:jmax-1 
     aP0(j,i)=(rho*cp*Dx(i)*Dy(j))/Dt;
    end
 end
 
 for i=2:imax-1
    for j=2:jmax-1 
     aP(j,i)=aE(j,i)+aW(j,i)+aN(j,i)+aS(j,i)+aP0(j,i);
    end
 end
 
% IC and BCs.  % % In this program not theta but T itself is taken to  be the non-dimensional temperature
T_w=0; T_inf=1;
T=zeros(jmax,imax);  T2=zeros(jmax,imax);
T(2:jmax-1,2:imax-1)=T_w;  T2(2:jmax-1,2:imax-1)=T_w; 
T(:,1)=T_inf;  T2(:,1)=T_inf; 
T(1,:)=T_w;  T2(1,:)=T_w;
T(:,imax)=T(:,imax-1);  T2(:,imax)=T2(:,imax-1);
T(jmax,:)=T_w; T2(jmax,:)=T_w; 
unsteadiness_nd=1;unsteadiness_nd2=1; n=0; N=0;
DTc=100;

while unsteadiness_nd2>=epsilon_st             % FOU Advection Scheme
 n=n+1;
 T_old_iter=T2;
     % Inner-Loop to Iterative solution (by GS method) at each time step
           for j=2:jmax-1 
           for i=2:imax-1               
                 b2(j,i)= (aP0(j,i)*T_old_iter(j,i));
           end
          end
    Error=1;
    while Error>=epsilon
        T_old=T2; 
        N=N+1;      
        for i=2:imax-1
               for j=2:jmax-1 
                  T2(j,i)=(aE(j,i)*T2(j,i+1))+((aW(j,i)))*(T2(j,i-1))+(aN(j,i)*T2(j+1,i))+(((aS(j,i)))*(T2(j-1,i)))+b2(j,i);
                  T2(j,i)=T2(j,i)/aP(j,i) ;
                end
        end
     Error=max(max(abs(T2-T_old)));   % parameter to check convergence.
    end
      T2(:,imax)= T2(:,imax-1);
      T2(jmax,:)=T_w; 
    unsteadiness2=max(max(abs(T2-T_old_iter)))/Dt;     % Steady state converegnce criterion.
    unsteadiness_nd2=unsteadiness2 *L*H/(alpha*DTc); 
    fprintf('Time step no.: %5d, Unsteadiness_nd = %8.4e \n', n , unsteadiness_nd2);
end
 
n=0;
while unsteadiness_nd >= epsilon_st                                % QUICK Advection scheme
 n=n+1;
 T_old_iter=T;
     % Inner-Loop to Iterative solution (by GS method) at each time step
         epsilon=0.0001;
          for j=2:jmax-1      % For deffered Temperatures
           for i=2:imax-1     
                  Tx_d_plus(j,1)=0;  
              if (i==imax-1)
                  Tx_d_plus(j,i)= 0; 
                else     
                  Tx_d_plus(j,i)=(w1q_x_plus(i) * T_old_iter(j,i+1))+((w2q_x_plus(i) - 1)*T_old_iter(j,i))+(w3q_x_plus(i) * T_old_iter(j,i-1));
                  
              end
           end
          end
          
         for j=2:jmax-1 
           for i=2:imax-1    
                Ty_d_plus(1,i)=0;       
            if (j==jmax-1)
                Ty_d_plus(j,i)=0;  
            else       
                  Ty_d_plus(j,i)=(w1q_y_plus(j) * T_old_iter(j+1,i))+((w2q_y_plus(j) - 1)*T_old_iter(j,i))+(w3q_y_plus(j) * T_old_iter(j-1,i));
                  
            end
           end
          end
          
          for j=2:jmax-1 
           for i=2:imax-1               
              Qadv_N(j,i)=((mx*Tx_d_plus(j,i-1))-((mx*Tx_d_plus(j,i))))*Dy(j);
              Qadv_N(j,i)=Qadv_N(j,i) + (((my*Ty_d_plus(j-1,i))-(my*Ty_d_plus(j,i)))*Dx(i));
              b(j,i)= Qadv_N(j,i)+(aP0(j,i)*T_old_iter(j,i));
           end
          end
                      
     Error=1;
    while Error >= epsilon
        T_old=T; 
        N=N+1;      
        for i=2:imax-1
               for j=2:jmax-1 
                  T(j,i)=(aE(j,i)*T(j,i+1))+((aW(j,i)))*(T(j,i-1))+(aN(j,i)*T(j+1,i))+(((aS(j,i)))*(T(j-1,i)))+b(j,i);
                  T(j,i)=T(j,i)/aP(j,i) ;
                end
        end
        Error=max(max(abs(T-T_old)));   % parameter to check convergence.
    end
      T(:,imax)= T(:,imax-1);
      T(jmax,:)=T_w; 
    unsteadiness=max(max(abs(T-T_old_iter)))/Dt;     % Steady state converegnce criterion.
    unsteadiness_nd=unsteadiness*L*H/(alpha*DTc); 
    fprintf('Time step no.: %5d, Unsteadiness_nd = %8.4e \n', n , unsteadiness_nd);
end


% ploting steady state tempeperature contour.
v=linspace(0,1,11);
figure; 
[C,h]=contourf(xc,yc,T,v);
clabel(C,h);
xlabel('X');
ylabel('Y');
title('Steady State Non Dimensional Temperature Contour for QUICK');
figure;
[C,h]=contourf(xc,yc,T2,v);
clabel(C,h);
xlabel('X');
ylabel('Y');
title('Steady State Non Dimensional Temperature Contour for FOU');
figure;
plot(T(:,13),yc,'r-s'); xlabel({'\theta'}); ylabel('Y'); title('Temperature profiles for QUICK'); grid on; hold on;
plot(T(:,26),yc,'g-*'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T(:,38),yc,'b-d'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T(:,50),yc,'k-o'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T(:,imax),yc,'y-p'); xlabel({'\theta'}); ylabel('Y'); legend('x/L=0.2','x/L=0.4','x/L=0.6','x/L=0.8','x/L=1') ; hold off;
figure;
plot(T2(:,13),yc,'r-s'); xlabel({'\theta'}); ylabel('Y'); title('Temperature profiles for FOU'); grid on; hold on;
plot(T2(:,26),yc,'g-*'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T2(:,38),yc,'b-d'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T2(:,50),yc,'k-o'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T2(:,imax),yc,'y-p'); xlabel({'\theta'}); ylabel('Y'); legend('x/L=0.2','x/L=0.4','x/L=0.6','x/L=0.8','x/L=1'); hold off;
