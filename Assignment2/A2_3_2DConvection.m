clear
% User-Input
rho = 10.0;   cp = 1.0; k=1; mu=1; 
L=6; H=1;
imax=62; jmax=22;
u=1; v=0; alpha=k/(rho*cp); Re=(rho*H*u)/mu; Pr=(mu*cp)/k;
mx=rho*u; my=rho*v;

epsilon_st=0.000001; 

% Geometrical Parameter and Stability criterion based time-step 
Dx = L/(imax-2); 
Dy = H/(jmax-2);
DTc=100;
Dt1=0.99/((u/Dx)+(v/Dy));
Dt2=0.5*0.99*Dx*Dx*Dy*Dy/(alpha*((Dx*Dx)+(Dy*Dy)));
Dt=0.9* min(Dt1,Dt2);

% IC and BCs
T_w=0; T_inf=1;
T=zeros(jmax,imax);  T2=zeros(jmax,imax);
T(2:jmax-1,2:imax-1)=T_w;  T2(2:jmax-1,2:imax-1)=T_w; 
T(:,1)=T_inf;  T2(:,1)=T_inf; 
T(1,:)=T_w;  T2(1,:)=T_w;
T(:,imax)=T(:,imax-1);  T2(:,imax)=T2(:,imax-1);
T(jmax,:)=T_w; T2(jmax,:)=T_w; 

x=linspace(0,L,imax-1);    % Coordinates of face centers
y=linspace(0,H,jmax-1); 

% Coordinates of Cell Centers; NEEDED FOR PLOTTING.
xc(1)=x(1); xc(imax)=x(imax-1); 
for i=2:imax-1
xc(i)=(x(i)+x(i-1))/2;
end
yc(1)=y(1); yc(jmax)=y(jmax-1); 
for i=2:jmax-1
yc(i)=(y(i)+y(i-1))/2;
end 

% Grid points ploting
figure;
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
 
unsteadiness_nd=1; unsteadiness_nd1=1; unsteadiness_nd2=1; n=0;
% Time-Marching for Explicit Unsteady State LAEs
while unsteadiness_nd>=epsilon_st             % FOU advection scheme
    n=n+1;
    T_old=T;
    w1=0; w2=1; w3=0;
    % Computation of conduction-flux and temperature.
    for j=2:jmax-1
       for i=1:imax-1 
         if(i==1)|(i==imax-1)
            qx_old(j,i)=-k*(T_old(j,i+1)-T_old(j,i))/(Dx/2.0);
            h_x(j,i)=((mx* T_old(j,i)))*cp;
        else
            qx_old(j,i)=-k*(T_old(j,i+1)-T_old(j,i))/Dx;
            h_x(j,i)=((mx *((w1*T_old(j,i+1))+(w2*T_old(j,i))+(w3*T_old(j,i-1)))))*cp;
        end
       end 
    end
    
    for j=1:jmax-1
       for i=2:imax-1 
         if(j==1)|(j==jmax-1)
            qy_old(j,i)=-k*(T_old(j+1,i)-T_old(j,i))/(Dy/2.0);
            h_y(j,i)=((my *T_old(j,i)))*cp;
        else
            qy_old(j,i)=-k*(T_old(j+1,i)-T_old(j,i))/Dy;
            h_y(j,i)=((my*((w1*T_old(j+1,i))+(w2*T_old(j,i))+(w3*T_old(j-1,i)))))*cp;
        end
       end 
    end
 
  for j=2:jmax-1
    for i=2:imax-1
        Qcond(j,i)=((qx_old(j,i-1)-qx_old(j,i))*Dy)+((qy_old(j-1,i)-qy_old(j,i))*Dx);
        Qadv(j,i)=((h_x(j,i)-h_x(j,i-1))*Dy)+((h_y(j,i)-h_y(j-1,i))*Dx);
        T(j,i)=T_old(j,i)+(((Dt/(rho*cp*Dx*Dy))) * (Qcond(j,i)-Qadv(j,i)));
    end    
  end        
  T(:,imax)= T(:,imax-1);
  T(jmax,:)=T_w;
  unsteadiness=max((max(abs(T-T_old))))/Dt;  % Steady state converegnce criterion.
  unsteadiness_nd=unsteadiness*L*H/(alpha*DTc); 
  fprintf('Time step no.: %5d, Unsteadiness_nd = %8.4e \n', n , unsteadiness_nd);
end
n=0;
while unsteadiness_nd2>=epsilon_st          % QUICK advection scheme
    n=n+1;
    T_old=T2;
    w1=3/8; w2=6/8; w3=-1/8;
    w11=1/3; w22=1; w33=-1/3;
    % Computation of conduction-flux and temperature.
    for j=2:jmax-1
       for i=1:imax-1 
         if(i==1)|(i==imax-1)
            h_x(j,i)=((mx* T_old(j,i)))*cp;
            qx_old(j,i)=-k*(T_old(j,i+1)-T_old(j,i))/(Dx/2.0);
          elseif (i==2)
            h_x(j,i)=((mx *((w11*T_old(j,i+1))+(w22*T_old(j,i))+(w33*T_old(j,i-1)))))*cp;
            qx_old(j,i)=-k*(T_old(j,i+1)-T_old(j,i))/(Dx);
          else
            h_x(j,i)=((mx *((w1*T_old(j,i+1))+(w2*T_old(j,i))+(w3*T_old(j,i-1)))))*cp;
            qx_old(j,i)=-k*(T_old(j,i+1)-T_old(j,i))/(Dx);
         end
       end 
    end
    
    for j=1:jmax-1
       for i=2:imax-1 
         if(j==1)|(j==jmax-1)
            h_y(j,i)=((my *T_old(j,i)))*cp;
            qy_old(j,i)=-k*(T_old(j+1,i)-T_old(j,i))/(Dy/2.0);
          elseif (j==2)
            h_y(j,i)=((my*((w11*T_old(j+1,i))+(w22*T_old(j,i))+(w33*T_old(j-1,i)))))*cp;
            qy_old(j,i)=-k*(T_old(j+1,i)-T_old(j,i))/(Dy);
         else
            h_y(j,i)=((my*((w1*T_old(j+1,i))+(w2*T_old(j,i))+(w3*T_old(j-1,i)))))*cp;
            qy_old(j,i)=-k*(T_old(j+1,i)-T_old(j,i))/(Dy);
        end
       end 
    end
 
  for j=2:jmax-1
    for i=2:imax-1
        Qcond(j,i)=((qx_old(j,i-1)-qx_old(j,i))*Dy)+((qy_old(j-1,i)-qy_old(j,i))*Dx);
        Qadv(j,i)=((h_x(j,i)-h_x(j,i-1))*Dy)+((h_y(j,i)-h_y(j-1,i))*Dx);
        T2(j,i)=T_old(j,i)+(((Dt/(rho*cp*Dx*Dy))) * (Qcond(j,i)-Qadv(j,i)));
    end    
  end        
  T2(:,imax)= T2(:,imax-1);
  T2(jmax,:)=T_w; 
  unsteadiness2=max(max(abs(T2-T_old)))/Dt;  % Steady state converegnce criterion.
  unsteadiness_nd2=unsteadiness2*L*H/(alpha*DTc); 
  fprintf('Time step no.: %5d, Unsteadiness_nd = %8.4e \n', n , unsteadiness_nd2);
end

% ploting steady state tempeperature contour.
v=linspace(0,1,11);
figure; 
[C,h]=contourf(xc,yc,T,v);
clabel(C,h);
xlabel('X');
ylabel('Y');
title('Steady State Temperature Contour for FOU');
figure;
[C,h]=contourf(xc,yc,T2,v);
clabel(C,h);
xlabel('X');
ylabel('Y');
title('Steady State Temperature Contour for QUICK');
figure;
plot(T(:,13),yc,'g-*'); xlabel({'\theta'}); ylabel('Y'); title('Temperature profiles for FOU'); grid on; hold on;
plot(T(:,26),yc,'b-d'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T(:,38),yc,'k-o'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T(:,50),yc,'y-p'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T(:,imax),yc,'r-s'); xlabel({'\theta'}); ylabel('Y'); legend('x/L=0.2','x/L=0.4','x/L=0.6','x/L=0.8','x/L=1') ; hold off;
figure;
plot(T2(:,13),yc,'g-*'); xlabel({'\theta'}); ylabel('Y'); title('Temperature profiles for QUICK'); grid on; hold on;
plot(T2(:,26),yc,'b-d'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T2(:,38),yc,'k-o'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T2(:,50),yc,'y-p'); xlabel({'\theta'}); ylabel('Y'); hold on;
plot(T2(:,imax),yc,'r-s'); xlabel({'\theta'}); ylabel('Y'); legend('x/L=0.2','x/L=0.4','x/L=0.6','x/L=0.8','x/L=1'); hold off;
