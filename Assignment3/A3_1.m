clear
clc
L1=1;L2=1;
imax=12;jmax=12;
rho=1;Uo=1;
mu=0.01;
Re=100;
epsilon_st=0.0001;
epsilon_mass =    0.00000001;

Dx=L1/(imax-2);Dy=L2/(jmax-2);
Dt=0.01;

x=linspace(0,L1,imax-1);    % Coordinates of face centers
y=linspace(0,L2,jmax-1);

% Coordinates of Cell Centers; 
xc(1)=x(1); xc(imax)=x(imax-1); 
for i=2:imax-1
xc(i)=(x(i)+x(i-1))/2;
end
yc(1)=y(1); yc(jmax)=y(jmax-1); 
for i=2:jmax-1
yc(i)=(y(i)+y(i-1))/2;
end 


u=zeros(jmax,imax-1);  v=zeros(jmax-1,imax);  p=zeros(jmax,imax);
p_dash=zeros(jmax,imax); u_pred=zeros(jmax,imax-1);  v_pred=zeros(jmax-1,imax);
u_pred(:,1)=0;u_pred(:,imax-1)=0;u_pred(1,:)=0; u_pred(jmax,:)=1;
v_pred(:,1)=0;v_pred(:,imax)=0;  v_pred(1,:)=0; v_pred(jmax-1,:)=0;


aE=(Dt*rho*Dy)/Dx;  aW=(Dt*rho*Dy)/Dx;
aN=(Dt*rho*Dx)/Dy;  aS=(Dt*rho*Dx)/Dy;
aP=aE+aW+aN+aS;

unsteadiness_nd=1;
while unsteadiness_nd>=epsilon_st
     u_old=u_pred;
     v_old=v_pred;
    
    for j=2:jmax-1                               
       for i=1:imax-2 
         mxu(j,i)=(rho*(u_old(j,i)+u_old(j,i+1)))/2;
         axu(j,i)=(max(mxu(j,i),0)*u_old(j,i))+(min(mxu(j,i),0)*u_old(j,i+1));
         dxu(j,i)=mu*((u_old(j,i+1)-u_old(j,i))/Dx);
       end 
    end
    
    for j=1:jmax-1                               
       for i=2:imax-2 
         myu(j,i)=(rho*(v_old(j,i)+v_old(j,i+1)))/2;
         ayu(j,i)=(max(myu(j,i),0)*u_old(j,i))+(min(myu(j,i),0)*u_old(j+1,i));
         dyu(j,i)=mu*((u_old(j+1,i)-u_old(j,i))/Dy);
       end 
    end
    
    for j=2:jmax-1                                
       for i=2:imax-2 
         Au(j,i)=((axu(j,i)-axu(j,i-1))*Dy)+((ayu(j,i)-ayu(j-1,i))*Dx);
         Du(j,i)=((dxu(j,i)-dxu(j,i-1))*Dy)+((dyu(j,i)-dyu(j-1,i))*Dx);
         Su(j,i)=(p(j,i)-p(j,i+1))*Dy;
         u_pred(j,i)=u_old(j,i)+((Dt/(rho*Dx*Dy))*(Du(j,i)-Au(j,i)+Su(j,i))) ;
       end 
    end
    
    for j=2:jmax-2                                
       for i=1:imax-1 
         mxv(j,i)=(rho*(u_old(j,i)+u_old(j+1,i)))/2;
         axv(j,i)=(max(mxv(j,i),0)*v_old(j,i))+(min(mxv(j,i),0)*v_old(j,i+1));
         dxv(j,i)=mu*((v_old(j,i+1)-v_old(j,i))/Dx);
       end 
    end
    
    for j=1:jmax-2                               
       for i=2:imax-1 
         myv(j,i)=(rho*(v_old(j,i)+v_old(j+1,i)))/2;
         ayv(j,i)=(max(myv(j,i),0)*v_old(j,i))+(min(myv(j,i),0)*v_old(j+1,i));
         dyv(j,i)=mu*((v_old(j+1,i)-v_old(j,i))/Dy);
       end 
    end
    
    for j=2:jmax-2                               
       for i=2:imax-1 
         Av(j,i)=((axv(j,i)-axv(j,i-1))*Dy)+((ayv(j,i)-ayv(j-1,i))*Dx);
         Dv(j,i)=((dxv(j,i)-dxv(j,i-1))*Dy)+((dyv(j,i)-dyv(j-1,i))*Dx);
         Sv(j,i)=(p(j,i)-p(j+1,i))*Dx;
         v_pred(j,i)=v_old(j,i)+((Dt/(rho*Dx*Dy))*(Dv(j,i)-Av(j,i)+Sv(j,i))) ;
       end 
    end
         
 while (1)      % 
  for j=2:jmax-1
    for i=2:imax-1
        Div(j,i)=((u_pred(j,i)-u_pred(j,i-1))*Dy)+((v_pred(j,i)-v_pred(j-1,i))*Dx);
    end    
  end  
  if (max(max(abs(Div(j,i)))) <= epsilon_mass)     
      break;

  else
      Error=1;
    while(Error>=epsilon_st)
     P_old=p_dash; 
     for j=2:jmax-1
       for i=2:imax-1
        p_dash(j,i)=((aE*p_dash(j,i+1))+(aW*p_dash(j,i-1))+(aN*p_dash(j+1,i))+(aS*p_dash(j-1,i))-(Div(j,i)))/aP;
       end    
     end 
        p_dash(:,imax) = p_dash(:,imax-1); 
        p_dash(jmax,:) = p_dash(jmax-1,:);
        p_dash(:,1) = p_dash(:,2); 
        p_dash(1,:) = p_dash(2,:); 
        Error = max(max(abs(p_dash-P_old)));
    end
    for j=2:jmax-1                              
       for i=2:imax-2 
          u_pred(j,i)=u_pred(j,i)+ ((Dt/(rho*Dx))*(p_dash(j,i)-p_dash(j,i+1))) ;
       end 
     end
    
    for j=2:jmax-2                              
       for i=2:imax-1 
          v_pred(j,i)=v_pred(j,i)+ ((Dt/(rho*Dy))*(p_dash(j,i)-p_dash(j+1,i))) ;
       end 
    end  
  end   
      
 end 
        p = p + p_dash;
        p(:,imax) = p(:,imax-1); 
        p(jmax,:) = p(jmax-1,:);
        p(:,1) = p(:,2); 
        p(1,:) = p(2,:); 
 
u_pred(:,imax-1)=0;  
v_pred(jmax-1,:)=0; 
unsteadiness_nd=max((max(max(abs(u_pred-u_old)))),(max(max(abs(v_pred-v_old)))));  
unsteadiness_nd
end

vu=linspace(-0.2,1,20); vp=linspace(-0.02,0.02,20);
figure; 
[C,h]=contourf(x,yc,u_pred,vu); clabel(C,h); xlabel('X'); ylabel('Y'); title('Steady State U Velocity Contour');
figure;
[C,h]= contourf(xc,yc,p,vp); clabel(C,h); xlabel('X'); ylabel('Y'); title('Steady State Pressure Contour');

ypA = [1 0.9766 0.9688 0.9609 0.9531 0.8516 0.7344 0.6172 0.5 0.4531 0.2813 0.1719 0.1016 0.0703 0.0625 0.0547 0];           % Benchmark Results
ucA100 = [1 0.84123 0.78871 0.73722 0.68717 0.23151 0.00332 -0.136641 -0.20581 -0.2109 -0.15662 -0.1015 -0.063434 -0.04775 -0.04192 -0.03717 0];
vcA100 =[0 -0.05906 -0.07391 -0.08864 -0.10313 -0.16914 -0.22445 -0.24533 0.05454 0.17527 0.17507 0.16077 0.12317 0.1089 0.10091 0.09233 0];
xpA = [1 0.9688 0.9609 0.9531 0.9453 0.9063 0.8594 0.8047 0.5 0.2344 0.2266 0.1563 0.0938 0.0781 0.0703 0.0625 0];

figure;
plot((u_pred(:,imax/2)+u_pred(:,(imax/2)+1))/2,yc,'r-s'); xlabel('U Velocity','FontSize',13); ylabel('Y','FontSize',13); title('U'); title('Variation of U-Velocity along the Verticle Centerline'); grid on; hold on;

plot(ucA100,ypA,'k-o'); hold off;

legend('Code','Benchmark Result','Location','southeast') ;
figure;
plot(xc,(v_pred(jmax/2,:)+v_pred((jmax/2)+1,:))/2,'g-*'); xlabel('X','FontSize',13,'FontWeight','bold'); ylabel('V Velolcity','FontSize',13,'FontWeight','bold'); title('U'); title('Variation of V-Velocity along the Horizontal Centerline'); grid on; hold on;
plot(xpA,vcA100,'k-o'); hold off;

legend('Code','Benchmark Result') ;

figure;
u_pred(:,imax)=0; v_pred(jmax,:)=0;
quiver(xc,yc,u_pred,v_pred,1.25); xlabel('X','FontSize',13,'FontWeight','bold'); ylabel('Y','FontSize',13,'FontWeight','bold'); title('U'); title('Velocity Vector'); 
