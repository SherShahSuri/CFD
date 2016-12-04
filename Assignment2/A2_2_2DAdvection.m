clear
% Grid Generation.
L=1; H=1;
imax=32; jmax=32;
xi=linspace(0,L,imax-1);
Beta=1.2;
Beta_p1=Beta+1;Beta_m1=Beta-1;Beta_p1_div_m1=(Beta_p1/Beta_m1).^(2*xi-1);
num=(Beta_p1*Beta_p1_div_m1)-Beta_m1;den=2*(1+Beta_p1_div_m1);
x=L*num./den;

yi=linspace(0,H,jmax-1);
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

Dx(2:imax-1)=x(2:imax-1)-x(1:imax-2);
dx(1:imax-1)=xc(2:imax)-xc(1:imax-1);
Dy(2:jmax-1)=y(2:jmax-1)-y(1:jmax-2);
dy(1:jmax-1)=yc(2:jmax)-yc(1:jmax-1);

 for i=2:imax-2          
     
                 w1s_x_plus(i)=0; w1q_x_plus(i)=((Dx(i))*((2*Dx(i))+Dx(i-1)))/((Dx(i)+Dx(i+1))*(Dx(i+1)+Dx(i-1)+2*Dx(i))); 
                w3s_x_plus(i)=(-1)*((Dx(i))/(Dx(i)+Dx(i-1)));  w3q_x_plus(i)=((Dx(i))*((Dx(i+1))))/((Dx(i)+Dx(i-1))*(Dx(i+1)+Dx(i-1)+2*Dx(i)));
                w2s_x_plus(i)=((2*Dx(i))+Dx(i-1))/(Dx(i)+Dx(i-1)); w2q_x_plus(i)=((Dx(i+1))*((2*Dx(i))+Dx(i-1)))/((Dx(i)+Dx(i+1))*(Dx(i-1)+Dx(i)));
     end
     
     for j=2:jmax-2 
        
                   w1s_y_plus(j)=0; w1q_y_plus(j)=((Dy(j))*((2*Dy(j))+Dy(j-1)))/((Dy(j)+Dy(j+1))*(Dy(j+1)+Dy(j-1)+2*Dy(j)));
                   w2s_y_plus(j)=((2*Dy(j))+Dy(j-1))/(Dy(j)+Dy(j-1));  w2q_y_plus(j)=((Dy(j+1))*((2*Dy(j))+Dy(j-1)))/((Dy(i)+Dy(j+1))*(Dy(j-1)+Dy(j)));
                   w3s_y_plus(j)=(-1)*((Dy(j))/(Dy(j)+Dy(j-1)));   w3q_y_plus(j)=((Dy(j))*((Dy(j+1))))/((Dy(j)+Dy(j-1))*(Dy(j+1)+Dy(j-1)+2*Dy(j)));
     end


% User Input.
rho = 1000.0; cp = 4180.0; k = 0.58;
u=1; v=1;
T0=50; T_wb=100.0; T_sb=0;
epsilon_st=0.000001; epsilon=0.0001; Dt=0.01;
mx=rho*u; my=rho*v;

% Coefficient of implicit LAEs

 for i=2:imax-1
    for j=2:jmax-1 
     aP0(j,i)=(rho*Dx(i)*Dy(j))/Dt;
    end
 end
 
 for i=2:imax-1
    for j=2:jmax-1 
     aP(j,i)=aP0(j,i)+(mx*Dy(j))+(my*Dx(i));
    end
 end
 
% IC and BCs. 
T=zeros(imax,jmax); T1=zeros(imax,jmax); T2=zeros(imax,jmax);
T(2:jmax,2:imax)=T0; T1(2:jmax,2:imax)=T0; T2(2:jmax,2:imax)=T0;
unsteadiness_nd=1;unsteadiness_nd1=1;unsteadiness_nd2=1; n=0; N=0; 
alpha=k/(rho*cp);
DTc=100;
T(:,1)=T_wb; T1(:,1)=T_wb; T2(:,1)=T_wb;
T(1,:)=T_sb; T1(1,:)=T_sb; T2(1,:)=T_sb;
n=0;
while unsteadiness_nd2>=epsilon_st                                       % FOU 
 n=n+1;
 T_old_iter=T2;
     % Inner-Loop to Iterative solution (by GS method) at each time step
    epsilon=0.0001;  
          
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
                  T2(j,i)=((mx*Dy(j)))*(T2(j,i-1))+(((my*Dx(i)))*(T2(j-1,i)))+b2(j,i);
                  T2(j,i)=T2(j,i)/aP(j,i) ;
                end
        end
         
        Error=max(max(abs(T2-T_old)));   % parameter to check convergence.
    end
      T2(:,imax)= T2(:,imax-1);
      T2(jmax,:)= T2(jmax-1,:);  
    unsteadiness2=max(max(abs(T2-T_old_iter)))/Dt;     % Steady state converegnce criterion.
    unsteadiness_nd2=unsteadiness2 *L*H/(alpha*DTc); 
    fprintf('Time step no.: %5d, Unsteadiness_nd = %8.4e \n', n , unsteadiness_nd2);
end
 
n=0;
while unsteadiness_nd >= epsilon_st                                %SOU
 n=n+1;
 T_old_iter=T;
     % Inner-Loop to Iterative solution (by GS method) at each time step
         epsilon=0.0001;
          for j=2:jmax-1 
           for i=2:imax-1     
              Tx_d_plus(j,1)=0;
              if (i==imax-1)
                Tx_d_plus(j,i)= 0; 
             else     
                  Tx_d_plus(j,i)=(w1s_x_plus(i) * T_old_iter(j,i+1))+((w2s_x_plus(i) - 1)*T_old_iter(j,i))+(w3s_x_plus(i) * T_old_iter(j,i-1));
                  
              end
           end
          end
          
          for j=2:jmax-1 
           for i=2:imax-1    
            Ty_d_plus(1,i)=0;     
            if (j==jmax-1)
                Ty_d_plus(j,i)= 0; 
             else       
                  Ty_d_plus(j,i)=(w1s_y_plus(j) * T_old_iter(j+1,i))+((w2s_y_plus(j) - 1)*T_old_iter(j,i))+(w3s_y_plus(j) * T_old_iter(j-1,i));
                 
             end
           end
          end
          
          for j=2:jmax-1 
           for i=2:imax-1               
              Qadv_N(j,i)=((mx*Tx_d_plus(j,i-1)) - ((mx*Tx_d_plus(j,i))))*Dy(j);
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
                  T(j,i)=((mx*Dy(j)))*(T(j,i-1))+(((my*Dx(i)))*(T(j-1,i)))+b(j,i);
                  T(j,i)=T(j,i)/aP(j,i) ;
                end
        end
        Error=max(max(abs(T-T_old)));   % parameter to check convergence.
    end
     T(:,imax)= T(:,imax-1);
     T(jmax,:)= T(jmax-1,:);  
    unsteadiness=max(max(abs(T-T_old_iter)))/Dt;     % Steady state converegnce criterion.
    unsteadiness_nd=unsteadiness*L*H/(alpha*DTc); 
    fprintf('Time step no.: %5d, Unsteadiness_nd = %8.4e \n', n , unsteadiness_nd);
end


%for i=2:imax-1;
 %           if (i==2)||(i==3)||(i==(imax-2))||(i==(imax-1))
  %              we1(i)=0;we2(i)=1;we3(i)=0; ww1(i)=0;ww2(i)=1;ww3(i)=0; wn1(i)=0;wn2(i)=1;wn3(i)=0; ws1(i)=0;ws2(i)=1;ws3(i)=0;
   %         else
    %            we1(i)=(Dx(i)*(2*Dx(i)+Dx(i-1)))/((Dx(i+1)+Dx(i))*(Dx(i+1) + 2*Dx(i) + Dx(i-1))); we2(i)=(Dx(i+1)*(2*Dx(i) + Dx(i-1)))/((Dx(i+1) + Dx(i))*(Dx(i) + Dx(i-1)));
     %           we3(i)=(Dx(i)*Dx(i+1))/((Dx(i-1) + Dx(i))*(Dx(i+1) + 2*Dx(i) + Dx(i-1))); ww1(i)=(Dx(i-1)*(2*Dx(i-1) + Dx(i-2)))/((Dx(i-1) + Dx(i))*(Dx(i) + 2*Dx(i-1) + Dx(i-2)));
      %          ww2(i)=(Dx(i)*(2*Dx(i-1) + Dx(i-2)))/((Dx(i) + Dx(i-1))*(Dx(i-2) + Dx(i-1))); ww3(i)=(Dx(i)*Dx(i-1))/((Dx(i-2) + Dx(i-1))*(Dx(i-2) + 2*Dx(i-1) + Dx(i)));
       %         wn1(i)=(Dy(i)*(2*Dy(i) + Dy(i-1)))/((Dy(i) + Dy(i+1))*(Dy(i+1) + 2*Dy(i) + Dy(i-1))); wn2(i)=(Dy(i+1)*(2*Dy(i) + Dy(i-1)))/((Dy(i) + Dy(i+1))*(Dy(i) + Dy(i-1)));
        %        wn3(i)=(Dy(i)*Dy(i+1))/((Dy(i) + Dy(i-1))*(Dy(i+1) + 2*Dy(i) + Dy(i-1))); ws1(i)=(Dy(i-1)*(2*Dy(i-1) + Dy(i-2)))/((Dy(i) + Dy(i-1))*(Dy(i) + 2*Dy(i-1) + Dy(i-2)));
         %       ws2(i)=(Dy(i)*(2*Dy(i-1) + Dy(i-2)))/((Dy(i) + Dy(i-1))*(Dy(i-2) + Dy(i-1))); ws3(i)=(Dy(i-1)*Dy(i))/((Dy(i-2) + Dy(i-1))*(Dy(i) + 2*Dy(i-1) + Dy(i-2 )));
          %  end
        %end   
        
        
  n=0;      
while unsteadiness_nd1>=epsilon_st                                      %QUICK
Told=T1; 
 % Calculation of the deferred Temperature terms
for i=2:imax-1
    for j=2:jmax-1
        %Using FOU  boundary values to simplicity
        if (i==2)||(i==3)||(i==(imax-2))||(i==(imax-1))||(j==2)||(j==3)||(j==(jmax-2))||(j==(jmax-1));
            Tde(i,j)=(w2q_x_plus(i-1)-1)*T(i,j); Tdw(i,j)=(w2q_x_plus(i-1)-1)*T(i-1,j);
            Tdn(i,j)=(w2q_y_plus(j-1)-1)*T(i,j); Tds(i,j)=(w2q_y_plus(j-1)-1)*T(i,j-1);
        else        
            Tde(i,j)=w1q_x_plus(i)*T(i+1,j)+(w2q_x_plus(i)-1)*T(i,j)+w3q_x_plus(i)*T(i-1,j); Tdw(i,j)=w1q_x_plus(i-1)*T(i,j)+(w2q_x_plus(i-1)-1)*T(i-1,j)+w3q_x_plus(i-1)*T(i-2,j);
            Tdn(i,j)=w1q_y_plus(j)*T(i,j+1)+(w2q_y_plus(j)-1)*T(i,j)+w3q_y_plus(j)*T(i,j-1); Tds(i,j)=w1q_y_plus(j-1)*T(i,j)+(w2q_y_plus(j-1)-1)*T(i,j-1)+w3q_y_plus(j-1)*T(i,j-2);   
        end
    end
end

%Computation of the Q advection term
 for i=2:imax-1
        for j=2:jmax-1
            Qadv=(mx*Tdw(i,j)-mx*Tde(i,j))*Dy(j)+(my*Tds(i,j)-my*Tdn(i,j))*Dx(i);
            b1(i,j)=aP0(i,j)*T(i,j)+Qadv;
        end
 end
 % Inner-Loop to Iterative solution (by GS method) at each time step
 epsilon=0.0001;   %Convergence Criterion
 N=0;  %Counter to iternation number
 Error=1; %some number greater than epsilon to start the loop below
    while Error>=epsilon
        T_old_iter=T1; % present iterative value stored as old one
        N=N+1; % increase in the iteration counter
        for i=2:imax-1
            for j=2:jmax-1
            T1(i,j)=((mx*Dy(j)))*T1(i-1,j)+((my*Dx(i)))*T1(i,j-1)+b1(i,j);
            T1(i,j)=T1(i,j)/aP(i,j);
            end
        end  
        Error=max(max(abs(T1-T_old_iter))); % parameter to check convergence
    end
    T1(imax,:)=T1(imax-1,:); T1(:,jmax)=T1(:,jmax-1);
    unsteadiness1=max(max(abs(T1-Told)));
    unsteadiness_nd1=unsteadiness1*L*H/(alpha*DTc); %STEP 5: Steady state converegnce criterion
    fprintf('Time step no.: %5d, Unsteadiness_nd = %8.4e \n', n , unsteadiness_nd1);
end


% ploting steady state tempeperature contour.
v=linspace(-5,103,16);
figure;
[C,h]=contourf(xc,yc,T,v);
clabel(C,h);
xlabel('X');
ylabel('Y');
title('Steady State Temperature Contour for SOU');

  Tc=T(:,(imax/2))/2;
figure;
plot(Tc,yc,'g-s');xlabel('Temperature(°C)'); ylabel('y'); title('Temperature profile at verticle Centerline for SOU'); grid on;

figure;
[C,h]=contourf(xc,yc,T1,v);
clabel(C,h);
xlabel('X');
ylabel('Y');
title('Steady State Temperature Contour for QUICK');
Tc1=T1(:,(imax/2))/2;
 
figure;
plot(Tc1,yc,'b-o'); xlabel('Temperature(°C)'); ylabel('y'); title('Temperature profile at verticle Centerline for QUICK'); grid on;

figure;
[C,h]=contourf(xc,yc,T2);
clabel(C,h);
xlabel('X');
ylabel('Y');
title('Steady State Temperature Contour for FOU');

Tc2=T2(:,(imax/2))/2;
   
figure;
plot(Tc2,yc,'r-s'); xlabel('Temperature(°C)'); ylabel('y'); title('Temperature profile at verticle Centerline for FOU'); grid on;


