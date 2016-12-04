clear
%STEP-1: User-Input
rho = 1000.0; cp = 4180.0; 
L=1;H=1; imax=32;jmax=32;
u=1;v=1;
T0=50; T_sb=0;    T_wb=100;
Q_vol_gen=0; epsilon_st=0.000001; 

%STEP-2: Geometrical Parameter and Stability criterion based time-step
DTc=100;
k=0.6;
Dx = L/(imax-2);Dy = H/(jmax-2); %Dt=0.005;

Q_gen=Q_vol_gen*Dx*Dy; % Total Heat Generation
mx=rho*u; my=rho*v;
alpha=k/(rho*cp) ;

%STEP-3: IC and BCs
T(2:imax,2:jmax)=T0; T(1,2:jmax)=T_wb; T(2:imax,1)=T_sb; 
%fprintf(' %d ' , T);

x=linspace(0,L,imax-1); % Coordinates of face centers
y=linspace(0,H,jmax-1);

% Coordinates of Cell Centers; NEEDED FOR PLOTTING
xc(1)=x(1); xc(imax)=x(imax-1); % 
for i=2:imax-1
xc(i)=(x(i)+x(i-1))/2;
end
yc(1)=y(1); yc(jmax)=y(jmax-1); % 
for i = 2:jmax-1
yc(i) = (y(i)+y(i-1))/2;
end

figure(1)
for i=1:imax-1
    for j=1:jmax-1
plot(x(i),y(j),'m-s')
hold on
plot(xc(i),yc(j),'k-o')
%xT=[xc yc T];
%savematfile('Temperature_ex5_1_st_0.dat','-ascii','xT')
    end
end
unsteadiness_nd=  1;   n=0;
%==== Time-Marching for Explicit Unsteady State LAEs: START ====
while unsteadiness_nd>=epsilon_st
%while n<=100
        n=n+1;
        T_old=T;
        Dt =0.99*(1/((u/Dx)+(v/Dy)));
    % STEP4: 4. Computation of advection-flux and temperature
    
      for i=1:imax-1
           for j=2:jmax-1
              hx(i,j)=mx*cp*T_old(i,j);
           end
      end
     
      for i=2:imax-1
            for j=1:jmax-1
               hy(i,j)=my*cp*T_old(i,j);
            end
      end
      
    num=0;
    
    for i=2:imax-1
        for j=2:jmax-1
        A(i,j) = ((hx(i,j)-hx(i-1,j))*(Dy)) + ((hy(i,j)-hy(i,j-1))*(Dx));
        T(i,j) = T_old(i,j)-((Dt/(rho*cp*Dx*Dy))*A(i,j));
        num = num  +  ((T(i,j)-T_old(i,j))*(T(i,j)-T_old(i,j)));
        end	
    end
    T(imax,:)= T(imax-1,:);
    T(:,jmax)=T_wb;
    unsteadiness = max(max(abs(T-T_old)/Dt));
    unsteadiness_nd=unsteadiness*L*H/(alpha*DTc); %STEP 5: Steady state converegnce criterion
    fprintf('Time step no.: %5d, Unsteadiness_nd = %8.4e \n', n , unsteadiness_nd);
    
end
T1(2:imax,2:jmax)=T0; T1(1,2:jmax)=T_wb; T1(2:imax,1)=T_sb;
unsteadiness_nd1=1;
n=0;

while unsteadiness_nd1>=epsilon_st  %SOU Advection Scheme
%while n<=100
        n=n+1;
        T_old=T1;
        Dt=0.001;
    % STEP4: 4. Computation of advection-flux and temperature
      hx(1,2:jmax-1)=mx*cp*T_old(1,2:jmax-1);
      hx(imax-1,2:jmax-1)=mx*cp*T_old(imax-1,2:jmax-1);
      for i=2:imax-2
           for j=2:jmax-1
              
              if(i==2)
                   hx(i,j)=mx*cp*((2*T_old(i,j))-T_old(i-1,j));
               else
                   hx(i,j)=mx*cp*(((3*T_old(i+1,j)) + (6*T_old(i,j))-T_old(i-1,j))/8);   
              end
           end
      end
     hy(2:imax-1,1)=my*cp*T_old(2:imax-1,1);
     hy(2:imax-1,jmax-1)=my*cp*T_old(2:imax-1,jmax-1);
      for i=2:imax-1
            for j=2:jmax-2
               if(j==2)
                   hy(i,j)=my*cp*((2*T_old(i,j))-T_old(i,j-1));
               else
                   hy(i,j)=my*cp*(((3*T_old(i,j+1)) + (6*T_old(i,j))-T_old(i,j-1))/8);
               end
            end
      end
      
    num=0;
    
    for i=2:imax-1
        for j=2:jmax-1
        A(i,j) = ((hx(i,j)-hx(i-1,j))*(Dy)) + ((hy(i,j)-hy(i,j-1))*(Dx));
        T1(i,j) = T_old(i,j)-((Dt/(rho*cp*Dx*Dy))*A(i,j));
        %num = num  +  ((T(i,j)-T_old(i,j))*(T(i,j)-T_old(i,j)));
        end	
    end
    T1(imax,:)= T1(imax-1,:);
    T1(:,jmax)=T_wb;
    
    unsteadiness = max(max(abs(T1-T_old)/Dt));
    unsteadiness_nd1=unsteadiness*L*H/(DTc*alpha); %STEP 5: Steady state converegnce criterion
    fprintf('Time step no.: %5d, Unsteadiness_nd = %8.4e \n', n , unsteadiness_nd1);
    
end
T2(2:imax,2:jmax)=T0; T2(1,2:jmax)=T_wb; T2(2:imax,1)=T_sb;
unsteadiness_nd2=1;
n=0;
while unsteadiness_nd2>=epsilon_st        %QUICK advection scheme
%while n<=100
        n=n+1;
        T_old=T2;
        Dt=0.001;
    % STEP4: 4. Computation of advection-flux and temperature
      hx(1,2:jmax-1)=mx*cp*T_old(1,2:jmax-1);
      hx(imax-1,2:jmax-1)=mx*cp*T_old(imax-1,2:jmax-1);
      for i=2:imax-1
           for j=2:jmax-1
               if(i==2)
                   hx(i,j)=mx*cp*(((T_old(i+1,j)) + (3*T_old(i,j))- T_old(i-1,j))/3);
               else
                   hx(i,j)=mx*cp*(((3*T_old(i+1,j)) + (6*T_old(i,j))-T_old(i-1,j))/8);
               end
           end
      end
     hy(2:imax-1,1)=mx*cp*T_old(2:imax-1,1);
     hy(2:imax-1,jmax-1)=my*cp*T_old(2:imax-1,jmax-1);
      for i=2:imax-1
            for j=2:jmax-2
               if(j==2)
                   hy(i,j)=my*cp*(((T_old(i,j+1)) + (3*T_old(i,j))-T_old(i,j-1))/3);
               else
                   hy(i,j)=my*cp*(((3*T_old(i,j+1)) + (6*T_old(i,j))-T_old(i,j-1))/8);
               end
            end
      end
      
    num=0;
    
    for i=2:imax-1
        for j=2:jmax-1
        A(i,j) = ((hx(i,j)-hx(i-1,j))*(Dy)) + ((hy(i,j)-hy(i,j-1))*(Dx));
        T2(i,j) = T_old(i,j)-((Dt/(rho*cp*Dx*Dy))*A(i,j));
        %num = num  +  ((T(i,j)-T_old(i,j))*(T(i,j)-T_old(i,j)));
        end	
    end
    T2(imax,:)= T2(imax-1,:);
    T2(:,jmax)=T_wb;
    
    unsteadiness = max(max(abs(T2-T_old)/Dt));
    unsteadiness_nd2=unsteadiness*L*H/(DTc*alpha); %STEP 5: Steady state converegnce criterion
    fprintf('Time step no.: %5d, Unsteadiness_nd = %8.4e \n', n , unsteadiness_nd2);
    
end



%****************************** OUTPUT PLOTTING *******************************
figure;
v=linspace(-6,104,15);
[C,h]=contourf  (xc,yc,T',v);
clabel(C,h);
xlabel('X length(m)');
ylabel('Y length(m)');
zlabel('Temperature');
title('Steady State Temperature Contour Advection FOU');
colorbar;
figure;
plot((100-(T(:,16)+T(:,17))/2),yc,'r-s'); xlabel('Temperature(°C)'); ylabel('y'); title('Temperature profile at verticle Centerline for FOU'); grid on;
figure;
v=linspace(-6,104,15);
[C,h]=contourf  (xc,yc,T1',v);
clabel(C,h);
xlabel('X length(m)');
ylabel('Y length(m)');
zlabel('Temperature');
title('Steady State Temperature Contour Advection SOU');
colorbar;
figure;
plot((100-(T1(:,16)+T1(:,17))/2),yc,'r-s'); xlabel('Temperature(°C)'); ylabel('y'); title('Temperature profile at verticle Centerline for SOU'); grid on;
figure;
v=linspace(-28,104,18);
[C,h]=contourf  (xc,yc,T2',v);
clabel(C,h);
xlabel('X length(m)');
ylabel('Y length(m)');
zlabel('Temperature');
title('Steady State Temperature Contour Advection QUICK');
colorbar;
figure;
plot((100-(T2(:,16)+T2(:,17))/2),yc,'r-s'); xlabel('Temperature(°C)'); ylabel('y'); title('Temperature profile at verticle Centerline for QUICK'); grid on;
%xT=[xc yc T];
%savematfile('Temperature_ex5_1_10.dat','-ascii','xT')

%xT_ana=[x_ana' T_ana'];
%savematfile('Temperature_ex5_1_anal.dat','-ascii','xT')

