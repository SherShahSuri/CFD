
% *********************************************************************************************
% Coefficient of LAEs based Solution methodology OF CFD DEVELOPMENT: Code developed by Prof.Atul Sharma,CFD Lab ,IIT Bombay.
% Distributed for Assignement # 1 of courese ME415: CFDHT
%**********************************************************************************************
clear
fprintf('\n************** ONE-DIMENSIONAL HEAT CONDUCTION ***************');
%STEP-1: User-Input
rho = 7750.0; cp = 500.0; k = 16.2;
L1=1;L2=1; imax=12;jmax=12;
T0=30; T_wb=100;T_sb=200;T_eb=300.0;T_nb=400; %T_inf=100.0;h=1000;
Q_vol_gen=50000; epsilon_st=0.0001; 

%STEP-2: Geometrical Parameter and Stability criterion based time-step
alpha=k/(rho*cp); DTc=T_nb-T_wb;
Dx = L1/(imax-2);Dy = L2/(jmax-2); Dt =0.99*0.25*(Dx*Dy/alpha);
Q_gen=Q_vol_gen*Dx*Dy; % Total Heat Generation

%STEP-3: IC and BCs
T(2:imax-1,2:imax-1)=T0; T(1,2:imax-1)=T_wb; T(imax,2:imax-1)=T_eb;T(2:imax-1,1)=T_sb;T(2:imax-1,imax)=T_nb; 
%fprintf(' %d ' , T);

x=linspace(0,L1,imax-1); % Coordinates of face centers
y=linspace(0,L2,imax-1);

% Coordinates of Cell Centers; NEEDED FOR PLOTTING
xc(1)=x(1); xc(imax)=x(imax-1); % 
for i=2:imax-1
xc(i)=(x(i)+x(i-1))/2;
end
yc(1)=y(1); yc(imax)=y(imax-1); % 
for i=2:imax-1
yc(i)=(y(i)+y(i-1))/2;
end

figure(1)
for i=1:imax-1
    for j=1:imax-1
plot(x(i),y(j),'m-s')
hold on
plot(xc(i),yc(j),'k-o')
%xT=[xc yc T];
%savematfile('Temperature_ex5_1_st_0.dat','-ascii','xT')
    end
end
unsteadiness_nd= 1;   n=0;
%==== Time-Marching for Explicit Unsteady State LAEs: START ====
while unsteadiness_nd>=epsilon_st
%while n<=100
    n=n+1;

%    T(imax)=(2*k*T(imax-1))+(h*Dx*T_inf);
%    T(imax)=T(imax)/(2*k+h*Dx);
        T_old=T;
    % STEP4: 4. Computation of conduction-flux and temperature
    for j=1:jmax-1
        for i=1: imax-1
         
         if ((i==1) &&(j==1))
             qx_old(i,j)=0;
             qy_old(i,j)=0;
         elseif ((i==1)&&((j<jmax-1)&&(j>1)))
             qx_old(i,j)=-k*(T_old(i+1,j)-T_old(i,j))/(Dx/2.0);
             qy_old(i,j)=-k*(T_old(i,j+1)-T_old(i,j))/Dy;
         elseif((i==1)&&(j==jmax-1))
             qx_old(i,j)=-k*(T_old(i+1,j)-T_old(i,j))/(Dx/2.0);
             qy_old(i,j)=0;
         end
         if(((i>1)&&(i<imax-1))&&((j==1)||(j==jmax-1)))
             qx_old(i,j)=-k*(T_old(i+1,j)-T_old(i,j))/Dx;
             qy_old(i,j)=-k*(T_old(i,j+1)-T_old(i,j))/(Dy/2.0);
         elseif(((i>1)&&(i<imax-1))&&((j>1)&&(j<jmax-1)))
             qx_old(i,j)=-k*(T_old(i+1,j)-T_old(i,j))/Dx;
             qy_old(i,j)=-k*(T_old(i,j+1)-T_old(i,j))/Dy; 
         end
         if((i==imax-1)&&(j==1))
             qx_old(i,j)=0;
             qy_old(i,j)=-k*(T_old(i,j+1)-T_old(i,j))/(Dy/2.0);
         elseif ((i==imax-1)&&((j>1)&&(j<jmax-1)))
             qx_old(i,j)=-k*(T_old(i+1,j)-T_old(i,j))/(Dx/2.0);
             qy_old(i,j)=-k*(T_old(i,j+1)-T_old(i,j))/Dy;
         elseif((i==imax-1)&&(j==jmax-1))
             qx_old(i,j)=-k*(T_old(i+1,j)-T_old(i,j))/(Dx/2.0);
             qy_old(i,j)=-k*(T_old(i,j+1)-T_old(i,j))/(Dy/2.0);
         end
        end
    end
   
    for i=2:imax-1
        for j=2:imax-1
        Q_cond_old(i,j)=((qx_old(i-1,j)-qx_old(i,j))*(Dy)) + ((qy_old(i,j-1)-qy_old(i,j))*(Dx));
        T(i,j)=T_old(i,j)+(Dt/(rho*cp*Dx*Dy))*(Q_cond_old(i,j)+Q_gen);
        end	
    end
    unsteadiness = max(max(abs(T-T_old)/Dt));
    unsteadiness_nd=unsteadiness*L1*L1/(alpha*DTc); %STEP 5: Steady state converegnce criterion
    fprintf('Time step no.: %5d, Unsteadiness = %8.4e \n', n , unsteadiness_nd);
    
end


%****************************** OUTPUT PLOTTING *******************************
figure(2)
v=linspace(100,500,15);
[C,h]=contourf(xc,yc,T,v);
clabel(C,h);
xlabel('X length(m)');
ylabel('Y length(m)');
title('Steady State Temperature Contour with Heat Generation');
colorbar
%xT=[xc yc T];
%savematfile('Temperature_ex5_1_10.dat','-ascii','xT')

%xT_ana=[x_ana' T_ana'];
%savematfile('Temperature_ex5_1_anal.dat','-ascii','xT')

