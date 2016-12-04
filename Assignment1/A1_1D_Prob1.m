
% *********************************************************************************************
% Coefficient of LAEs based Solution methodology OF CFD DEVELOPMENT: Code developed by Prof.Atul Sharma,CFD Lab ,IIT Bombay.
% Distributed for Assignement # 1 of courese ME415: CFDHT
%**********************************************************************************************
clear
fprintf('\n************** ONE-DIMENSIONAL HEAT CONDUCTION ***************');
%STEP-1: User-Input
rho = 7750.0; cp = 500.0; k = 16.2;
L=0.01; imax=12;
T0=30; T_wb=0.0;T_eb=100.0; %T_inf=100.0;h=1000;
Q_vol_gen=10^8; epsilon_st=0.0001; 

%STEP-2: Geometrical Parameter and Stability criterion based time-step
alpha=k/(rho*cp); DTc=T_eb-T_wb;
Dx = L/(imax-2);   Dt =0.99*0.5* (Dx*Dx/alpha); 
Q_gen=Q_vol_gen*Dx; % Total Heat Generation

%STEP-3: IC and BCs
T(2:imax-1)=T0; T(1)=T_wb; T(imax)=T_eb; 

x=linspace(0,L,imax-1); % Coordinates of face centers
y=zeros(1,imax-1); 

% Coordinates of Cell Centers; NEEDED FOR PLOTTING
xc(1)=x(1); xc(imax)=x(imax-1); % 
for i=2:imax-1
xc(i)=(x(i)+x(i-1))/2;
end
yc=zeros(1,imax); 
figure(1)
plot(x,y,'m-s')
hold on
plot(xc,yc','k-o')
xT=[xc T];
%savematfile('Temperature_ex5_1_st_0.dat','-ascii','xT')

  unsteadiness_nd=1;   n=0;
%==== Time-Marching for Explicit Unsteady State LAEs: START ====
while unsteadiness_nd>=epsilon_st
%while n<=100
    n=n+1;

%    T(imax)=(2*k*T(imax-1))+(h*Dx*T_inf);
%    T(imax)=T(imax)/(2*k+h*Dx);
        T_old=T;
    % STEP4: 4. Computation of conduction-flux and temperature
    for i=1:imax-1
         if(i==1)  ||(i==imax-1)
            qx_old(i)=-k*(T_old(i+1)-T_old(i))/(Dx/2.0);
        else
            qx_old(i)=-k*(T_old(i+1)-T_old(i))/Dx;
        end
    end
    for i=2:imax-1
        Q_cond_old(i)=(qx_old(i-1)-qx_old(i));
        T(i)=T_old(i)+(Dt/(rho*cp*Dx))*(Q_cond_old(i)+Q_gen);
    end	
    unsteadiness=max(abs(T-T_old))/Dt;
    unsteadiness_nd=unsteadiness*L*L/(alpha*DTc); %STEP 5: Steady state converegnce criterion
    fprintf('Time step no. %5d, Unsteadiness = %8.4e\n', n , unsteadiness);
end


%****************************** OUTPUT PLOTTING *******************************
figure(2)
hold on
title('TEMPERATURE PROFILE PLOT with Heat Generation','color','black','fontsize',7);
plot(xc,T,'mo');
xlabel('Length (m)'); % x-axis label
ylabel('Temperature (°C)'); % y-axis label

xT=[xc T];
%savematfile('Temperature_ex5_1_10.dat','-ascii','xT')
x_ana=linspace(0,L,100);
t1=(T_eb-T_wb).*x_ana./L; t2=Q_vol_gen.*x_ana.*(L-x_ana);t2=t2/(2*k);
T_ana=T_wb+t1+t2;
plot(x_ana,T_ana)
xT_ana=[x_ana' T_ana'];
legend('Numerical Solution','Analytical Solution');
%savematfile('Temperature_ex5_1_anal.dat','-ascii','xT')

