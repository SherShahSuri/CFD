% *********************************************************************************************
% MATHEMATICAL APPROACH OF CFD DEVELOPMENT: Code developed by Prof.Atul Sharma,CFD Lab ,IIT Bombay.
% Distributed for Assignement # 1 of courese ME415: CFDHT
%**********************************************************************************************
clear
L=0.01; imax=12;xi=linspace(0,1,imax-1);
Beta=1.2;
Beta_p1=Beta+1;Beta_m1=Beta-1;Beta_p1_div_m1=(Beta_p1/Beta_m1).^(2*xi-1);
num=(Beta_p1*Beta_p1_div_m1)-Beta_m1;den=2*(1+Beta_p1_div_m1);
x=L*num./den;
y=zeros(1,imax-1); 

xc(2:imax-1)=(x(2:imax-1)+x(1:imax-2))/2;
xc(1)=x(1);xc(imax)=x(imax-1);yc=zeros(1,imax);
figure(1)
plot(x,y,'m-s')
plot(xc,yc,'k-o')

Dx(2:imax-1)=x(2:imax-1)-x(1:imax-2);
dx(1:imax-1)=xc(2:imax)-xc(1:imax-1);

fprintf('\n************** ONE-DIMENSIONAL HEAT CONDUCTION ***************');
%STEP-1: User-Input
rho = 7750.0; cp = 500.0; k = 16.2;
T0=30; T_wb=0.0;T_inf=100.0;h=1000;
Q_vol_gen=10^8; epsilon_st=0.0001; Dt=1;
 Q_gen=Q_vol_gen.*Dx; % Total Heat Generation
%STEP-2: Coefficient of implicit LAEs
  for i=1:imax-1
     aE(i)=k/dx(i);
 end
 for i=2:imax-1
     aP0(i)=rho*cp*Dx(i)/Dt;
     aP(i)=  aP0(i)+aE(i)+aE(i-1);
 end
%STEP-3: IC and BCs
T(2:imax-1)=T0; T(1)=T_wb; %T(imax)=T_eb; 
%    T(imax)=(k*T(imax-1))+(h*dx(imax-1)*T_inf);
%    T(imax)=T(imax)/(k+h*dx(imax-1));
% xT=[xc' T];
%savematfile('Temperature_ex5_2_0.dat','-ascii','xT')
unsteadiness_nd=1; n=0; alpha=k/(rho*cp); DTc=T_inf-T_wb;
 
%==== Time-Marching for Implicit Unsteady State LAEs: START ====
while unsteadiness_nd>=epsilon_st
%while n<=49
    n=n+1;
    T(imax)=(k*T(imax-1))+(h*dx(imax-1)*T_inf);
    T(imax)=T(imax)/(k+h*dx(imax-1));
    T_old=T;
    for i=2:imax-1
         b(i)=aP0(i)*T_old(i)+Q_gen(i);
    end
    % Inner-Loop for Iterative solution (by GS method) at each time step
    epsilon=0.0001;   %Convergence Criterion
    N=0;  % Counter for iternation number
    Error=1; % some number greater than epsilon to start the while loop below
    while Error>=epsilon
        T_old_iter=T; % present iterative value stored as old one
        N=N+1; % increase in the iteration counter
        for i=2:imax-1
             T(i)=aE(i)*T(i+1)+aE(i-1)*T(i-1)+b(i);T(i)=T(i)/aP(i);
        end
    Error=max(abs(T-T_old_iter)); % parameter to check convergence
    end

    unsteadiness=max(abs(T-T_old))/Dt;
    unsteadiness_nd=unsteadiness*L*L/(alpha*DTc); %STEP 5: Steady state converegnce criterion
    fprintf('Time step number. %5d, Unsteadiness_nd = %8.4e\n', n , unsteadiness_nd);
end
figure(2)
plot(xc,T,'ko');
hold on
xT=[xc T];
title('Steady State Temperature profile with Heat generation');
xlabel('Length (m)'); 
ylabel('Temperature (°C)');
%savematfile('Temperature_ex5_2_5.dat','-ascii','xT')

x_ana=linspace(0,L,100);
t1=(T_wb-T_inf)*h.*x_ana./(k+h*L); 
t2=(2*k+h*L)/(k+h*L); t2=Q_vol_gen.*x_ana.*(t2*L-x_ana);t2=t2/(2*k);
T_ana=T_wb-t1+t2;
   plot(x_ana,T_ana)
xT_ana=[x_ana' T_ana'];
legend('Numerical Solution','Analytical Solution');
%save('Temperature_ex5_2_anal.dat','-ascii','xT')
