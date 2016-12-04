% *********************************************************************************************
% MATHEMATICAL APPROACH OF CFD DEVELOPMENT: Code developed by Prof.Atul Sharma,CFD Lab ,IIT Bombay.
% Distributed for Assignement # 1 of courese ME415: CFDHT
%**********************************************************************************************
clear
L1=1;L2=1;imax=12;jmax=12;Ei=linspace(0,1,imax-1);Fi=linspace(0,1,jmax-1);
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
    for j=1:jmax-1
plot(x(i),y(j),'m-s')
hold on
plot(xc(i),yc(j),'k-o')
    end
end

Dx(2:imax-1)=x(2:imax-1)-x(1:imax-2);
dx(1:imax-1)=xc(2:imax)-xc(1:imax-1);

Dy(2:jmax-1)=y(2:jmax-1)-y(1:jmax-2);
dy(1:jmax-1)=yc(2:jmax)-yc(1:jmax-1);

fprintf('\n************** ONE-DIMENSIONAL HEAT CONDUCTION ***************');
%STEP-1: User-Input
rho = 7750.0; cp = 500.0; k = 16.2;
T0=30; T_wb=100.0;T_inf=30.0;h=100;qw=10000;
Q_vol_gen=50000; epsilon_st=0.0001; Dt=1;
for i=2:imax-1
    for j=2:jmax-1
        Q_gen(i,j)=Q_vol_gen*Dx(i)*Dy(j); % Total Heat Generation
    end
end

%STEP-2: Coefficient of implicit LAEs
  for i=1:imax-1
      for j=1:jmax-1
      aE(i,j)=k*dy(j)/dx(i);
      aF(i,j)=k*dx(i)/dy(j);
      
      end      
 end
 for i=2:imax-1
     for j=2:jmax-1
         aP0(i,j)=rho*cp*Dx(i)*Dy(j)/Dt;
         aP(i,j)=aP0(i,j)+aE(i,j)+aE(i-1,j)+ aF(i,j) + aF(i,j-1);
     end
 end
%STEP-3: IC and BCs
T(2:imax,1:jmax)=T0; T(1:imax,1)=T_wb; %T(imax)=T_eb; 

%    T(imax)=(k*T(imax-1))+(h*dx(imax-1)*T_inf);
%    T(imax)=T(imax)/(k+h*dx(imax-1));
% xT=[xc' T];
%savematfile('Temperature_ex5_2_0.dat','-ascii','xT')
unsteadiness_nd=1; n=0; alpha=k/(rho*cp); DTc=T_wb-T_inf;
 
%==== Time-Marching for Implicit Unsteady State LAEs: START ====
while unsteadiness_nd>=epsilon_st
%while n<=49
    n=n+1;
    
    T(imax,1:jmax)=T(imax-1,1:jmax);
    T(1:imax,jmax)=T(1:imax,jmax-1)+ (qw/k)*dx(imax-1);
    T(1,1:jmax)=   k*T(2,1:jmax)+ h*dy(jmax-1)*T_inf;
    T(1,1:jmax)=  T(1,1:jmax)/((h*dy(jmax-1)+ k));
    T_old=T;
    for i=2:imax-1
        for j=2:jmax-1
             b(i,j)=aP0(i,j)*T_old(i,j)+Q_gen(i,j);
        end
    end
    % Inner-Loop for Iterative solution (by GS method) at eac+-*+-**+h time step
    epsilon=0.0001;   %Convergence Criterion
    N=0;  % Counter for iternation number
    Error=1; % some number greater than epsilon to start the while loop below
    while Error>=epsilon
        T_old_iter=T; % present iterative value stored as old one
        N=N+1; % increase in the iteration counter
        for i=2:imax-1
            for j=2:jmax-1
                  T(i,j)=aE(i,j)*T(i+1,j)+aE(i-1,j)*T(i-1,j)+ aF(i,j)*T(i,j+1)+ aF(i,j-1)*T(i,j-1) + b(i,j);
                  T(i,j)=T(i,j)/aP(i,j);
            end
        end
    Error=max(abs(T-T_old_iter)); % parameter to check convergence
    end
    unsteadiness=max(max(abs(T-T_old)))/Dt;
    unsteadiness_nd=unsteadiness*L1*L1*0.0001/(alpha*DTc); %STEP 5: Steady state converegnce criterion
    fprintf('Time step number. %5d, Unsteadiness_nd = %8.4e\n', n , unsteadiness_nd);
end
figure(2)
v=linspace(80,1600,20);
[C,h]=contourf(xc,yc,T,v);
clabel(C,h);
xlabel('X length(m)');
ylabel('Y length(m)');
title('Steady State Temperature Contour with Heat generation');
colorbar
% xT=[xc yc T];
%savematfile('Temperature_ex5_2_5.dat','-ascii','xT')

%save('Temperature_ex5_2_anal.dat','-ascii','xT')
