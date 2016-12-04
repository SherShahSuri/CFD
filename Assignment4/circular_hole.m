clc
clear all
close all;

L=10; D=4; 
rho=7750; cp=500;  k=16.2; 
T_hole=1; T_outer_boundary=0;
dt=100;i_n=9; j_m=4; epsilon_st=0.000001; dz_ep=1/(i_n-1);  dz_tau=1/(j_m-1);
x=zeros(j_m,i_n); y=zeros(j_m,i_n);
    for j=1:i_n
        x(1,j)= 5+2*cos(-(j-1)*2*pi/(i_n-1));
        y(1,j)= 5+2*sin(-(j-1)*2*pi/(i_n-1));
    end
    for i=1:j_m
        x(i,i_n) = 5+((i-1)*((10/2)-2)/(j_m-1))+2;
        y(i,i_n) = 5;
    end
    x(:,1)=x(:,i_n); y(:,1)=5;
    x(j_m,i_n)=10; x(4,8)=10; x(4,7)=5; x(4,6)=0; x(4,5)=0; x(4,4)=0; x(4,3)=5; x(4,2)=10; x(4,1)=10;
    y(4,9)=5; y(4,8)=10; y(4,7)=10; y(4,6)=10; y(4,5)=5; y(4,4)=0; y(4,3)=0; y(4,2)=0; y(4,1)=5;
    xold=x; yold=y;
     error=1;
    while error>=epsilon_st
        for i=2:j_m-1
            for j=2:i_n-1
                  A=(((xold(i+1,j)-xold(i-1,j))/(2*dz_tau))^2)+(((yold(i+1,j)-yold(i-1,j))/(2*dz_tau))^2);
               
                  B=((xold(i,j+1)-xold(i,j-1))/(2*dz_ep))*((xold(i+1,j)-xold(i-1,j))/(2*dz_tau))+((yold(i,j+1)-yold(i,j-1))/(2*dz_ep))*((yold(i+1,j)-yold(i-1,j))/(2*dz_tau));
                
                  C=(((xold(i,j+1)-xold(i,j-1))/(2*dz_ep))^2)+(((yold(i,j+1)-yold(i,j-1))/(2*dz_ep))^2);
               
                x(i,j)=(-2*B*((x(i+1,j+1)-x(i-1,j+1)-x(i+1,j-1)+x(i-1,j-1))/(4*dz_ep*dz_tau))+A*((x(i,j+1)+x(i,j-1))/(dz_ep^2))+C*((x(i+1,j)+x(i-1,j))/(dz_tau^2)))/((2*A/(dz_ep^2))+(2*C/(dz_tau^2)));
                y(i,j)=(-2*B*((y(i+1,j+1)-y(i-1,j+1)-y(i+1,j-1)+y(i-1,j-1))/(4*dz_ep*dz_tau))+A*((y(i,j+1)+y(i,j-1))/(dz_ep^2))+C*((y(i+1,j)+y(i-1,j))/(dz_tau^2)))/((2*A/(dz_ep^2))+(2*C/(dz_tau^2)));
                
                error=max(max(abs((x-xold))+abs((y-yold))));
                xold=x; yold=y;
            end
        end     
    end
x_face=x; y_face=y;
for i=2:j_m
    for j=2:i_n
        x_cell(i,j)=(x(i-1,j)+x(i-1,j-1)+x(i,j-1)+x(i,j))/4;
        y_cell(i,j)=(y(i-1,j)+y(i-1,j-1)+y(i,j-1)+y(i,j))/4;
    end
end
for j=2:i_n
    x_cell(1,j)=(x(1,j-1)+x(1,j))/2;
    y_cell(1,j)=(y(1,j-1)+y(1,j))/2;
end
for j=2:i_n
    x_cell(j_m +1,j)=(x(j_m,j-1)+x(j_m,j))/2;
    y_cell(j_m+1,j)=(y(j_m,j-1)+y(j_m,j))/2;
end
x_cell(:,1)=x_cell(:,i_n); y_cell(:,1)=y_cell(:,i_n);
x_cell(:,10)=x_cell(:,2); y_cell(:,10)=y_cell(:,2);

% Temperatures.
T=zeros(j_m+1,i_n+1);
T(1,:)=T_hole; T(j_m+1,:)=T_outer_boundary;
Told=T;
Tnew=Told;
imax=4; jmax=9;
% Vertical faces Parameters.
dztau_e=zeros(imax,jmax); dzep_e=zeros(imax,jmax); dse_ztau=zeros(imax,jmax); dse_zep=zeros(imax,jmax);
for i=2:imax
    for j=1:jmax
        dztau_ex_=x_cell(i,j+1)-x_cell(i,j); 
        dztau_ey_=y_cell(i,j+1)-y_cell(i,j);
        dztau_e(i,j)=sqrt(dztau_ex_^2+dztau_ey_^2);
        ztau1_e=dztau_ex_/dztau_e(i,j); 
        ztau2_e=dztau_ey_/dztau_e(i,j);
        dzep_ex=x_face(i,j)-x_face(i-1,j);
        dzep_ey=y_face(i,j)-y_face(i-1,j); 
        dzep_e(i,j)=sqrt(dzep_ex^2+dzep_ey^2);
        zep1_e=dzep_ex/dzep_e(i,j);
        zep2_e=dzep_ey/dzep_e(i,j); 
        dse_y=-dzep_ex;
        dse_x=dzep_ey;
        dse_ztau(i,j)=(-dse_y*zep1_e+dse_x*zep2_e)/(ztau1_e*zep2_e-ztau2_e*zep1_e);
        dse_zep(i,j)=(dse_y*ztau1_e-dse_x*ztau2_e)/(ztau1_e*zep2_e-ztau2_e*zep1_e);
    end
end
% Horizontal faces Parameters.
dztau_n=zeros(imax,jmax); dzep_n=zeros(imax,jmax); dsn_ztau=zeros(imax,jmax); dsn_zep=zeros(imax,jmax);
for i=1:imax
    for j=2:jmax
        dztau_nx=x_cell(i+1,j)-x_cell(i,j); 
        dztau_ny=y_cell(i+1,j)-y_cell(i,j);
        dztau_n(i,j)=sqrt(dztau_nx^2+dztau_ny^2);
        ztau1_n=dztau_nx/dztau_n(i,j); 
        ztau2_n=dztau_ny/dztau_n(i,j);
        dzep_nx=x_face(i,j)-x_face(i,j-1); 
        dzep_ny=y_face(i,j)-y_face(i,j-1); 
        dzep_n(i,j)=sqrt(dzep_nx^2+dzep_ny^2);
        zep1_n=dzep_nx/dzep_n(i,j); 
        zep2_n=dzep_ny/dzep_n(i,j);
        dsn_y=dzep_nx; 
        dsn_x=-dzep_ny; 
        dsn_ztau(i,j)=(-dsn_y*zep1_n+dsn_x*zep2_n)/(ztau1_n*zep2_n-ztau2_n*zep1_n);
        dsn_zep(i,j)=(dsn_y*ztau1_n-dsn_x*ztau2_n)/(ztau1_n*zep2_n-ztau2_n*zep1_n);
    end
end
% Volume.
dv=zeros(imax+1,jmax+1);
for i=2:imax
    for j=2:jmax
        d1=sqrt((x_face(i-1,j)-x_face(i,j-1))^2+(y_face(i-1,j)-y_face(i,j-1))^2);
        d2=sqrt((x_face(i-1,j-1)-x_face(i,j))^2+(y_face(i-1,j-1)-y_face(i,j))^2);
        dv(i,j)=(1/4)*sqrt(4*(d1^2)*(d2^2)-(dzep_e(i,j)^2+dzep_e(i,j-1)^2-dzep_n(i,j)^2-dzep_n(i-1,j)^2));
    end
end
dv(:,1)=dv(:,9);dv(:,10)=dv(:,2);
error=1; count=0;
while error>=epsilon_st
% Vertex temperatures.
for i=2:imax-1
    for j=1:jmax
        Tv(i,j)=(dv(i,j+1)*Told(i,j+1)+dv(i,j)*Told(i,j)+dv(i+1,j)*Told(i+1,j)+dv(i+1,j+1)*Told(i+1,j+1))/(dv(i,j+1)+dv(i,j)+dv(i+1,j)+dv(i+1,j+1));
    end
end
Tv(1,:)=T_hole; Tv(imax,:)=T_outer_boundary;
for i=2:imax
    for j=2:jmax
        Q1zep=k*(((Told(i,j+1)-Told(i,j))*dse_ztau(i,j)/(dztau_e(i,j)))+((Told(i+1,j)-Told(i,j))*dsn_ztau(i,j)/(dztau_n(i,j)))-((Told(i,j)-Told(i,j-1))*dse_ztau(i,j-1)/(dztau_e(i,j-1)))-((Told(i,j)-Told(i-1,j))*dsn_ztau(i-1,j)/(dztau_n(i-1,j))));
        Q1ztau=k*(((Tv(i,j)-Tv(i-1,j))*dse_zep(i,j)/dzep_e(i,j))+((Tv(i,j)-Tv(i,j-1))*dsn_zep(i,j)/dzep_n(i,j))-((Tv(i,j-1)-Tv(i-1,j-1))*dse_zep(i,j-1)/dzep_e(i,j-1))-((Tv(i-1,j)-Tv(i-1,j-1))*dsn_zep(i-1,j)/dzep_n(i-1,j)));
        Qcond=Q1zep+Q1ztau;
        Tnew(i,j)=Told(i,j)+(dt/(rho*dv(i,j)*cp))*(Qcond);
    end
end 
Tnew(:,1)=Tnew(:,jmax); Tnew(:,jmax+1)=Tnew(:,2);
error=max(max(Tnew-Told));

Told=Tnew;
end

for i=2:imax-1
    for j=1:jmax
        Tv(i,j)=(dv(i,j+1)*Told(i,j+1)+dv(i,j)*Told(i,j)+dv(i+1,j)*Told(i+1,j)+dv(i+1,j+1)*Told(i+1,j+1))/(dv(i,j+1)+dv(i,j)+dv(i+1,j)+dv(i+1,j+1));
    end
end
Tv(1,:)=T_hole; Tv(imax,:)=T_outer_boundary;
figure;
v=linspace(0,1,15);
contourf(x_face,y_face,Tv);
[C,h]=contourf(x_face,y_face,Tv,v);
clabel(C,h);
colorbar;

% HeatTransfer
Q1=0;
for j=2:jmax
    Q1zep=-k*((Tnew(2,j)-Tnew(1,j))*dsn_ztau(1,j)/(dztau_n(1,j)));
    Q1ztau=-k*((Tv(1,j)-Tv(1,j-1))*dsn_zep(1,j)/dzep_n(1,j));
    Q1=Q1+Q1zep+Q1ztau;
end
disp(Q1);
Q2=0;
for j=2:jmax
    Q1zep=-k*(Tnew(imax+1,j)-Tnew(imax,j))*dsn_ztau(imax,j)/(dztau_n(imax,j));
    Q1ztau=-k*((Tv(imax,j)-Tv(imax,j-1))*dsn_zep(imax,j)/dzep_n(imax,j));
    Q2=Q2+Q1zep+Q1ztau;
end
disp(Q2);
