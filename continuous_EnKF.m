function [xvec,x_estvec]=continuous_EnKF(mass,drag,ODm,holeD,Kc,rho,k,h1,n,ODj,h2,...
    ODc,hC,ODb,hB,aa,bb,Au1,Au2,AuC,AuB,Hook,num_members,x_tr,rNom,r,w,z,num_iterations,maxtime)

Zcov = eye(1); %for creating measurement noise covariance matrix
xvec        = [];
x_estvec    = [];
yvec        = [];

for j = 1:1
    Zcov(j,j) = z(j)^2;
end
p1=10;

tspan = [0 maxtime];
y0    = ones(1,10)*0.001;

for i = 1:num_iterations
    F_hook = Hook(i);
    disp(['Iterasi ke ',num2str(i)])
    if i == 1
        x_tr = y0;
        x_est = ones(num_members,10)*5;
    end
    
    [t,y_tr] = ode45(@(t,y)truefun_axial12nov19(t,y,mass,drag,ODm,holeD,Kc,rho,k,h1,n,ODj,h2,...
        ODc,hC,ODb,hB,aa,bb,Au1,Au2,AuC,AuB,F_hook,r), tspan, x_tr);
    x_tr = y_tr(end,:)+ w.*randn(1,p1); %compute true value of state at next time step
    C = [0 1 0 0 0 0 0 0 0 0];
    for j=1:num_members
        W(j,:)=w.*randn(1,p1);                          %create process noise
        Z(j,:)=z.*randn(1,1);                          %create measurement noise
        [t,y_est] = ode45(@(t,y)fun_axial12nov19(t,y,mass,drag,ODm,holeD,Kc,rho,k,h1,n,ODj,h2,...
    ODc,hC,ODb,hB,aa,bb,Au1,Au2,AuC,AuB,F_hook,rNom), tspan, x_est(j,:));
        x_est(j,:)=y_est(end,:)+ W(j,:);      %forecast state
        y(j,:)= C*x_tr' + Z(j,:);                 %make measurement
        y_for(j,:)=C*x_est(j,:)';              %forecast measurement
    end
    x_estbar=mean(x_est,1);
    y_forbar=mean(y_for,1);
    
    for j=1:p1
        Ex(:,j)= x_est(:,j)-x_estbar(j);
    end
    
    for j=1:1
        Ey(:,j)= y_for(:,j)-y_forbar(j);
    end
    
    Pxy = Ex'*Ey/(num_members-1);
    Pyy = Ey'*Ey/(num_members-1)+Zcov;                     %The addition of Zcov to Pyy is not done in Gillijns et. al but I use it here in case num_members=2 or Pyy is nearly singular
    K = Pxy*inv(Pyy);
    perbaharui = y-y_for;
    x_est = x_est+(K*perbaharui')';
    x_estbar=mean(x_est,1);
    xvec = [xvec; x_tr];
    x_estvec = [x_estvec; x_estbar];
    
end