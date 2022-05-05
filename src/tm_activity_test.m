% CARLsim
cap_u_cs = 0.1119;%0.135;%0.1119;
tau_u_cs = 60;%24;
tau_x_cs = 80;%40;%80;%40;%70;%200;%120;
tau_d_cs = 24;%28;%20;%28;%38.5;%23;%18;%15;
g_cs = 1.01;%1.0238;%1.0315;%1.0325;
u_cs = 0;
x_cs = 1;
A_cs = 100;%20*200;%4;%2;
i_cs = 0;
% Keivan
cap_u_km = 0.005;%0.2;%0.1119;%0.135;%0.1119;
tau_u_km = 100;%25;%60;%24;
tau_x_km = 1000;%500;%80;%40;%80;%40;%70;%200;%120;
tau_d_km = 3.0;%5;%24;%28;%20;%28;%38.5;%23;%18;%15;
g_km = 1;%1.0238;%1.0315;%1.0325;
u_km = 0;
x_km = 1;
A_km = 0;%1;%20*200;%4;%2;
w_km = 10;
i_km = 0;
Vm = 70;
pst = 0;
spk = 1;
ts=0.5; % time step
t_total = 1000;
simdur = linspace(0,t_total,t_total);
u_all_cs=[]; u_all_km=[];
x_all_cs=[]; x_all_km=[];
i_all_cs=[]; i_all_km=[];

for t=1:length(simdur)
    % CARLsim's methods
	u_cs=u_cs+((-u_cs/tau_u_cs)+(cap_u_cs*(1-u_cs)).*spk);
	x_cs=x_cs+(((1-x_cs)/tau_x_cs)-u_cs.*x_cs.*spk);
	i_cs=i_cs+ts*(-i_cs/tau_d_cs+A_cs.*u_cs.*x_cs-spk);
	i_cs=i_cs+ts*(-i_cs/tau_d_cs+A_cs.*u_cs.*x_cs-spk);
	i_cs=i_cs*g_cs;
    
    % Keivan's methods
    deltat = pst-t;
    tau1r = tau_d_km / (tau_d_km-tau_x_km);
    y_ = A_km * exp(deltat/tau_d_km);
    x_ = 1 + (x_km-1+tau1r*A_km)*exp(deltat/tau_x_km) - tau1r*y_;
    u_ = u_km * exp(deltat/tau_u_km);
    u_km = u_ + (cap_u_km * (1-u_));
    A_km = y_ + u_km*x_;
    x_km = x_ - u_km*x_;
    i_km=g_km*A_km*Vm;
    pst = t;
    
	u_all_cs(end+1)=u_cs;x_all_cs(end+1)=x_cs;i_all_cs(end+1)=i_cs;
    u_all_km(end+1)=u_km;x_all_km(end+1)=x_km;i_all_km(end+1)=i_km;
end

%subplot(2,1,1);
figure
%plot(simdur,u_all_cs,'-',simdur,x_all_cs,'-',simdur,i_all_cs,'-','LineWidth',3);
%legend('u parameter','x parameter','i parameter')
%title("CARLsim Synaptic Current - U:"+cap_u_cs+"; t_u:"+tau_u_cs+"; t_x:"+tau_x_cs+"; t_d:"+tau_d_cs+"; g:"+g_cs+"; A:"+A_cs)
%grid on;
%subplot(2,1,2); 
plot(simdur,u_all_km,'-',simdur,x_all_km,'-',simdur,i_all_km,'-','LineWidth',3);
legend('u parameter','x parameter','i parameter')
title("Keivan Synaptic Current - U:"+cap_u_km+"; t_u:"+tau_u_km+"; t_x:"+tau_x_km+"; t_d:"+tau_d_km+"; g:"+g_km)
grid on;