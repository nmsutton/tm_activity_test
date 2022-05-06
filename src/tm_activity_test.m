% CARLsim
cap_u_cs = 0.2;%0.1119;%0.135;%0.1119;
tau_u_cs = 20;%60;%24;
tau_x_cs = 140;%80;%40;%80;%40;%70;%200;%120;
tau_d_cs = 15;%24;%28;%20;%28;%38.5;%23;%18;%15;
g_cs = 1.01;%1.0238;%1.0315;%1.0325;
u_cs = 0;
x_cs = 1;
A_cs = 20;%20*200;%4;%2;
i_cs = 0;
% Keivan
cap_u_km = 0.2229;%0.25;%0.005;%0.2;%0.1119;%0.135;%0.1119;
tau_u_km = 18.86;%25;%100;%25;%60;%24;tf
tau_x_km = 678.4;%600;%1000;%500;%80;%40;%80;%40;%70;%200;%120;tr
tau_d_km = 5.966;%5;%3.0;%5;%24;%28;%20;%28;%38.5;%23;%18;%15;
g_km = 6.333*1.475;%*3;%12;%3;%1;%1.0238;%1.0315;%1.0325;
u_km = 0;
x_km = 1;
A_km = 0;%1;%20*200;%4;%2;
%w_km = 10;
i_km = 0;
Vm = 70;%-80;
pst = 0;
spk = 0;
first_spk_km = 0; % first spike start time
ts=0.5; % time step
ts_km=0.025;%.5;%0.025; % time step
steps=40;%2;%40;
t_total = 201;%1000;
simdur = linspace(0,t_total,t_total);
u_all_cs=[]; u_all_km=[];
x_all_cs=[]; x_all_km=[];
i_all_cs=[]; i_all_km=[]; a_all_km=[];

for t=1:length(simdur)
    % CARLsim's methods
	u_cs=u_cs+((-u_cs/tau_u_cs)+(cap_u_cs*(1-u_cs)).*spk);
	x_cs=x_cs+(((1-x_cs)/tau_x_cs)-u_cs.*x_cs.*spk);
	i_cs=i_cs+ts*(-i_cs/tau_d_cs+A_cs.*u_cs.*x_cs-spk);
	i_cs=i_cs+ts*(-i_cs/tau_d_cs+A_cs.*u_cs.*x_cs-spk);
	i_cs=i_cs*g_cs;

    if mod(t,200)==0
    %    spk = 1;
    else
    %    spk = 0;
    end
    
    % Keivan's methods
    u_km=u_km+((-u_km/tau_u_km)+(cap_u_km*(1-u_km)).*spk);
    x_km=x_km+(((1-x_km-A_km)/tau_x_km)-u_km.*x_km.*spk);
    A_km=A_km+((-A_km/tau_d_km)+u_km.*x_km.*spk);
    i_km=g_km*A_km*Vm;

    %A_km=A_km+ts*((-i_km/tau_d_km)+u_km.*x_km.*spk);
    %A_km=A_km+ts*((-i_km/tau_d_km)+u_km.*x_km.*spk);
    %i_km=i_km*g_km;
    %i_km=g_km*A_km*Vm;
    %i_km=i_km+ts*(g_km*A_km*Vm);
    %i_km=i_km+ts*(g_km*A_km*Vm);
    for s=1:length(steps)
    %    i_km=i_km+ts_km*(g_km*A_km*Vm);  
    %    u_km=u_km+ts*(((-u_km/tau_u_km)+(cap_u_km*(1-u_km)).*spk));
    %    x_km=x_km+ts*((((1-x_km-A_km)/tau_x_km)-u_km.*x_km.*spk));
    %    A_km=A_km+ts*(((-A_km/tau_d_km)+u_km.*x_km.*spk));
    %    i_km=i_km+ts*(g_km*A_km*Vm);        
    end
    %i_km=g_km*A_km*Vm;
    %deltat = pst-t;
    %g = g_km*A_km;
    %i_km=g*exp(deltat/tau_d_km)*Vm;

    if mod((t+first_spk_km),20)==0
        spk = 1;
        %pst = t;
    elseif t == 1
        spk = 1;
    else
        spk = 0;
    end
    if t == 50
    %    spk = 1;
    %    pst = t;
    else
    %    spk = 0; 
    end
    %{
    deltat = pst-t;
    tau1r = tau_d_km / (tau_d_km-tau_x_km);
    y_ = A_km * exp(deltat/tau_d_km);
    x_ = 1 + (x_km-1+tau1r*A_km)*exp(deltat/tau_x_km) - tau1r*y_;
    u_ = u_km * exp(deltat/tau_u_km);
    u_km = u_ + (cap_u_km * (1-u_));
    A_km = y_ + u_km*x_;
    x_km = x_ - u_km*x_;
    if 0
        i_km=g_km*A_km*Vm;
    end
    if 1
        g = g_km*A_km;
        i_km=g*exp(-deltat/tau_d_km)*Vm;
    end
    if mod(t,20)==0
        pst = t;
    end
    %pst = t;
    %}
    
	u_all_cs(end+1)=u_cs;x_all_cs(end+1)=x_cs;i_all_cs(end+1)=i_cs;
    u_all_km(end+1)=u_km;x_all_km(end+1)=x_km;a_all_km(end+1)=A_km;i_all_km(end+1)=i_km;
end

if 0
%subplot(2,1,1);
figure
plot(simdur,u_all_cs,'-',simdur,x_all_cs,'-',simdur,i_all_cs,'-','LineWidth',3);
%plot(simdur,u_all_cs,'-',simdur,x_all_cs,'-','LineWidth',3);
legend('u parameter','x parameter','i parameter')
title("CARLsim Synaptic Current - U:"+cap_u_cs+"; t_u:"+tau_u_cs+"; t_x:"+tau_x_cs+"; t_d:"+tau_d_cs+"; g:"+g_cs+"; A:"+A_cs)
grid on;
end
if 1
%subplot(2,1,2); 
figure
plot(simdur,u_all_km,'-',simdur,x_all_km,'-',simdur,a_all_km,'-',simdur,i_all_km,'-','LineWidth',3);
%plot(simdur,u_all_km,'-',simdur,x_all_km,'-',simdur,a_all_km,'-','LineWidth',3);
%plot(simdur,u_all_km,'-',simdur,x_all_km,'-',simdur,i_all_km,'-','LineWidth',3);
%plot(simdur,u_all_km,'-',simdur,x_all_km,'-','LineWidth',3);
%plot(simdur,i_all_km,'-','LineWidth',3);
legend('u parameter','x parameter','a parameter','i parameter')
title("Keivan Synaptic Current - U:"+cap_u_km+"; t_u:"+tau_u_km+"; t_x:"+tau_x_km+"; t_d:"+tau_d_km+"; g:"+g_km)
grid on;
end