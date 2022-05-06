% CARLsim
cap_u_cs = 0.25;%0.2;
tau_u_cs = 25;%20;
tau_x_cs = 600;%140;
tau_d_cs = 5;%15;
g_cs = 3;%1.01;
u_cs = 0;
x_cs = 1;
A_cs = 200;%6.333;%200;
i_cs = 0;
ts_cs=0.5; % time step
% Keivan
cap_u_km = 0.25;
tau_u_km = 25;
tau_x_km = 600;
tau_d_km = 5;
g_km = 3;%6.333*3;
u_km = 0;
x_km = 1;
A_km = 0;
i_km = 0;
ts_km=1.0; % time step
Vm = 200;%-80;
% general
spk = 0;
first_spk_km = 10; % first spike start time
ISI=50;%;20; % inter spike interval in milliseconds
t_total = 600;%1000;
simdur = linspace(0,t_total,t_total);
u_all_cs=[]; u_all_km=[];
x_all_cs=[]; x_all_km=[];
i_all_cs=[]; i_all_km=[]; a_all_km=[];

for t=1:length(simdur)
    % CARLsim's methods
    for i=1:(1/ts_cs)
        u_cs=u_cs+ts_cs*((-u_cs/tau_u_cs)+(cap_u_cs*(1-u_cs)).*spk);
        x_cs=x_cs+ts_cs*(((1-x_cs)/tau_x_cs)-u_cs.*x_cs.*spk);
        i_cs=i_cs+ts_cs*(-i_cs/tau_d_cs+A_cs.*u_cs.*x_cs-spk)*g_cs;
    end   

    % Keivan's methods
    i_km = 0;
    for i=1:(1/ts_km)       
        u_km=u_km+ts_km*((-u_km/tau_u_km)+(cap_u_km*(1-u_km)).*spk);
        x_km=x_km+ts_km*(((1-x_km-A_km)/tau_x_km)-u_km.*x_km.*spk);
        A_km=A_km+ts_km*((-A_km/tau_d_km)+u_km.*x_km.*spk);
        i_km=i_km+ts_km*(g_km*A_km*Vm);
    end

    % create spikes at ISI rate
    if mod((t+first_spk_km),ISI)==0
        spk = 1;
    elseif t == 1
        %spk = 1;
    else
        spk = 0;
    end
    
	u_all_cs(end+1)=u_cs;x_all_cs(end+1)=x_cs;i_all_cs(end+1)=i_cs;
    u_all_km(end+1)=u_km;x_all_km(end+1)=x_km;a_all_km(end+1)=A_km;i_all_km(end+1)=i_km;
end

if 1
subplot(2,1,1);
%figure
plot(simdur,u_all_cs,'-',simdur,x_all_cs,'-',simdur,i_all_cs,'-','LineWidth',3);
legend('u parameter','x parameter','i parameter')
title("CARLsim Synaptic Current - U:"+cap_u_cs+"; t_u:"+tau_u_cs+"; t_x:"+tau_x_cs+"; t_d:"+tau_d_cs+"; g:"+g_cs+"; A:"+A_cs)
grid on;
end
if 1
subplot(2,1,2); 
%figure
plot(simdur,u_all_km,'-',simdur,x_all_km,'-',simdur,a_all_km,'-',simdur,i_all_km,'-','LineWidth',3);
legend('u parameter','x parameter','a parameter','i parameter')
title("Keivan Synaptic Current - U:"+cap_u_km+"; t_u:"+tau_u_km+"; t_x:"+tau_x_km+"; t_d:"+tau_d_km+"; g:"+g_km+"; V:"+Vm)
grid on;
end