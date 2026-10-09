%solving diff equation generated using bond graph technique 

clc 
clear

%initial conditions 

x0 = [80,0]; %initial condition for momentum and displacement 

timeSpan = [0,10]; % 0 to 10s

[t,x] = ode45(@MassSpringDamper_Oct8,timeSpan,x0);

figure
plot(t,x(:,1))
xlabel('Time (s)')
ylabel('Momentum (kg*m*s^-1')
title('Momentum over Time')

figure
plot(t,x(:,2))
xlabel('Time (s)')
ylabel('Displacement (m)')
title('Displacement over Time')
