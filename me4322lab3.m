
clc;
clear; 

% Define masses of each part of the system (in Kilograms):

% m of platform 
m1 = 0.57664;
% mass of main arm support 
m2 = 0.09792;
% mass of support arms 
m3 = 0.03312;
% mass of rack 
m4 = 0.00423;
% mass of translation lever 
m5 = 0.00118;

mp = 0.00110;

% define input weight (in Newtons): 

Fin = 845;

% define lengths of each part of the system (in Meters) 

l1 = 0.048;
l2 = 0.19; 
l3 = 0.064;
l4 = 0.127;
l5 = 0.02;
l6 = 0.025;
l7 = 0.12; 
rp = 0.006; 

% dampers: 
d1 = 0.01;
d2 = 0.01;
d3 = 0.001;
d4 = 0.001;


%spring
% assume all springs are made out of steel
G = 80 * 10^9;
ns = 9;
n1 = 6;
n2 = 35;
Ds = 0.001;
D1 = 0.002;
D2 = 0.00055;
dss = 0.008;
ds1 = 0.02;
ds2 = 0.000314;
ks = (Ds^4 * G) / (8*dss^3 * ns);
k1 = (D1^4 * G) / (8*ds1^3 * n1);
k2 = (D2^4 * G) / (8*ds2^3 * n2);

% lumping mass: 

j1 = 1/3 * m2 * l2^2;
j2 = 1/3 * m3 * l4^2;
%jpiv = 1/2 * m5 * (l6/l5)^2;
jp = 1/2 * mp * rp^2;

meq1 = j1 / l4^2;
jeq1 = meq1 * l1 ^2 + j2;
meq2 = (jeq1 / l2 ^ 2) * 2;
jeq2 = meq2 * l6^2;
meq3 = (jeq2 / l5 ^2) + m4;
jeq = meq3 * rp^2 + jp;

%jeq1 = ((j1)/4) + j2 * l3^2;
%meq1 = (1/l4)^2 * jeq1;
%jeq2 = ((j1)/4) + meq1 * (l2/2)^2;
%meq2 = (1/l2)^2 * jeq2;
%meq3 = 2 * meq2;
%jeq3 = l5 ^ 2 * meq3;
%meq4 = l6^2 * jeq3;
%jeq = meq4 + rp^2 + jp;

% lumping dampers
deq1 = d1 / l4^2;
deq2 = deq1 * (l2/2)^2 + d2; 
deq3 = (deq2 / l2^2);
deq4 = deq3 * (l5)^2;
deq5 = deq4 / l6^2;
deq = deq5 * (rp/l7)^2; 

% lumping springs 

keq1 = (ks/2) * l1^2;
keq2 = (ks/2) * l3^2;
keq3 = (keq2) / (l4^2);
keq4 = ((((l2 / 2) ^ 2) * keq3) + keq1) / (l2^2);
keq5 = 2 * keq4;
keq6 = keq5 + k1;
keq7 = keq6 * l6^2;
keq8 = (keq7/(l6 ^ 2)) + k2;
keq = keq8 * rp^2;

% lump force 
feq1 = (Fin/4)*((1/l2));
feq2 = (Fin/4)*((1/l4));
feq3 = (2 / (l2)) * feq2 + feq1;
feq4 = 2 * feq3;
teq1 = l5 * feq4;
feq5 = ((1/l6) * teq1)/l7;
teq = feq5 * rp;

syms x(t) % x(t) is the x in Mx'' + kx = F 

eqn1 = jeq*diff(x,t,2) + deq*diff(x,t) + keq*x == teq;

Dx = diff(x,t);

%specifying initial conditions 

initialCon = [x(0) == 0, Dx(0)==10];

% solving for X(t)
solutionX = dsolve(eqn1,initialCon);

% solving for V(t) - Velocity as function of time 
solutionV = diff(solutionX);

% solving for A(t) - Acceleration as function of time 
solutionA = diff(solutionV);

fprintf('Jeq = %.6f\n', jeq);
fprintf('Keq = %.6f\n', keq);
fprintf('Deq = %.6f\n', deq);
fprintf('teq = %.6f\n', teq);

%plotting X(t) between times 0 and 10s 
%fplot(solutionX,[0 10])

%fplot(solutionV,[0 10]);

%fplot(solutionA,[0 10]);

%laplaceTrans = laplace(eqn1)
%

% Plotting X(t) between times 0 and 10s.
fplot(solutionX, [0,10]); xlabel('Time (seconds)'); ylabel('Angle (degrees)'); title('Scale Position vs Time');

%fplot(solutionV, [0,10]); xlabel('Time (seconds)'); ylabel('Angular Velocity (degrees/sec)'); title('Scale Angular Velocity vs Time');


%fplot(solutionA, [0,10]); xlabel('Time (seconds)'); ylabel('Angular Acceleration (degrees/sec^2)'); title('Scale Angular Acceleration vs Time');

