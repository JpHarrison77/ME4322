% Homework 1
% Six Bar linkage
% Static Equilibrium 

clc;
clear;

% define the joints 
A = [1.4 0.485 0];
B = [1.67 0.99 0];
C = [0.255 1.035 0];
D = [0.285 0.055 0];
E = [0.195 2.54 0];
F = [-0.98 2.57 0];
G = [0.05 0.2 0];

% define the length of bars 
lAB = norm(B-A);
lBC = norm(C-B);
lCD = norm(C-D);
lDE = norm(D-E);
lEF = norm(F-E);
lFG = norm(G-F);
lCE = norm(C-E);

% define weight of each link 
WAB = [0 -129.16415 0];
WBC = [0 -312.3123372 0];
WDE = [0 -538.8038514 0];
WEF = [0 -260.1114633 0];
WFG = [0 -566.1259767 0];

% define center of mass for each link
S1 = (A+B)/2;
S2 = (B+C)/2;
S3 = (D+C+E)/3;
S4 = (E+F)/2;
S5 = (F+G)/2;

% defining RPM of input link 

% given requirements 

NumofParts = 12500;
TimeOfPart = 9;

TimeInMin = TimeOfPart * 60; 

InputRPM = NumofParts / TimeInMin; 

InputRadperSec = InputRPM * (2*pi/60);

disp('Required Input Speed (RPM): ');
disp(InputRPM);
disp('Required Input Angular Velocity (rad/s): ');
disp(InputRadperSec);


syms FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin

ForceA = [FAx FAy 0];
ForceB = [FBx FBy 0];
ForceC = [FCx FCy 0];
ForceD = [FDx FDy 0];
ForceE = [FEx FEy 0];
ForceF = [FFx FFy 0];
ForceG = [FGx FGy 0];
inputTorque = [0 0 Tin];
AppliedForce = [0 -200 0];

% Static Equilibrium Conditions for Link AB
% Sum of forces = 0
% Fa + Fb + WeightofAB = 0

eqn1 = ForceA + ForceB + WAB == 0;

% Sum of moments = 0 with respect to COM of Link AB 
% S1A x FA + S1B x FB + InputTorque = 0

eqn2 = cross(A-S1,ForceA) + cross(B-S1,ForceB) + inputTorque == 0;

% Static Equilibrium Conditions for Link BC
% Sum of forces = 0
% Fb + Fc + WeightofBC = 0

eqn3 = -ForceB + ForceC + WBC == 0;

% Sum of moments = 0 with respect to COM of Link BC 
% S2B x FB + S2B x FC + WeightofBC = 0
eqn4 = cross(B-S2, -ForceB) + cross(C-S2, ForceC) == 0;

% static equilibrium conditions for link DE
% Sum of forces = 0
% -Fc + Fd + Fe + WeightofDE = 0

eqn5 = -ForceC + ForceD + ForceE + WDE == 0;

% Sum of moments = 0 with respect to COM of Link DE 
% S3C x FC + S3D x FD + S3E x FE + WeightofDE = 0
eqn6 = cross(C-S3, -ForceC) + cross(D-S3, ForceD) + cross(E-S3, ForceE) == 0;

% Sum of forces = 0 for Link EF
eqn7 = -ForceE + ForceF + WEF == 0;

% Sum of moments = 0 with respect to COM of Link EF 
% S4E x FE + S4F x FF = 0
eqn8 = cross(E-S4, -ForceE) + cross(F-S4, ForceF) == 0;

% Sum of forces = 0 for Link FG
eqn9 = -ForceF + ForceG + AppliedForce + WFG == 0;

% Sum of moments = 0 for Link FG
eqn10 = cross(F - S5, -ForceF) + cross(G - S5, ForceG) == 0;

% SOLVING THE 10 EQUATIONS

eqnMatrix = [eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10];

StaticSolution = solve(eqnMatrix,[FAx FAy FBx FBy FCx FCy FDx FDy FEx FEy FFx FFy FGx FGy Tin]);

Force_Ax = double(StaticSolution.FAx);
Force_Ay = double(StaticSolution.FAy);
Force_Bx = double(StaticSolution.FBx);
Force_By = double(StaticSolution.FBy);
Force_Cx = double(StaticSolution.FCx);
Force_Cy = double(StaticSolution.FCy);
Force_Dx = double(StaticSolution.FDx);
Force_Dy = double(StaticSolution.FDy);
Force_Ex = double(StaticSolution.FEx);
Force_Ey = double(StaticSolution.FEy);
Force_Fx = double(StaticSolution.FFx);
Force_Fy = double(StaticSolution.FFy);
Force_Gx = double(StaticSolution.FGx);
Force_Gy = double(StaticSolution.FGy);
InputTorque = double(StaticSolution.Tin);

% Display the results of the forces and input torque 
disp('Force A: ');
disp([Force_Ax, Force_Ay]);
disp('Force B: ');
disp([Force_Bx, Force_By]);
disp('Force C: ');
disp([Force_Cx, Force_Cy]);
disp('Force D: ');
disp([Force_Dx, Force_Dy]);
disp('Force E: ');
disp([Force_Ex, Force_Ey]);
disp('Force F: ');
disp([Force_Fx, Force_Fy]);
disp('Force G: ');
disp([Force_Gx, Force_Gy]);
disp('Input Torque: ');
disp(InputTorque);

% Angular Velocity Calculations

% Loop ABCDA
syms wBC wCD

omega_AB = [0 0 InputRadperSec];
omega_BC = [0 0 wBC];
omega_CD = [0 0 wCD];

eqn11 = cross(omega_AB, B - A) + cross(omega_BC, C - B) + cross(omega_CD, D - C) == 0;

loop1Solution = solve(eqn11(1:2) == [0; 0], [wBC wCD]);

angularVelocity_BC = double(loop1Solution.wBC);
angularVelocity_CD = double(loop1Solution.wCD);

omegaBC = [0 0 angularVelocity_BC];
omegaCD = [0 0 angularVelocity_CD];

% Second Loop D-C-E-F-G-D
syms wEF wFG

omega_EF = [0 0 wEF];
omega_FG = [0 0 wFG];

eqn12vec = cross(omegaCD, D - C) + cross(omega_EF, F - E) + cross(omega_FG, G - F) == 0;

solW = solve(eqn12vec(1:2) == [0; 0], [wEF wFG]);

angularVelocity_EF = double(solW.wEF);
angularVelocity_FG = double(solW.wFG);

omegaEF = [0 0 angularVelocity_EF];
omegaFG = [0 0 angularVelocity_FG];

% Angular Acceleration

% Loop 1 ABCDA
syms aBC aCD

alpha_AB = [0 0 0];
alpha_BC = [0 0 aBC];
alpha_CD = [0 0 aCD];

accLoop1 = cross(alpha_AB, B - A) + cross(omega_AB, cross(omega_AB, B - A)) + cross(alpha_BC, C - B) + cross(omegaBC, cross(omegaBC, C - B)) + cross(alpha_CD, D - C) + cross(omegaCD, cross(omegaCD, D - C));

loop1AccSolution = solve(accLoop1(1:2) == [0; 0], [aBC aCD]);

alphaBC = double(loop1AccSolution.aBC);
alphaCD = double(loop1AccSolution.aCD);

alphaBC_vector = [0 0 alphaBC];
alphaCD_vector = [0 0 alphaCD];

% Loop 2 D-C-E-F-G-D
syms aEF aFG

alpha_EF = [0 0 aEF];
alpha_FG = [0 0 aFG];

accLoop2 = cross(alphaCD_vector, E - D) + cross(omegaCD, cross(omegaCD, D - C)) + cross(alpha_EF, F - E) + cross(omegaEF, cross(omegaEF, F - E)) + cross(alpha_FG, G - F) + cross(omegaFG, cross(omegaFG, G - F));

loop2AccSolution = solve(accLoop2(1:2) == [0; 0], [aEF aFG]);

alphaEF = double(loop2AccSolution.aEF);
alphaFG = double(loop2AccSolution.aFG);

alphaEF_vector = [0 0 alphaEF];
alphaFG_vector = [0 0 alphaFG];

% Acceleration at CoM of each Link

aB_A = cross(alpha_AB, B - A) + cross(omega_AB, cross(omega_AB, B - A));

aS1_A = cross(alpha_AB, S1 - A) + cross(omega_AB, cross(omega_AB, S1 - A));

aS2_A = cross(alphaBC_vector, S2 - B) + cross(omegaBC, cross(omegaBC, S2 - B)) + aB_A;

aS3_D = cross(alphaCD_vector, S3 - D) + cross(omegaCD, cross(omegaCD, S3 - D));

aF_G = cross(alphaFG_vector, F - G) + cross(omegaFG, cross(omegaFG, F - G));

aS4_G = aF_G + cross(alphaEF_vector, S4 - F) + cross(omegaEF, cross(omegaEF, S4 - F));

aS5_G = cross(alphaFG_vector, S5 - G) + cross(omegaFG, cross(omegaFG, S5 - G));


% Velocity at Joint

vB_A = cross(omega_AB, B - A);
vC_B = cross(omegaBC, C - B);
vE_D = cross(omegaCD, E - D);
V_S4_F = cross(omegaEF, S4 - F);
V_F_G = cross(omegaFG, F - G);
vS4_G = V_S4_F + V_F_G;


% mass of Links in KG

MassAB = 13.16658; 
MassBC = 31.83612;
MassDCE = 54.92394;
MassEF = 26.51493;
MassFG = 57.70907;

% Mass MOMENT OF INTERTIA 

J_AB = 1;
J_BC = 1;
J_DCE = 1;
J_EF = 1;
J_FG = 1;

% Dynamic Equilibrium

syms NFAx NFAy NFBx NFBy NFCx NFCy NFDx NFDy NFEx NFEy NFFx NFFy NFGx NFGy NTin

NForceA = [NFAx NFAy 0];
NForceB = [NFBx NFBy 0];
NForceC = [NFCx NFCy 0];
NForceD = [NFDx NFDy 0];
NForceE = [NFEx NFEy 0];
NForceF = [NFFx NFFy 0];
NForceG = [NFGx NFGy 0];
NInputTorque = [0 0 NTin];

eqn15 = NForceA + NForceB + WAB == MassAB * aS1_A;
eqn16 = cross(A - S1, NForceA) + cross(B - S1, NForceB) + NInputTorque == J_AB * alpha_AB;

eqn17 = -NForceB + NForceC + WBC == MassBC * aS2_A;
eqn18 = cross(B - S2, NForceB) + cross(C - S2, NForceC) == J_BC * alphaBC_vector;

eqn19 = -NForceC + NForceD + NForceE + WDE == MassDCE * aS3_D;
eqn20 = cross(C - S3, NForceC) + cross(D - S3, NForceD) + cross(E - S3, NForceE) == J_DCE * alphaCD_vector;

eqn21 = -NForceE + NForceF + WEF == MassEF * aS4_G;
eqn22 = cross(E - S4, NForceE) + cross(F - S4, NForceF) == J_EF * alphaEF_vector;

eqn23 = -NForceF + NForceG + WFG + AppliedForce == MassFG * aS5_G;
eqn24 = cross(F - S5, NForceF) + cross(G - S5, NForceG) + cross(G - S5, AppliedForce) == J_FG * alphaFG_vector;

% solving equations
NeqnMatrix = [eqn15, eqn16, eqn17, eqn18, eqn19, eqn20, eqn21, eqn22, eqn23, eqn24];

DynamicSolution = solve(NeqnMatrix, [NFAx, NFAy, NFBx, NFBy, NFCx, NFCy, NFDx, NFDy, NFEx, NFEy, NFFx, NFFy, NFGx, NFGy, NTin]);

% Extract Forces from the dynamic Solution 
NForce_Ax = double(DynamicSolution.NFAx);
NForce_Ay = double(DynamicSolution.NFAy);
NForce_Bx = double(DynamicSolution.NFBx);
NForce_By = double(DynamicSolution.NFBy);
NForce_Cx = double(DynamicSolution.NFCx);
NForce_Cy = double(DynamicSolution.NFCy);
NForce_Dx = double(DynamicSolution.NFDx);
NForce_Dy = double(DynamicSolution.NFDy);
NForce_Ex = double(DynamicSolution.NFEx);
NForce_Ey = double(DynamicSolution.NFEy);
NForce_Fx = double(DynamicSolution.NFFx);
NForce_Fy = double(DynamicSolution.NFFy);
NForce_Gx = double(DynamicSolution.NFGx);
NForce_Gy = double(DynamicSolution.NFGy);
NInputTorque = double(DynamicSolution.NTin);

% Circle intersection technique

initial_theta = atan2(B(2) - A(2), B(1) - A(1));
if initial_theta < 0
    inputAngle = 2*pi + initial_theta;
else
    inputAngle = initial_theta;
end

thetaDeg = 0:360;

omegaAB_all = InputRadperSec * ones(size(thetaDeg));
omegaBC_all = nan(size(thetaDeg));
omegaCD_all = nan(size(thetaDeg));
omegaEF_all = nan(size(thetaDeg));
omegaFG_all = nan(size(thetaDeg));

new_B_x = nan(size(thetaDeg));
new_B_y = nan(size(thetaDeg));
new_C_x = nan(size(thetaDeg));
new_C_y = nan(size(thetaDeg));
new_E_x = nan(size(thetaDeg));
new_E_y = nan(size(thetaDeg));
new_F_x = nan(size(thetaDeg));
new_F_y = nan(size(thetaDeg));

for k = 1:numel(thetaDeg)

    theta = thetaDeg(k);

    % --- Joint B ---
    B_new = A + [lAB*cos(inputAngle + deg2rad(theta)), ...
                 lAB*sin(inputAngle + deg2rad(theta)), 0];

    % --- Joint C ---
    [Cx, Cy] = circcirc(B_new(1), B_new(2), lBC, D(1), D(2), lCD);

    if any(isnan(Cx)) || any(isnan(Cy))
        fprintf('C cannot be determined at %d degrees\n', theta);
        continue
    end

    C_candidates = [Cx(:), Cy(:), zeros(2,1)];
    [~, idxC] = min(vecnorm(C_candidates - C, 2, 2));
    C_new = C_candidates(idxC, :);

    % --- Joint E ---
    [Ex, Ey] = circcirc(B_new(1), B_new(2), lBC, C_new(1), C_new(2), lCE);

    if any(isnan(Ex)) || any(isnan(Ey))
        fprintf('E cannot be determined at %d degrees\n', theta);
        continue
    end

    E_candidates = [Ex(:), Ey(:), zeros(2,1)];
    [~, idxE] = min(vecnorm(E_candidates - E, 2, 2));
    E_new = E_candidates(idxE, :);

    % --- Joint F ---
    [Fx, Fy] = circcirc(E_new(1), E_new(2), lEF, G(1), G(2), lFG);

    if any(isnan(Fx)) || any(isnan(Fy))
        fprintf('F cannot be determined at %d degrees\n', theta);
        continue
    end

    F_candidates = [Fx(:), Fy(:), zeros(2,1)];
    [~, idxF] = min(vecnorm(F_candidates - F, 2, 2));
    F_new = F_candidates(idxF, :);

    % --- Solve angular velocities at this configuration ---
    syms wBC wCD wEF wFG

    eqnVel1 = cross([0 0 InputRadperSec], B_new - A) + ...
              cross([0 0 wBC], C_new - B_new) + ...
              cross([0 0 wCD], D - C_new) == 0;

    solVel1 = solve(eqnVel1(1:2), [wBC wCD]);
    omegaBC_all(k) = double(solVel1.wBC);
    omegaCD_all(k) = double(solVel1.wCD);

    eqnVel2 = cross([0 0 omegaCD_all(k)], E_new - C_new) + cross([0 0 wEF], F_new - E_new) + cross([0 0 wFG], G - F_new) == 0;

    solVel2 = solve(eqnVel2(1:2), [wEF wFG]);
    omegaEF_all(k) = double(solVel2.wEF);
    omegaFG_all(k) = double(solVel2.wFG);

    % --- Store results ---
    new_B_x(k) = B_new(1); 
    new_B_y(k) = B_new(2);
    new_C_x(k) = C_new(1); 
    new_C_y(k) = C_new(2);
    new_E_x(k) = E_new(1); 
    new_E_y(k) = E_new(2);
    new_F_x(k) = F_new(1); 
    new_F_y(k) = F_new(2);
end

figure;

h1 = plot(thetaDeg, omegaAB_all, 'k', 'LineWidth', 1.5); hold on;
h2 = plot(thetaDeg, omegaBC_all, 'r', 'LineWidth', 1.5);
h3 = plot(thetaDeg, omegaCD_all, 'b', 'LineWidth', 1.5);
h4 = plot(thetaDeg, omegaEF_all, 'g', 'LineWidth', 1.5);
h5 = plot(thetaDeg, omegaFG_all, 'm', 'LineWidth', 1.5);

grid on;
xlabel('Input Angle (deg)');
ylabel('Angular Velocity (rad/s)');
title('Angular Velocities of All Links');

legend([h1 h2 h3 h4 h5], {'AB', 'BC', 'CD / DCE', 'EF', 'FG'}, 'Location', 'best');