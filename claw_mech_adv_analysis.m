clear; clc; close all
%%
% Trigger Dimensions
F_in = 10; %total input force from locking mechanism
y_trigger = 2.5; %position of force on trigger (in)
x_trigger = 1; %Position of contact to plates (in)

%%
%Each following sections perform the same static analysis but with varying 
%claw geometries. This code could easilly be cleaned up and condensed using
%a function or something similar. However, comparing the different configs 
%was an afterthought and I wanted to code it quickly. Doing this actually 
%helped us optimize our design quite a bit.
%%
%Original Dimensions(V1):
x_a = 0;
y_a_min = 1.3; %Minimum slider height (in)
y_a_max =3.6; %Maximum slider height (in)
BC = 2; %Distance between pin B and C (in) 
AC = 2.4; %Distance between pin A and C (in)
CD = 3.5; %Distance between pin C and D (in)

%%Static Analysis:
%Geometry:
y_a = (y_a_min:.001:y_a_max); %Slider height range vector
alpha = zeros(1,length(y_a));
beta = zeros(1,length(y_a));
MA = zeros(1,length(y_a));
for i = 1:length(y_a)
    BA = sqrt(y_a(i)^2+x_a^2);
    psi = asind(x_a/BA);
    theta = acosd((BC^2-BA^2-AC^2)/(-2*BA*AC));
    phi = asind(AC*sind(theta)/BC);
    alpha(i) = psi + phi;
    beta(i) = (90-alpha(i));
    y_c = BC*cosd(alpha(i));
    y_d = CD*sind(alpha(i))-y_c;
    x_c = BC*sind(alpha(i));
    
    %Forces
    F_trigger = F_in*y_trigger/x_trigger;
    F = F_trigger/3; %force through each of the 3 slider pinsC_y = -F
    A_x = F*(x_c-x_a)/(y_a(i)-y_c);
    C_x = -A_x;
    C_y = -F;
    F_ac = sqrt(A_x^2+F^2); %axial force through member AC
    D_x = (-C_x*y_c-C_y*x_c)/y_d;
    F_out = D_x*3; %total inward horizontal clamping force
    
    %Mechanical Advantage:
    %Clamp Trigger MA not currently included
    MA(i) = F_out/F_in;
end

% Mechanical Advantage (V1):
plot(beta,MA,'r')
hold on
grid on
xlabel('Claw Angle (degrees)')
ylabel('Mechanical Advantage')
title('Mechanical Advantage as a Function of Claw Angle:')
legend('(V1)')
xlim([0.0 100.0])
ylim([0.0 14.0])
hold on


%%
%LAD Configuration (V2):
%Claw Dimensions
x_a = 0;
BC = .9; %Distance between pin B and C (in) 
AC = 1.25; %Distance between pin A and C (in)
CD = 3.5; %Distance between pin C and D (in)
y_a_min = .85; %Minimum slider height (in)
y_a_max =2.09; %Maximum slider height (in)

%%Static Analysis:
%Geometry:
y_a = (y_a_min:.001:y_a_max); %Slider height range vector
alpha = zeros(1,length(y_a));
beta = zeros(1,length(y_a));
MA = zeros(1,length(y_a));
for i = 1:length(y_a)
    BA = sqrt(y_a(i)^2+x_a^2);
    psi = asind(x_a/BA);
    theta = acosd((BC^2-BA^2-AC^2)/(-2*BA*AC));
    phi = asind(AC*sind(theta)/BC);
    alpha(i) = psi + phi;
    beta(i) = abs(90-alpha(i));
    y_c = BC*cosd(alpha(i));
    y_d = CD*sind(alpha(i))-y_c;
    x_c = BC*sind(alpha(i));
    
    %Forces
    F_trigger = F_in*y_trigger/x_trigger;
    F = F_trigger/3; %force through each of the 3 slider pinsC_y = -F
    A_x = F*(x_c-x_a)/(y_a(i)-y_c);
    C_x = -A_x;
    C_y = -F;
    F_ac = sqrt(A_x^2+F^2); %axial force through member AC
    D_x = (-C_x*y_c-C_y*x_c)/y_d;
    F_out = D_x*3; %total inward horizontal clamping force
    
    %Mechanical Advantage:
    %Clamp Trigger MA not currently included
    MA(i) = F_out/F_in;
end

%%MA wrt slider height
plot(beta,MA,'b')
hold on

grid on
xlabel('Claw Angle (degrees)')
ylabel('Mechanical Advantage')
title('Mechanical Advantage as a Function of Claw Angle:')
% plot(y_a_max,MA(length(y_a)),'r*')
% plot(y_a_min,MA(1),'r*')
hold off
xlim([0.0 100.0])
ylim([0.0 14.0])

legend('V1)','(V2)')
hold on

%%
%LAD Configuration (V3)
%Claw Dimensions
x_a = 0;
BC = 1.2; %Distance between pin B and C (in) 
AC = 1.4; %Distance between pin A and C (in)
CD = 2.9; %Distance between pin C and D (in)
y_a_min = .7; %Minimum slider height (in)
y_a_max =2.35; %Maximum slider height (in)

%%Static Analysis:
%Geometry:
y_a = (y_a_min:.001:y_a_max); %Slider height range vector
alpha = zeros(1,length(y_a));
beta = zeros(1,length(y_a));
MA = zeros(1,length(y_a));
for i = 1:length(y_a)
    BA = sqrt(y_a(i)^2+x_a^2);
    psi = asind(x_a/BA);
    theta = acosd((BC^2-BA^2-AC^2)/(-2*BA*AC));
    phi = asind(AC*sind(theta)/BC);
    alpha(i) = psi + phi;
    beta(i) = abs(90-alpha(i));
    y_c = BC*cosd(alpha(i));
    y_d = CD*sind(alpha(i))-y_c;
    x_c = BC*sind(alpha(i));
    
    %Forces
    F_trigger = F_in*y_trigger/x_trigger;
    F = F_trigger/3; %force through each of the 3 slider pinsC_y = -F
    A_x = F*(x_c-x_a)/(y_a(i)-y_c);
    C_x = -A_x;
    C_y = -F;
    F_ac = sqrt(A_x^2+F^2); %axial force through member AC
    D_x = (-C_x*y_c-C_y*x_c)/y_d;
    F_out = D_x*3; %total inward horizontal clamping force
    
    %Mechanical Advantage:
    %Clamp Trigger MA not currently included
    MA(i) = F_out/F_in;
end

%%MA wrt slider height
plot(beta,MA,'m')
hold on

grid on
xlabel('Claw Angle (degrees)')
ylabel('Mechanical Advantage')
title('Mechanical Advantage as a Function of Claw Angle:')
% plot(y_a_max,MA(length(y_a)),'k*')
% plot(y_a_min,MA(1),'k*')
legend('Claw Config:(V1)','Claw Config:(V2)','Claw Config:(V3)')
hold on
xlim([0.0 100.0])
ylim([0.0 14.0])