clear; clc; close all
%%
%Static pin analysis (Cantilever beam model)
S_yield = 42351.02; %estimated pin yield strength (psi, lb/in^2)
d = 0.065; %pin shaft diam (in)
I = pi*d^4/64; %Area Moment of Inertia (in^4)
y_cent = d/2; %Centroid dist (in)
FS = 4;  %Factor of Safety
S_allow = S_yield/FS; %Allowable stress (psi, lb/in^2)

M_max = S_allow*I/y_cent; %max moment in pin at surface plate (lb/in)

L_min = .05; %minimum pin extension based on CAD (in)
L_max = 0.75; %maximum pin extension based on CAD (in)

%%
%Plotting
L= linspace(L_min,L_max,100); 
F_max = zeros(1,100);
for i = 1:100
    F_max(i) = M_max/L(i);
end
plot(L,F_max,'r') 
title("Maximum Allowable Force On One Pin")
ylabel("Force (lb)")
xlabel("Pin Extension (in)")
grid on