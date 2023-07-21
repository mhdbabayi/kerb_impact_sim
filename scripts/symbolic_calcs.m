clc
clear
syms theta0 A B theta theta1 theta2 beta real 
assume(beta , 'positive')
beam = exp(-beta*theta)*(A*sin(beta*theta) +B*cos(beta*theta))
normal_force_int = int(beam*cos(theta) , theta);
tangent_force_int = int(beam*sin(theta), theta);
disp("normal force:")
pretty(normal_force_int)
disp("tangential force")
pretty(tangent_force_int)
%%
clear all
clc
trapz([1:10])
tyre_radius = 0.788/2;
terrain_radius = 5;
penetration = 0.08;
d = tyre_radius + terrain_radius - penetration;
r = @(theta) (d * cos(theta) - terrain_radius * sqrt(1 - d^2*sin(theta).^2/(terrain_radius^2)));
t = deg2rad(linspace(0 , 90 , 90));
all_r = r(t);
sep_idx= find(all_r > tyre_radius , 1 , 'first');
tic
trapz(t(1:sep_idx) , all_r(1:sep_idx));
toc
tic
b  = 10^3 - 5^3;
toc