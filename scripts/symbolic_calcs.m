clc
clear
syms theta0 A B theta theta1 theta2 beta real 
assume(beta , 'positive')
beam = exp(-beta*theta)*(A*sin(beta*theta) +B*cos(beta*theta));
normal_force_int = int(beam*cos(theta) , theta);
tangent_force_int = int(beam*sin(theta), theta);
disp("normal force:")
pretty(normal_force_int)
disp("tangential force")
pretty(tangent_force_int)