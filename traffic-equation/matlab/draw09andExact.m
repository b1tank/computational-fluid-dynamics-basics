% this script is used to import the solution grid using a safety factor of
% 0.9 and generate figures plotting the solution at t=0s, t=20s, t=40s and
% t=60s for four kinds of schemes and the exact solution.
%
% Besides, the speeds with which the vehicles pass the light using a safety 
% factor of 0.9 are also plotted in a separate figure for four kinds of 
% schemes and the exact solution.
clear;
delimiterIn = ' ';
filename1 = 'C:\Users\Zhichao\OneDrive\PJ-3\f95\19';
s1 = importdata(filename1,delimiterIn);
filename2 = 'C:\Users\Zhichao\OneDrive\PJ-3\f95\29';
s2 = importdata(filename2,delimiterIn);
filename3 = 'C:\Users\Zhichao\OneDrive\PJ-3\f95\39';
s3 = importdata(filename3,delimiterIn);
filename4 = 'C:\Users\Zhichao\OneDrive\PJ-3\f95\49';
s4 = importdata(filename4,delimiterIn);
filename5 = 'C:\Users\Zhichao\OneDrive\PJ-3\f95\59';
s5 = importdata(filename5,delimiterIn);

figure(1) % plotting the solution using a safety factor of 0.9 and Exacte Riemann method
hold on
plot(s1(1,2:152),s1(2,2:152),s1(1,2:152),s1(24,2:152),s1(1,2:152),s1(46,2:152),s1(1,2:152),s1(68,2:152));
title('Exact Riemann, safety factor of 0.9');
xlabel('x (ft)');
ylabel('Density of cars (vehicles/ft)');
legend('t = 0 s','t = 20 s','t = 40 s','t = 60 s');
xlim([-2000,7000]);
ylim([0,0.17]);
hold off

figure(2) % plotting the solution using a safety factor of 0.9 and Lax-Wendroff scheme
hold on
plot(s2(1,2:152),s2(2,2:152),s2(1,2:152),s2(24,2:152),s2(1,2:152),s2(46,2:152),s2(1,2:152),s2(68,2:152));
title('Lax-Wendroff, safety factor of 0.9');
xlabel('x (ft)');
ylabel('Density of cars (vehicles/ft)');
legend('t = 0 s','t = 20 s','t = 40 s','t = 60 s');
xlim([-2000,7000]);
ylim([0,0.17]);
hold off

figure(3) % plotting the solution using a safety factor of 0.9 and Roe's Approximate Riemann Solver
hold on
plot(s3(1,2:152),s3(2,2:152),s3(1,2:152),s3(24,2:152),s3(1,2:152),s3(46,2:152),s3(1,2:152),s3(68,2:152));
title('Roe''s Approximate Riemann Solver, safety factor of 0.9');
xlabel('x (ft)');
ylabel('Density of cars (vehicles/ft)');
legend('t = 0 s','t = 20 s','t = 40 s','t = 60 s');
xlim([-2000,7000]);
ylim([0,0.17]);
hold off

figure(4) % plotting the solution using a safety factor of 0.9 and Roe's Approximate Riemann Solver with Entropy Fix
hold on
plot(s4(1,2:152),s4(2,2:152),s4(1,2:152),s4(24,2:152),s4(1,2:152),s4(46,2:152),s4(1,2:152),s4(68,2:152));
title('Roe''s Approximate Riemann Solver with Entropy Fix, safety factor of 0.9');
xlabel('x (ft)');
ylabel('Density of cars (vehicles/ft)');
legend('t = 0 s','t = 20 s','t = 40 s','t = 60 s');
xlim([-2000,7000]);
ylim([0,0.17]);
hold off

figure(5) % plotting the exact solution
hold on
plot(s5(1,2:152),s5(2,2:152),s5(1,2:152),s5(24,2:152),s5(1,2:152),s5(46,2:152),s5(1,2:152),s5(68,2:152));
title('Exact Solution');
xlabel('x (ft)');
ylabel('Density of cars (vehicles/ft)');
legend('t = 0 s','t = 20 s','t = 40 s','t = 60 s');
xlim([-2000,7000]);
ylim([0,0.17]);
hold off

v_m=50;
rho_m=0.1;

figure(6) % plotting the speeds with which the vehicles pass the light where x = 0.
hold on
plot(s1(2:68,1),v_m.*(1.-s1(2:68,22)./rho_m),s2(2:68,1),v_m.*(1.-s2(2:68,22)./rho_m),s3(2:68,1),v_m.*(1.-s3(2:68,22)./rho_m),...
    s4(2:68,1),v_m.*(1.-s4(2:68,22)./rho_m),s5(2:68,1),v_m.*(1.-s5(2:68,22)./rho_m));
title('Safety factor of 0.9');
xlabel('t (s)');
ylabel('Speed with which the vehicles pass the light (ft/s)');
xlim([0 60]);
ylim([0 30]);
legend('Exact Riemann','Lax-Wendroff','Roe''s Approximate Riemann Solver','Roe''s Approximate Riemann Solver with Entropy Fix','Exact Solution');
hold off

