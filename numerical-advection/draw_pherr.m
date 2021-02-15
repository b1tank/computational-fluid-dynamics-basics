filename1 = 'A:\$_Courses\!_Computational_Fluid_Dynamics\PJ-2\f95\p';

delimiterIn = ' ';
p_n0 = importdata(filename1,delimiterIn);
p_n = p_n0.data;

p_t = [1,0,0,0,0;0.75,0.005636,0.000805,-0.002787,4.51326e-6;0.5,0.009732,0,-0.004796,0;0.25,0.01222,-0.002415,-0.006008,-5.82941e-6];

figure
hold on
plot(p_n(:,1),p_n(:,2:5),'-*');
plot(p_t(:,1),p_t(:,2:5),'-o');
% set(gca,'yscale','log');
title('Comparison of phase errors');
xlabel('\nu');
ylabel('Phase error');
legend('numerical: 1','numerical: |\nu|','numerical: \nu^2','numerical: 1/3+2/3*\nu^2','theoretical: 1','theoretical: |\nu|','theoretical: \nu^2','theoretical: 1/3+2/3*\nu^2');
hold off