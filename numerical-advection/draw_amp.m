filename1 = 'A:\$_Courses\!_Computational_Fluid_Dynamics\PJ-2\f95\a';

delimiterIn = ' ';
a_n0 = importdata(filename1,delimiterIn);
a_n = a_n0.data;

a_t = [1,1,1,1,1;0.75,0.991639,0.996391,0.999955,0.997183;0.5,0.985624,0.995185,0.999965,0.995185;0.25,0.981997,0.996391,0.999989,0.993992];

figure
hold on
plot(a_n(:,1),a_n(:,2:5),'-*');
plot(a_t(:,1),a_t(:,2:5),'-o');
title('Comparison of amplitudes');
xlabel('\nu');
ylabel('Amplitude');
legend('numerical: 1','numerical: |\nu|','numerical: \nu^2','numerical: 1/3+2/3*\nu^2','theoretical: 1','theoretical: |\nu|','theoretical: \nu^2','theoretical: 1/3+2/3*\nu^2');
hold off