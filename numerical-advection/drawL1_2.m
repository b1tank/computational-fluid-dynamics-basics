filename1 = 'A:\$_Courses\!_Computational_Fluid_Dynamics\PJ-2\f95\L';

delimiterIn = ' ';
L1_n0 = importdata(filename1,delimiterIn);
L1_n = L1_n0.data;

figure
hold on
plot(L1_n(:,1),L1_n(:,2:5),'-*');

% set(gca,'yscale','log');
title('L_1 norms');
xlabel('\nu');
ylabel('L_1 norms');
legend('1','|\nu|','\nu^2','1/3+2/3*\nu^2');
hold off