% This script is to import the solution grid using a safety factor of 0.3
% and generate figures plotting the paths of the first and last cars for
% four kinds of schemes.
clear;
delimiterIn = ' ';
filename1 = 'C:\Users\Zhichao\OneDrive\PJ-3\f95\13';
s1 = importdata(filename1,delimiterIn);
filename2 = 'C:\Users\Zhichao\OneDrive\PJ-3\f95\23';
s2 = importdata(filename2,delimiterIn);
filename3 = 'C:\Users\Zhichao\OneDrive\PJ-3\f95\33';
s3 = importdata(filename3,delimiterIn);
filename4 = 'C:\Users\Zhichao\OneDrive\PJ-3\f95\43';
s4 = importdata(filename4,delimiterIn);

drawpath(s1)
title('Exact Riemann: first and last vehicles'' paths (Safety Factor of 0.3)');
drawpath(s2)
title('Lax-Wendroff: first and last vehicles'' paths (Safety Factor of 0.3)');
drawpath(s3)
title('Roe''s Approximate: first and last vehicles'' paths (Safety Factor of 0.3)');
drawpath(s4)
title('Entropy Fix: first and last vehicles'' paths (Safety Factor of 0.3)');


