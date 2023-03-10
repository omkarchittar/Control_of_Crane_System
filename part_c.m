%requires Symbolic Math Toolbox
clear all
syms m1 g m2 M L1 L2
A = [0 1 0 0 0 0; 0 0 -m1*g/M 0 -m2*g/M 0; 0 0 0 1 0 0; 0 0 -((M*g)+(m1*g))/(M*L1) 0 -g*m2/(M*L1) 0; 0 0 0 0 0 1; 0 0 -m1*g/(M*L2) 0 -((M*g)+(m2*g))/(M*L2) 0];
B = [0; 1/M; 0; 1/(L1*M); 0; 1/(L2*M)];
Rank = rank([B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B]);

%value of D denotes the determinant and same is shared in the report.
D = det([B A*B (A^2)*B (A^3)*B (A^4)*B (A^5)*B]); 