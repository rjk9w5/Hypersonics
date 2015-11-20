clear all
clc
close all
v1 = 3000; % m/s
T1 = 300; % K
tv = 3390;

M = 14.0067*2;
Ru = 8314;
Rn2 = Ru/M;

fh= @(T,R) 7/2*R*T + R * T * tv/T/(exp(tv./T)  -1);

h1 = fh(T1,Rn2);

addpath('/home/ryan/Documents/Aero/Hypersonic/Numerical/Eqair/')
eps0 = 0.001;
eps1 = 0.1;

err = 1e-8;

% Old values
u2 = eps0 * v1;
T2 = eps0*(T1 + v1 ^ 2 / Rn2 * (1 - eps0));
h2 = h1 + v1 ^ 2 / 2 * (1 - eps0 ^ 2);

h2_ = fh(T2,Rn2);

fepso = h2 - h2_;


delta = 1;
it = 0;
figure
hold on

while err < delta
	u2 = eps1 * v1;
	T2 = eps1*(T1 + v1 ^ 2 / Rn2 * (1 - eps1));
	h2 = h1 + v1 ^ 2 / 2 * (1 - eps1 ^ 2);

	h2_ = fh(T2,Rn2);

	feps = h2 - h2_;

	tmp = eps1;

	eps1 = eps1 - feps/((feps - fepso)/(eps1-eps0));

	eps0 = tmp;
	fepso = feps;
	delta = abs(feps);
	it = it + 1;
	plot(it,delta)
end
T21 = T2/T1

rmpath('/home/ryan/Documents/Aero/Hypersonic/Numerical/Eqair/')


% temperature ratio over normal shock