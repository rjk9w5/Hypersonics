clear all
close all
clc
v1 = 4e3; % km/s
alt = 60e3; % km

atmos = dlmread('/home/ryan/Documents/Aero/Data/atmos_data.csv',',',1,0);

index = find(atmos(:,1)==alt);

% Temperature [Kelvin]	Pressure [pascal]	Density [kg/m3]	Speed of sound [m/s]
fs = atmos(index,2:end-1);
T1 = fs(1)
p1 = fs(2)
r1 = fs(3)
a1 = fs(4)
h1 = 1004.5*T1;

% CPG with gamma = 1.4

gmma = 1.4;
M1 = v1 / a1
B = pi/2;
p2p1 = @(M) 1 + 2 * gmma / (gmma + 1) * (M ^ 2 * sin(B)^2 - 1);

T2T1 = @(M) (2 * gmma * M ^ 2 * sin(B)^2 - (gmma - 1))*((gmma - 1) * M ^ 2 * sin(B) ^ 2  + 2 ) ...
	/ ((gmma + 1) ^ 2 * M ^ 2 * sin(B) ^ 2);

h2h1 = T2T1;

r2r1 = @(M) ((gmma + 1) * M ^ 2 * sin(B) ^ 2) / ((gmma - 1) * M ^ 2 * sin(B) ^ 2 + 2);


cpg_p21 = p2p1(M1)
cpg_T21 = T2T1(M1)
cpg_r21 = r2r1(M1)
cpg_h21 = cpg_T21 % for a cpg, cp is constant, therefore the temerpature ratio is the same as the enthalpy ratio 

% euilibrium flow using iterative method and TGAS 
addpath('/home/ryan/Documents/Aero/Hypersonic/Numerical/Eqair/')
eps0 = 0.001;
eps1 = .1;

err = 1e-8;

% Old values
r2 = r1 / eps0;
u2 = eps0 * v1;
p2 = p1 + r1 * v1 * (1 - eps0);
h2 = h1 + v1 ^ 2 / 2 * (1 - eps0 ^ 2);
h2_ = tgas4(p2,r2);

fepso = h2 - h2_;


delta = 1;
it = 0;
figure
hold on

while err < delta
	r2 = r1 / eps1;
	u2 = eps1 * v1;
	p2 = p1 + r1 * v1 ^ 2 * (1 - eps1);
	h2 = h1 + v1 ^ 2 / 2 * (1 - eps1 ^ 2);

	[p, s, rho, e, a, h2_, T] = tgas4(2,p2,r2);

	feps = h2 - h2_;

	tmp = eps1;

	eps1 = eps1 - feps/((feps - fepso)/(eps1-eps0));

	eps0 = tmp;
	fepso = feps;
	delta = abs(feps);
	it = it + 1;
	plot(it,delta)
end
T2 = T;
eq_p21 = p2/p1
eq_T21 = T2/T1
eq_r21 = r2/r1
eq_h21 = h2/h1

rmpath('/home/ryan/Documents/Aero/Hypersonic/Numerical/Eqair/')

pdiff = @(x1,x2) abs(x1 - x2) / x1 * 100;

pd_p21 = pdiff(eq_p21,cpg_p21)
pd_T21 = pdiff(eq_T21,cpg_T21)
pd_r21 = pdiff(eq_r21,cpg_r21)
pd_h21 = pdiff(eq_h21,cpg_h21)
