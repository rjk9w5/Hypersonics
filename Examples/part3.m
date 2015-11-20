clear all
close all
clc
v1 = 6069; % m/s
alt = 30480; % km

atmos = dlmread('/home/ryan/Documents/Aero/Data/atmos_data.csv',',',1,0);
addpath('/home/ryan/Documents/Aero/Hypersonic/Numerical/')

index1 = max(find(atmos(:,1)<alt));
index2 = min(find(atmos(:,1)>alt));
m = atmos(index1,2:end-1) - atmos(index2,2:end-1) ./ (atmos(index1,1) - atmos(index2,1));

fs = m .* (alt - atmos(index1,1)) + atmos(index1,2:end-1);


% Temperature [Kelvin]	Pressure [pascal]	Density [kg/m3]	Speed of sound [m/s]
fs = atmos(index1,2:end-1);
T1 = fs(1)
p1 = fs(2)
r1 = fs(3)
a1 = fs(4)
h1 = 1004.5*T1;

% CPG with gamma = 1.4

gmma = 1.4;
M1 = v1 / a1
delta = 20*pi/180;
B = T_B_M(delta,M1,gmma);
p2p1 = @(M) 1 + 2 * gmma / (gmma + 1) * (M ^ 2 * sin(B) ^ 2 - 1);

T2T1 = @(M) (2 * gmma * M ^ 2 * sin(B) ^ 2 - (gmma - 1))*((gmma - 1) * M ^ 2 * sin(B) ^ 2  + 2 ) ...
	/ ((gmma + 1) ^ 2 * M ^ 2 * sin(B) ^ 2);

r2r1 = @(M) ((gmma + 1) * M ^ 2 * sin(B) ^ 2) / ((gmma - 1) * M ^ 2 * sin(B) ^ 2 + 2);


cpg_p21 = p2p1(M1)
cpg_T21 = T2T1(M1)
cpg_r21 = r2r1(M1)
cpg_h21 = cpg_T21;
cpg_B_T = (B - delta)*180/pi

% euilibrium flow using iterative method and TGAS 
addpath('/home/ryan/Documents/Aero/Hypersonic/Numerical/Eqair/')
eps0 = .001;
eps1 = .1;

err = 1e-7;

% Old values
r2 = r1 / eps0;
ee = 1-eps0;
B = atan((ee - sqrt(ee ^ 2 - 4 * eps0 * tan(delta) ^ 2))/(2 * eps0 * tan(delta)));
u1n = v1*sin(B);
u2n = eps0 * u1n;
p2 = p1 + r1 * u1n ^ 2 * (1 - eps0);
h2 = h1 + u1n ^ 2 / 2 * (1 - eps0 ^ 2);
h2_ = tgasM(p2,r2);

fepso = h2 - h2_;


del = 1;
it = 0;
figure
hold on

while err < del
	r2 = r1 / eps1;
	ee = 1 - eps1;
	B = atan((ee - sqrt(ee ^ 2 - 4 * eps1 * tan(delta) ^ 2))/(2 * eps1 * tan(delta)));
	u1n = v1 * sin(B);
	u2n = eps1 * u1n;
	p2 = p1 + r1 * u1n ^ 2 * (1 - eps1);
	h2 = h1 + u1n ^ 2 / 2 * (1 - eps1 ^ 2);

	[p, s, rho, e, a, h2_, T] = tgasM(2,p2,r2);

	feps = h2 - h2_;

	tmp = eps1;

	eps1 = eps1 - feps/((feps - fepso)/(eps1-eps0));

	eps0 = tmp;
	fepso = feps;
	del = abs(feps);
	it = it + 1;
	plot(it,delta)
end
T2 = tgas3(p2,r2);
eq_p21 = p2/p1
eq_T21 = T2/T1
eq_r21 = r2/r1
eq_h21 = h2/h1;
eq_B_T = (B - delta)*180/pi

rmpath('/home/ryan/Documents/Aero/Hypersonic/Numerical/Eqair/')
rmpath('/home/ryan/Documents/Aero/Hypersonic/Numerical/')

pdiff = @(x1,x2) abs(x1 - x2) / x1 * 100;

pd_p21 = pdiff(eq_p21,cpg_p21)
pd_T21 = pdiff(eq_T21,cpg_T21)
pd_r21 = pdiff(eq_r21,cpg_r21)
pd_b_t21 = pdiff(eq_B_T,cpg_B_T)
