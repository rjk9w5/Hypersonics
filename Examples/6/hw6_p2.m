% 
% 

addpath('/home/ryan/workspace/Hypersonics/HyperOct/Subroutines/')
P = 12044.6;
r = 0.193674;
T = 216.65;
M = 18;
V = M*sqrt(1.4*287*T);
st = state_init('Temperature', T, 'Pressure', P, 'Density', r, 'Mach', M, 'Velocity', V);


[p2p1, T2T1, r2r1, M, ~] = ns(st,'CPG');

p2 = p2p1*P;
T2 = T2T1*T;
r2 = r2r1*r;

TT = 1 + (1.4 - 1) / 2 * M ^ 2;
M;
pt2 = p2*(TT)^(1.4/(1.4-1));
Tt2 = T2*TT;
rt2 = r2*TT^(1/(1.4-1));
estate = state_init('Pressure', pt2, 'Temperature', Tt2, 'Density', rt2);

qdot = qdot_stag(st, estate, 2000, 0.15);

[~, ~, ~, ~, estate_e] = ns(st,'Equillibrium');

TT = 1 + (estate_e.gma - 1) / 2 * estate_e.M ^ 2;
estate_e.M;
estate_e.p = estate_e.p*(TT)^(1.4/(1.4-1));
estate_e.T = estate_e.T*TT;
estate_e.r = estate_e.r*TT^(1/(1.4-1));
[p, s, rho, e, a, h, T] = tgasM(2, estate_e.p, estate_e.p/(287*2000));
wstate = state_init('Pressure', p, 'Temperature', T, 'Density', rho, 'Enthalpy', h, 'Entropy', s, 'Energy', e, 'SoS', a);

qdot = equill_qdot_stag(st, estate_e, wstate, 0.15)