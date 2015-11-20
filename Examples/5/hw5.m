% Homework 5
close all
clear all
clc
warning("off")
addpath('~/Documents/Hypersonics/Subroutines')
state = state_init();
state.T = 216.650;     % temperature [K]
state.h = state.cp*state.T;
state.p = 12044.6;     % pressure [Pa]
state.r = 0.193674;     % density [kg/m^3]
state.a = 295.070;
state.M = 18;     % Mach number
state.V = state.M*state.a;     % velocity magnitude [m/s]

% Parameters
Rn = 0.15;
thc1 = 30*pi/180;
thc2 = 5*pi/180;
c1 = 0.5;
c2 = 0.5;

% chordwise lengths
x12 = Rn*(1-sin(thc1));
x23 = c1 - x12;
x34 = c2 + c1;

% curve-wise lengths
s12 = (pi/2-thc1)*Rn;
s23 = s12 + x23/cos(thc1);
s34 = s23 + c2/cos(thc2);

smsh = linspace(0,s34,1000);
[~,it1] = min(abs(smsh-s12));
[~,it2] = min(abs(smsh-s23));
[~,it3] = min(abs(smsh-s34));

for i = 1:1:it1
    nm_cp(i) = lincl_cp(3, (pi/2-smsh(i)/Rn), state);
    mnmCPG_cp(i) = lincl_cp(4, (pi/2-smsh(i)/Rn), state);
    mnmEquil_cp(i) = lincl_cp(5, (pi/2-smsh(i)/Rn), state);
end

nm_cp(it1+1:it2) = lincl_cp(3, thc1);
mnmCPG_cp(it1+1:it2) = lincl_cp(4, thc1, state);
mnmEquil_cp(it1+1:it2) = lincl_cp(5, thc1, state);

mnmCPG_cp(it2+1:it3) = lincl_cp(4, thc2, state);
mnmEquil_cp(it2+1:it3) = lincl_cp(5, thc2, state);
nm_cp(it2+1:it3) = lincl_cp(3, thc2);

figure
plot(smsh/Rn,nm_cp)
title("Newtonian Method")

figure
plot(smsh/Rn,mnmEquil_cp)
title("Modified Newtonian Method - Equillibrium")

figure
plot(smsh/Rn,mnmCPG_cp)
title("Modified Newtonian Method - CPG")

% Make a profile plot
xmsh = linspace(0,1,1000);
[~,it1] = min(abs(xmsh-x12));
[~,it2] = min(abs(xmsh-c1));
[~,it3] = min(abs(xmsh-(c1+c2)));

prfu(1:it1) = sqrt(Rn^2 - (xmsh(1:it1)-Rn).^2);
prfu(it1+1:it2) = (xmsh(it1+1:it2)-x12)*tan(thc1)+prfu(it1);
prfu(it2+1:it3) = (xmsh(it2+1:it3)-c2)*tan(thc2) + prfu(it2);
prfl = -prfu;
cl(1:11) = 0;
figure
hold on
plot(xmsh,prfu,'k')
plot(xmsh,prfl,'k')
plot(xmsh(1:10:end),zeros(1,100),'-.k')
title("Geometry Profile")
xlabel("Chord Scale")

pause(6)
