% Tangent Method 
clear all
close all
clc
addpath('~/Documents/Hypersonics/Subroutines')
warning("off")

st = state_init("Mach", 12, "SpecRatio", 1.4);

c1 = 0.5774;
c2 = -0.2887;

L = 1;
n = 5;
dx = L/n;
x = 0:dx:L;

y = c1.*x + c2.*x.^2;
dy = diff(y);
dydx = dy/dx;
li = atan(dydx);

cd_nm = 0;
cd_tw = 0;
cd_pm = 0;
for i = 1:1:n
    cp_nm(i) = lincl_cp(3, li(i), st);
    prat_guess(i) = cp_nm(i)*st.gma*st.M/2 + 1;
    cd_nm = cd_nm + cp_nm(i)*dy(i)/L;
    
    cp_tw(i) = lincl_cp(9, li(i), st);
    cd_tw = cd_tw + cp_tw(i)*dy(i)/L;
end

cp_nm;




cd_nm = 2*cd_nm
cd_tw = 2*cd_tw

% Shock Wave expansion Method

fcp = @(p21) 2 / (st.gma * st.M ^ 2) * (p21 - 1);

[prat(1),~,~,~,~,B] = obs(st, li(1), "CPG");

Mn = st.M*sin(B);
Mt = st.M*cos(B);
M1n = ((st.gma - 1) * Mn ^ 2 + 2) / (2 * st.gma * Mn ^ 2 - (st.gma - 1));

M(1) = sqrt(M1n ^ 2 + Mt ^ 2);

cp_pm(1) = fcp(prat(1));

dli = -diff(li);

for i = 2:1:n
    nu(i-1) = P_M_angle(M(i-1));
    dnu_dli = @(Mi) nu(i-1) - dli(i-1) - P_M_angle(Mi);

    M(i) = NewtonMethod(prat_guess(i),dnu_dli);
    
    prat(i) = isentropic(M(i), st.gma, "Pressure") / isentropic(M(i-1), st.gma, "Pressure");
    
    pp = prod(prat);
    
    cp_pm(i) = fcp(pp);
end

for i = 1:1:n
    cd_pm = cd_pm + cp_pm(i)*dy(i)/L;
end

cd_pm = 2*real(cd_pm)

pdiff = @(a, g) abs(a - g) / a * 100;

NM = pdiff(real(cd_pm), real(cd_nm))
TW = pdiff(real(cd_pm), real(cd_tw))

figure
subplot(2,1,1)
hold on
stairs(x,[cp_nm,cp_nm(end)],'r')
stairs(x, [cp_tw, cp_tw(end)],'b')
stairs(x, [cp_pm, cp_pm(end)],'g')
legend("Newtonian", "Tangent-Wedge", "Prandtl-Meyer")
title("C_p distribution")
ylabel("C_p")

subplot(2,1,2)
plot(x,y,'k')
xlabel("x")
ylabel("y")
