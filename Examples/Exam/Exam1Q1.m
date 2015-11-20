% This is the driver to solve problem 1

% Initialize state variable for freestream conditions
init_state = state_init()
init_state.T = 237;
init_state.p = 559;
init_state.M = 25;
init_state.a = sqrt(1.4*286.9*init_state.T);
init_state.V = init_state.M * init_state.a;
init_state.h = 1004.5*init_state.T;
init_state.r = init_state.p/(286.9*init_state.T);

% Geometry Variables
theta = [15; 5; 25] * pi/180;
alpha = 10 * pi/180;
c = [1.25;0.75];
S = sum(c)
q_inf = .5*init_state.r*init_state.v^2;

% Flow Deflection angles
delta = [theta(1)-alpha; theta(2)+alpha; theta(3)-theta(2)];

% Solve surfaces with CPG assumptions
[cpg_face(1), cpg_bd(1)] = cpg_obs(init_state, delta(1), 1.4);
[cpg_face(2), cpg_bd(2)] = cpg_obs(init_state, delta(2), 1.4);
[cpg_face(3), cpg_bd(3)] = cpg_obs(cpg_face(2), delta(3), 1.4);

% Solve upper surface with Equillibrium Assumption
[equil_face(1), equil_bd(1)] = equil_obs(init_state, delta(1));
[equil_face(2), equil_bd(2)] = equil_obs(init_state, delta(2));
[equil_face(3), equil_bd(3)] = equil_obs(equil_face(2), delta(3));


fpdiff = @(val1, val2) abs(val1 - val2)/val1 * 100;

pdiffs = zeros(3,4);

for f = 1:3
% Calculate % Difference between CPG and Equillibrium Assumptions
	pdiffs(f,2) = fpdiff(equil_face(f).p,cpg_face(f).p);
	pdiffs(f,4) = fpdiff(equil_face(f).r,cpg_face(f).r);
	pdiffs(f,3) = fpdiff(equil_face(f).T,cpg_face(f).T);
	pdiffs(f,1) = fpdiff(equil_bd(f),cpg_bd(f));
	
% Calculate the compressibility and molecular weight of the mixture on each face
	[z(f), M(f)] = postprocessor_d(equil_face(f));
	
% Calculate the cp values for each face
	if f < 3
		cp_HSP(f) = lincl_cp(2,delta(f),init_state,1.4);
	elseif f == 3
		% The lcp function will return the cp wrt face 2. It is desired to have
		% it wrt to the freestream conditions. Because of the nature the the
		% HSP, the results from the cpg assumption are used to get the Mach
		% number over face 2
		cp_HSP(3) = lincl_cp(2,delta(3),cpg_face(2),1.4);
		
		% Calculate the dynamic pressure for face 2 (freestream to face 3)
		q1 = .5*cpg_face(2).r*cpg_face(2).v^2;
		p2 = cp_HSP(2)*q_inf + init_state.p;
		
		% Back out the value of the pressure on face 3 from the previously 
		% calculated value of cp using HSP
		tmp = cp_HSP(f)*q1 + p2;
		
		% Re-calculate the cp for HSP wrt to the freestream for the vehicle
		cp_HSP(3) = (tmp - init_state.p)/q_inf;
	end
	cp_equil(f) = (equil_face(f).p - init_state.p)/q_inf;
	cp_cpg(f) = (cpg_face(f).p - init_state.p)/q_inf;	
end

pdiffs

cp_equil
cp_cpg
cp_HSP

fcl = @(cp) cos(alpha) * (-cp(1) + cp(2)*c(1)/S + ...
	cp(3)*c(2)/S) - sin(alpha)*(cp(1)*tan(theta(1)) + ...
	cp(2)*tan(theta(2))*c(1)/S + cp(3)*tan(theta(3))*c(2)/S)
    
fcd = @(cp) sin(alpha) * (-cp(1) + cp(2)*c(1)/S + ...
	cp(3)*c(2)/S) + cos(alpha)*(cp(1)*tan(theta(1)) + ...
	cp(2)*tan(theta(2))*c(1)/S + cp(3)*tan(theta(3))*c(2)/S)

% Lift and Drag coefficients with Equillibrium Assumption
cl_equil = cos(alpha) * (-cp_equil(1) + cp_equil(2)*c(1)/S + ...
	cp_equil(3)*c(2)/S) - sin(alpha)*(cp_equil(1)*tan(theta(1)) + ...
	cp_equil(2)*tan(theta(2))*c(1)/S + cp_equil(3)*tan(theta(3))*c(2)/S)

cd_equil = sin(alpha) * (-cp_equil(1) + cp_equil(2)*c(1)/S + ...
	cp_equil(3)*c(2)/S) + cos(alpha)*(cp_equil(1)*tan(theta(1)) + ...
	cp_equil(2)*tan(theta(2))*c(1)/S + cp_equil(3)*tan(theta(3))*c(2)/S)
	
ld_equil = cl_equil/cd_equil
% Lift and Drag coefficients with Equillibrium Assumption
cl_cpg = cos(alpha) * (-cp_cpg(1) + cp_cpg(2)*c(1)/S + ...
	cp_cpg(3)*c(2)/S) - sin(alpha)*(cp_cpg(1)*tan(theta(1)) + ...
	cp_cpg(2)*tan(theta(2))*c(1)/S + cp_cpg(3)*tan(theta(3))*c(2)/S)

cd_cpg = sin(alpha) * (-cp_cpg(1) + cp_cpg(2)*c(1)/S + ...
	cp_cpg(3)*c(2)/S) + cos(alpha)*(cp_cpg(1)*tan(theta(1)) + ...
	cp_cpg(2)*tan(theta(2))*c(1)/S + cp_cpg(3)*tan(theta(3))*c(2)/S)

ld_cpg = cl_cpg/cd_cpg
% Lift and Drag coefficients with Equillibrium Assumption
cl_HSP = cos(alpha) * (-cp_HSP(1) + cp_HSP(2)*c(1)/S + ...
	cp_HSP(3)*c(2)/S) - sin(alpha)*(cp_HSP(1)*tan(theta(1)) + ...
	cp_HSP(2)*tan(theta(2))*c(1)/S + cp_HSP(3)*tan(theta(3))*c(2)/S)

cd_HSP = sin(alpha) * (-cp_HSP(1) + cp_HSP(2)*c(1)/S + ...
	cp_HSP(3)*c(2)/S) + cos(alpha)*(cp_HSP(1)*tan(theta(1)) + ...
	cp_HSP(2)*tan(theta(2))*c(1)/S + cp_HSP(3)*tan(theta(3))*c(2)/S)

ld_HSP = cl_HSP/cd_HSP

cl_diff(1) = fpdiff(cl_equil,cl_cpg);
cl_diff(2) = fpdiff(cl_equil,cl_HSP);

cd_diff(1) = fpdiff(cd_equil,cd_cpg);
cd_diff(2) = fpdiff(cd_equil,cd_HSP);

ld_diff(1) = fpdiff(ld_equil,ld_cpg);
ld_diff(2) = fpdiff(ld_equil,ld_HSP);

cl_diff
cd_diff
ld_diff
