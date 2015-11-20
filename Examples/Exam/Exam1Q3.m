theta = [15;25]*pi/180;
c = [2;1];
alpha = 5*pi/180;
S = sum(c);
delta = [theta(1) - alpha;theta(1) + alpha; theta(2) - theta(1)];

cN = 2*c(1)/S * (sin(delta(1))^2 + sin(delta(2))^2);

cA = 2*c(1)*tan(theta(1))/S * (sin(delta(1))^2 + sin(delta(2))^2) + 4*c(2)*tan(theta(1))/S * sin(delta(3))^2;

cl = cN*cos(alpha) - cA*sin(alpha);
cd = cN*sin(alpha) + cA*cos(alpha);

LD = cl/cd
