function Op = propper2(t2, t1, gamma, omega, alpha, beta, omega0);
% gamma = 4;
% omega = 2;
% omega0 = 2;
% alpha = 2;
% beta = 1;

Rmat = @(theta) [cos(theta), sin(theta); -sin(theta), cos(theta)];
Bmat = @(t) [alpha*t, gamma+beta*t; -beta*t, alpha*t];
fun = @(t,y) Rmat(omega*t)*Bmat(t)*Rmat(omega*t)'*y; 
vopt = odeset ("RelTol", 1e-12, "AbsTol", 1e-12, 'InitialStep', 1e-2, 'MaxStep',1e-2);

[t, y1] = ode45(fun, [t1, t2], [1; 0], vopt);
[t, y2] = ode45(fun, [t1, t2], [0; 1], vopt);
Op = [y1(end,:)', y2(end,:)'];
Op = Rmat(omega0*t2)'*Op*Rmat(omega0*t1);
