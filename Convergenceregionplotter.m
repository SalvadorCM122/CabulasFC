% Plotting regions of convergence.

% Create a complex meshgrid
edge = 4;
res = 40;
[X, Y] = meshgrid(linspace(-edge,edge,res), linspace(-edge,edge,res));
% Create the complex variable Z=h\lambda
Z = X + 1i*Y;

% Specify convergence criteria - R(Z)<1  

% Runge-Kutta 4
R4 = 1 + Z + (Z.^2)/2 + (Z.^3)/6 + (Z.^4)/24;

% Runge-Kutta 2
R2 = 1+Z+Z.^2/2;

% Create a contour plot at |R(z)|=1
contour(X, Y, abs(R4), [1 1], 'k', 'LineWidth', 2);
hold on;
contour(X, Y, abs(R2), [1 1], 'k', 'LineWidth', 2);
xlabel('Re(z)'); ylabel('Im(z)');
title('Stability Region of RK[1-4]');
grid on; axis equal;
