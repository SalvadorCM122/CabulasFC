% <strong>RK Solver 2</strong>
% 
% This is the Runge-Kutta Solver for 2 coupled equations.
% 
% <strong>Arguments:</strong>
% 
% <strong>Btable</strong> Butcher Tableau for the RK algorithm to use.
% It's an NxN matrix for a N-1 order RK. 
% 
% <strong>funcs</strong> 1x2 cell array of function handles. It'of the
% form funcs={fx,fy} where fx,fy are the functions of (t,x,v) in a system
%
% dy/dt = fy(t,x,y);          dx/dt = fx(t,x,v)
%
% <strong>conds</strong>  1x5 line vector containing the conditions of the
% problem, in the order: step-h, initial time-t0, final time-tf, initial x
% value x0, initial y value-y0. So conds=[h,t0,tf,x0,y0].
%
% <strong>Outputs:</strong>
%
% <strong>t</strong>,<strong>x</strong>,<strong>y</strong>  line vectors with the values for each variable.

% RKSolver2 by R.Campos (2025)
% ricardo.jpcampos@ua.pt

function [t,x,y] = RKSolver2(Btable,funcs,conds)
    order = size(Btable,2) - 1; 
    C = Btable(:,1);             
    A = Btable(:,2:end);         
    D = Btable(end,2:end);       
    fx = funcs{1};
    fy = funcs{2};
    h = conds(1,1);
    t0 = conds(1,2);
    tf = conds(1,3);
    x0 = conds(1,4);
    y0 = conds(1,5);     
    t = t0:h:tf;                 
    steps = length(t);           
    x = zeros(1, steps);
    y = zeros(1, steps);
    x(1) = x0;
    y(1) = y0;
    R = zeros(order,2);          
    for k = 1:steps-1
        R(1,1) = fx(t(k), x(k), y(k));
        R(1,2) = fy(t(k), x(k), y(k));
        for n = 2:order
            R(n,2) = fy(t(k)+h*C(n,1), x(k)+h*sum(A(n,:).*R(:,1)'), y(k)+h*sum(A(n,:).*R(:,2)'));
            R(n,1) = fx(t(k)+h*C(n,1), x(k)+h*sum(A(n,:).*R(:,1)'), y(k)+h*sum(A(n,:).*R(:,2)'));
        end
        x(1,k+1) = x(1,k)+h*sum(D.*R(:,1)');
        y(1,k+1) = y(1,k)+h*sum(D.*R(:,2)');
    end
end