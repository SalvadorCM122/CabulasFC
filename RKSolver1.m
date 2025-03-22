function [t,x] = RKSolver1(Btable,func,conds)
    order = size(Btable,2) - 1; 
    C = Btable(:,1);             
    A = Btable(:,2:end);         
    D = Btable(end,2:end);       
    f = func;                    
    x0 = conds(1,1);             
    t0 = conds(1,2);            
    tf = conds(1,3);            
    h = conds(1,4);              
    t = t0:h:tf;                 
    steps = length(t);           
    x = zeros(1, steps);         
    x(1) = x0;                   
    R = zeros(order,1);          

    for k = 1:steps-1
        R(1,1) = f(t(k), x(k));
        for n = 2:order
            R(n,1) = f(t(k) + h * C(n,1), x(k) + h * sum(A(n,:) .* R'));
        end
        x(1, k+1) = x(1, k) + h * sum(D .* R');
    end
end
