%% Método de Shooting

clc, clear all, close all

%Constantes
T = 5.0E4 ; w = 1.0E5 ; a = 5.0E-8; L=3; h=0.01;

%Arrays
x = 0:h:L ; N = length(x) ; y = zeros(1,N) ; y(1) = 0 ; y(L) = 0;

%Relacionado a Shoooting
guess(1) = -2; guess(2) = -2.2 ; B=0 ; tol = 1E-4;

for is = 1:150
    dy = zeros(1,N);
    dy(1) = guess(is);

    for k=1:N-1
        dy(k+1) = dy(k) + h*(2*a*T*y(k) + a*w*x(k)*(L-x(k)));
        y(k+1) = y(k) + h*dy(k);
    end

    result(is) = y(end);

    if (is>1)
        m = ( result(is) - result(is-1) ) / (guess(is) - guess(is-1));
        guess(is+1) = guess(is) +  (B - result(is)) / m;

        %Critério de Paragem

        if abs(B - result(is)) < tol
            fprintf('Numero de iterações %d',is)   
            break
        end
    end
end

plot(x,y)
xlabel("x (m)")
ylabel("y (m)")

%% Shooting com ODE45
clc; close all; clear all;

%constantes
m=1.5; K=2; alpha=-0.2; h=0.001;
t0=0; t_end=5;
x0=1.9; v0=0; 

B=-1.5; tol=1E-4;

guess(1)=-0.2;
guess(2)=-0.18;


for iw=1:150 %número máximo de tentativas
    alpha=guess(iw); %cada tentativa atualiza a frequência para tentar acertar

    % Iteração para obter os valores (método de Euler-Cromer)
    options=odeset('RelTol',3E-14,'AbsTol',[1E-13 1E-13]);
    [t, sol] = ode45(@(t, sol) func(t, sol, m, K,alpha), t0:h:t_end, [x0 v0], options);
    
    x=sol(:,1);v=sol(:,2);

    current_min = min(x); % Pegamos o valor mínimo (amplitude negativa)
    result(iw) = current_min;
     
    % Critério de parada
    if abs(current_min - B) < tol
        fprintf('Convergido! alpha = %.6f produz amplitude mínima de %.6f\n', alpha, current_min);
        break;
    end
    
    % Método da secante (após ter pelo menos dois pontos)
    if iw > 1
        m_sec = (result(iw) - result(iw-1)) / (guess(iw) - guess(iw-1));
        guess(iw+1) = guess(iw) + (B - result(iw)) / m_sec;
        
        % Prevenção contra divisão por zero ou valores absurdos
        if abs(m_sec) < 1e-10 || isnan(guess(iw+1)) || isinf(guess(iw+1))
            guess(iw+1) = guess(iw) * 0.95; % Pequeno ajuste se o método da secante falhar
        end
    end
end

function derivadas=func(t, sol, m, K,alpha)
        derivadas=zeros(2,1);
        x=sol(1); y=sol(2);
        derivadas(1)=y; %MUDAR A EXPRESSÃO DA DERIVADA - dx/dt=v - primeira derivada ou função
        derivadas(2)=-K*x/m*(1+3*alpha*x/2); %MUDAR A EXPRESSÃO DA DERIVADA - dv/dt - segunda derivada ou função
end

%% Shooting com 2 Eulers

clc; close all; clear all

%constantes
g=9.8; h=0.1; massa=0.5; D=0.2;
alcance=50; alfa=deg2rad(10);

%arrays
t=0:h:20; Nt=length(t);
x=zeros(1,Nt); 
y=zeros(1,Nt); y(1)=20;
vx=zeros(1,Nt);
vy=zeros(1,Nt);

B=50; tol=1E-4;
guess(1)=39; guess(2)=39.9;

for is = 1:150
    v=guess(is);
    vx(1)=v*cos(alfa);
    vy(1)=v*sin(alfa);

    for k=1:Nt-1
        vx(k+1) = vx(k) + h*(-D*vx(k)/massa);
        x(k+1) = x(k) + h*vx(k);

        vy(k+1) = vy(k) + h*(-g-D*vy(k)/massa);
        y(k+1) = y(k) + h*vy(k);

        if y(k+1)<=0
            idx=k+1;
            x(k+1:Nt) = x(k+1);
            y(k+1:Nt) = 0;
            break
        end
    end
    
    x_int=interp1(y(idx-1:idx),x(idx-1:idx),0);
    
    result(is) = x_int;

    if (is>1)
        m = ( result(is) - result(is-1) ) / (guess(is) - guess(is-1));
        guess(is+1) = guess(is) +  (B - result(is)) / m;

        %Critério de Paragem

        if abs(B - result(is)) < tol
            fprintf('Numero de iterações %d',is)   
            break
        end
    end
end

plot(x,y)
xlabel('x(m)');
ylabel('y(m)');

%% Método de Jacobi

clc; clear all; close all;

% Definição da matriz A
A = -eye(N,N);

% Definição do vetor B
B = zeros(N,1);

% Inicialização do método de Jacobi

x_old = zeros(N,1);
x_new = zeros(N,1);
max_iter = 1000;
tol = 10E-7;

for k = 1:max_iter
    for i = 1:N
        sigma = 0;
        for j = 1:N
            if j ~= i
                sigma = sigma + A(i,j) * x_old(j);
            end
        end
        x_new(i) = (B(i) - sigma) / A(i,i);
    end

    % Critério de paragem
    if norm(x_new - x_old, inf) < tol
        fprintf('Convergiu em %d iterações.\n', k);
        break;
    end
    x_old = x_new;
end

% Solução com linsolve
x_linsolve = linsolve(A, B);

% Exibir resultados
disp('Solução aproximada com Jacobi:');
disp(x_new);

disp('Solução exata com linsolve:');
disp(x_linsolve);

%% Diferenças finitas centradas

f' = (T(k+1) - T(k-1)) / 2h
f'' = (T(k-1) - 2*T(k) + T(k+1)) / h^2

%% Diferenças finitas avançadas

f' = T(k+1) - T(k) / h
f'' = (T(k-1) - 2*T(k) + T(k+1)) / h^2
