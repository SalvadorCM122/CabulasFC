
% ESTE É O FICHEIRO DE DESENVOLVIMENTO.

%NOTAS DE COMO USAR ESTE FICHEIRO:

%----------------------------Valores a se mudar--------------------------%

% Para cada método basta mudar os valores iniciais, as constantes e a 
% função f. A função f é a função que se encontra na 
% equação: dx/dt = f(x,t)

%----------------------------Métodos implementados-----------------------%

%  Método de Euler (Primeira e Segunda Ordem)
%  Método de Euler-Cromer (Segunda Ordem)
%  Método de Crank-Nicolson (Primeira e Segunda Ordem)
%  Método de Range-Kutta (Segunda e Quarta Ordem)

%------------------------------------------------------------------------%

%% Método de Euler - Primeira Ordem

clc, clear all, close all

% Condições iniciais & finais
t0 = 0 ; t_end = 10;
x0 = 1;
h = 0.01;

%Constantes
a = 1 ; b = 2;

%Definir a função f(t, x)
f = @(t, x) x+t; 

% Método de Euler - 1ª ordem
[t, x] = euler_method(f, t0, x0, h, t_end);

% Plot da solução
plot(t, x, '-');
xlabel('t');
ylabel('x');
title('Solução pelo método de Euler');
grid on;


function [t, x] = euler_method(f, t0, x0, h, t_end)

    % Criar os arrays para armazenar os dados
    t = t0:h:t_end; % Criar o array do tempo
    x = zeros(1,length(t)); % Inicializar variável
    x(1) = x0; % Condicção inicial da variável

    % Iteração método de Euler   
    for n = 1:length(t)-1
        x(n+1) = x(n) + h * f(t(n), x(n));
    end

end


%% Método de Euler - Segunda Ordem

clc, clear all, close all

% Condições iniciais & finais
t0 = 0 ; t_end = 10;
x0 = 1; v0 = 0;
h = 0.01;

%Constantes 
K = 1 ; m = 1 ; w = K/m;

% Definir a função f(t, x, v)
f = @(t, x, v) -w*x;

% Método de Euler - 2ª ordem
[t, x, v] = metodo_euler(f, t0, x0, v0, h, t_end);

%Cálculo de Energia mecânica
Em = 0.5*m*v.^2 + 0.5*1*x.^2;

% Plot the solution
figure(1)
plot(t, x, '-', t, v, '-');
xlabel('t');
ylabel('x & v');
title('Solução pelo método de Euler');
legend('x(t)', 'v(t) = dx/dt');
grid on;

figure(2)

plot(t, Em, '-');
xlabel('t');
ylabel('Em');
title('Solução pelo método de Euler');
legend('Energia mecânica');
grid on;

function [t, x, v] = metodo_euler(f, t0, x0, v0, h, t_end)

    % Criar os arrays para armazenar os dados
    t = t0:h:t_end;  % Criar o array do tempo
    x = zeros(1,length(t));  % Inicializar variável
    v = zeros(1,length(t));  % Inciializar derivada da variável
    x(1) = x0;  % Condicção inicial da variável
    v(1) = v0;  % Condicção inicial da derivada da variável

    % Iteração para obter os valores (método de Euler)
    for n = 1:length(t)-1
        x(n+1) = x(n) + h * v(n);  % Update da variável
        v(n+1) = v(n) + h * f(t(n), x(n), v(n));  % Update da derivada   
    end
end

%% Método de Euler-Cromer

clc, clear all, close all

% Condições iniciais & finais
t0 = 0 ; t_end = 10;
x0 = 1; v0 = 0;
h = 0.01;

%Constantes 
K = 1 ; m = 1 ; w = K/m;

% Definir a função f(t, x, v)
f = @(t, x, v) -w*x;

% Método de Euler - 2ª ordem
[t, x, v] = metodo_euler_cromer(f, t0, x0, v0, h, t_end);

%Cálculo de Energia mecânica
Em = 0.5*m*v.^2 + 0.5*1*x.^2;

% Plot the solution
figure(1)
plot(t, x, '-', t, v, '-');
xlabel('t');
ylabel('x & v');
title('Solução pelo método de Euler');
legend('x(t)', 'v(t) = dx/dt');
grid on;

figure(2)

plot(t, Em, '-');
xlabel('t');
ylabel('Em');
title('Solução pelo método de Euler');
legend('Energia mecânica');
grid on;

function [t, x, v] = metodo_euler_cromer(f, t0, x0, v0, h, t_end)

    % Criar os arrays para armazenar os dados
    t = t0:h:t_end;  % Criar o array do tempo
    x = zeros(1,length(t));  % Inicializar variável
    v = zeros(1,length(t));  % Inciializar derivada da variável
    x(1) = x0;  % Condicção inicial da variável
    v(1) = v0;  % Condicção inicial da derivada da variável

    % Iteração para obter os valores (método de Euler)
    for n = 1:length(t)-1
        v(n+1) = v(n) + h * f(t(n), x(n), v(n));  % Update da derivada
        x(n+1) = x(n) + h * v(n+1);  % Update da variável           
    end
end

 
%% Método de Euler-Ímplícito
% O método de euler-ímplicito não tem uma solução geral, é necessário resolver primeiro o sistema analiticamente (consultar folha de apontamentos)

% v(k+1)=v(k)+f(k+1)*h  (1)
% y(k+1)=y(k)+v(k+1)*h  (2)
% se resolver normalmente substituir 1 em 2 e deixar y(k+1) só em função de coeficientes k
% se fizer por linsolve achar matriz dos coeficentes k+1 (A) e matriz dos coeficoentes k (igual a crank nicholson)

clc, clear all, close all

% Condições iniciais & finais 
t0 = 0 ; t_end = 10;
x0 = 1; v0=0;
h = 0.01;

%Constantes
K = k1 ; m = k2; w=sqrt(K/m);

% Criar os arrays para armazenar os dados
t = t0:h:t_end;  % Criar o array do tempo
N=length(t);
x = zeros(1,N);  % Inicializar variável
v = zeros(1,N);  % Inciializar derivada da variável
x(1) = x0;  % Condicção inicial da variável
v(1) = v0;  % Condicção inicial da derivada da variável

% Iteração para obter os valores (método de Euler-Cromer) - usando sistema de equações
for n = 1:N-1 
    x(n+1)=(x(n)+v(n)*h)/(1+w^2*h^2);  % Update da derivada
    v(n+1)=v(n)-w^2*x(k+1)*h;  % Update da variável
end

Et=1/2*m*v.^2+1/2*K*x.^2; %cálculo da energia mecânica

% Plot da solução x
figure(1)
plot(t, x, '-');
xlabel('t');
ylabel('x');
title('Solução pelo método de Euler-implicito');
grid on;

% Plot da solução v
figure(2)
plot(t, v, '-');
xlabel('t');
ylabel('v');
title('Solução pelo método de Euler-implicito');
grid on;

% Plot da solução Et
figure(2)
plot(t, Et, '-');
xlabel('t');
ylabel('Et');
title('Solução pelo método de Euler-implicito');
grid on;

%% Método de Euler-Ímplícito - usando linsolve

%substituir apenas o for loop por:

A= [1 -h; w^2*h 1]; % Matriz do sistema linear
for k=1:N-1    
    b= [x(k); v(k)]; % Vetor do lado direito da equação    
    
    Z= linsolve(A,b);  % Resolver o sistema linear

    % Atualizar os valores de x e v
    x(k+1)=Z(1,1);
    v(k+1)= Z(2,1);    
end

%% Método de Crank-Nicholson

% v(k+1)=v(k)+[f(k)+f(k+1)]*h/2
% y(k+1)=y(k)+(v(k)+v(k+1))*h/2
% Achar matriz coefientes k+1 (A) e matriz coeficientes k (B)

% O método de Crank-Nicholson não tem uma solução geral, é necessário resolver primeiro o sistema analiticamente (consultar folha de apontamentos)

%Exemplo para o OHS

%v=dx/dt
%funcao=dv/dt

clc, clear all, close all

% Condições iniciais & finais 
t0 = 0 ; t_end = 10;
x0 = 1; v0=0;
h = 0.01;

%Constantes k_
K = k1 ; m = k2; w=sqrt(K/m);

% Criar os arrays para armazenar os dados
t = t0:h:t_end;  % Criar o array do tempo
N=length(t);
x = zeros(1,N);  % Inicializar variável
v = zeros(1,N);  % Inciializar derivada da variável
x(1) = x0;  % Condicção inicial da variável
v(1) = v0;  % Condicção inicial da derivada da variável

A= [1 -h/2; w^2*h/2 1]; % Matriz do sistema linear
for k=1:N-1    
    b= [x(k)+h/2*v(k); v(k)-w^2*h/2*x(k)]; % Vetor do lado direito da equação    
    
    Z= linsolve(A,b);  % Resolver o sistema linear

    % Atualizar os valores de x e v
    x(k+1)=Z(1,1);
    v(k+1)= Z(2,1);    
end

Et=1/2*m*v.^2+1/2*K*x.^2; %cálculo da energia mecânica

% Plot da solução x
figure(1)
plot(t, x, '-');
xlabel('t');
ylabel('x');
title('Solução pelo método de Euler-implicito');
grid on;

% Plot da solução v
figure(2)
plot(t, v, '-');
xlabel('t');
ylabel('v');
title('Solução pelo método de Euler-implicito');
grid on;

% Plot da solução Et
figure(2)
plot(t, Et, '-');
xlabel('t');
ylabel('Et');
title('Solução pelo método de Euler-implicito');
grid on;

%% Método Runge-Kutta 2ª Ordem

% Exemplo de tabela de Butcher para RK2
% 0   |
% 1/2  | 1/2
% ----|---------
%     |  0   1

clc, clear all, close all

% Condições iniciais & finais
t0 = 0 ; t_end = 10;
x0 = 1; v0 = 0;
h = 0.1;

%Constantes 
K = 16 ; m = 1; w = sqrt(K/m) ;

% Funções das derivadas x e v
fv = @(t, x, v) -w^2 * x;  % dv/dt = -w^2 * x
fx = @(t, x, v) v;          % dx/dt = v

% Método Runge-Kutta 2ª ordem
[t, x, v] = runge_kutta(fv, fx, t0, x0, v0, h, t_end);

% Cálculo de Energia mecânica
Em = 0.5 * m * v.^2 + 0.5 * K * x.^2;  % Energia total

% Plot das soluções
figure(1)
plot(t, x, '-', t, v, '-');
xlabel('t');
ylabel('x & v');
title('Solução Runge-Kutta 2ª ordem');
legend('x(t)', 'v(t) = dx/dt');
grid on;

figure(2)
plot(t, Em, '-');
xlabel('t');
ylabel('Em');
title('Energia Mecânica');
legend('Energia mecânica');
grid on;

% Função do método Runge-Kutta 2ª ordem
function [t, x, v] = runge_kutta(fv, fx, t0, x0, v0, h, t_end)

    % Criar os arrays para armazenar os dados
    t = t0:h:t_end;  % Criar o array do tempo
    N = length(t);   % Número de passos
    x = zeros(1, N); % Inicializar a variável
    v = zeros(1, N); % Inicializar a derivada da variável
    x(1) = x0;       % Condição inicial da variável
    v(1) = v0;       % Condição inicial da derivada da variável

    % Iteração do método Runge-Kutta 2ª ordem
    for k = 1:N-1
        
        % Parte 1
        r1v = fv(t(k), x(k), v(k));  % dv/dt no ponto atual
        r1x = fx(t(k), x(k), v(k));  % dx/dt no ponto atual

        % Midpoint
        r2v = fv(t(k) + h/2, x(k) + r1x * h/2, v(k) + r1v * h/2);  % dv/dt no ponto intermediário
        r2x = fx(t(k) + h/2, x(k) + r1x * h/2, v(k) + r1v * h/2);  % dx/dt no ponto intermediário

        % Atualização de x e v
        x(k+1) = x(k) + h * r2x;  % Atualização de x usando o midpoint
        v(k+1) = v(k) + h * r2v;  % Atualização de v usando o midpoint
    end
end

%% Método de Runge-Kutta de 3ªordem

% 0   | 
% 1/2 | 1/2 
% 3/4  | 0   3/4   
% ----|----------------
%     |  2/9   1/3   4/9

clc, clear all, close all

% Condições iniciais & finais
t0 = 0 ; t_end = 50;
x0 = 1; v0 = 1;
h = 0.01;

%Constantes 
K = 1 ; m = 1; w = sqrt(K/m) ; alfa=-0.1;

% Funções das derivadas x e v
fv = @(t, x, v) -K/m*(x+2*alfa*x^3);  
fx = @(t, x, v) v;          % dx/dt = v

% Método Runge-Kutta 3ª ordem
[t, x, v] = runge_kutta_3(fv, fx, t0, x0, v0, h, t_end);

% Cálculo de Energia mecânica
Em=1/2*m*v.^2+K/2*x.^2.*(1+alfa*x.^2);  % Energia total

% Plot das soluções
figure(1)
plot(t, x, '-', t, v, '-');
xlabel('t');
ylabel('x & v');
title('Solução Runge-Kutta 4ª ordem');
legend('x(t)', 'v(t) = dx/dt');
grid on;

figure(2)
plot(t, Em, '-');
xlabel('t');
ylabel('Em');
title('Energia Mecânica');
legend('Energia mecânica');
grid on;

function [t, x, v] = runge_kutta_3(fv, fx, t0, x0, v0, h, t_end)

    % Criar os arrays para armazenar os dados
    t = t0:h:t_end;  % Criar o array do tempo
    N = length(t);   % Número de passos
    x = zeros(1, N); % Inicializar a variável
    v = zeros(1, N); % Inicializar a derivada da variável
    x(1) = x0;       % Condição inicial da variável
    v(1) = v0;       % Condição inicial da derivada da variável

    % Iteração do método Runge-Kutta 2ª ordem
    for k = 1:N-1
        
    % Parte 1
        k1v = fv(t(k), x(k), v(k));
        k1x = fx(t(k), x(k), v(k));

        % Parte 2
        k2v = fv(t(k) + h/2, x(k) + k1x * h/2, v(k) + k1v * h/2);
        k2x = fx(t(k) + h/2, x(k) + k1x * h/2, v(k) + k1v * h/2);

        % Parte 3
        k3v = fv(t(k) + 3*h/4, x(k) + k2x * 3*h/4, v(k) + k2v * 3*h/4);
        k3x = fx(t(k) + 3*h/4, x(k) + k2x * 3*h/4, v(k) + k2v * 3*h/4);

        % Update de x e v
        x(k+1) = x(k) + (h/9) * (2*k1x + 3*k2x + 4*k3x);
        v(k+1) = v(k) + (h/9) * (2*k1v + 3*k2v + 4*k3v);
    end

end
%% Runge-Kutta 3ª ordem só com uma função

clc; clear all; close all;

% Parâmetros do problema
h0 = 0.35; % cm
D = 0.35; % cm
d = 0.025; % cm
g = 9.81; % m/s²

% Tempo de simulação
t0 = 0; % s
tf = 30; % s
dt = 0.1; % Passo de tempo
t = t0:dt:tf; % Vetor de tempo
N = length(t); % Número de passos


fh=@(t,h) -sqrt(g*2) * (d/D)^2*sqrt(h);


% Inicialização das variáveis
h = zeros(1, N);
h(1) = h0;

for k = 1:N-1
        
        % Parte 1
        k1x = fh(t(k), h(k));

        % Parte 2
        k2x = fh(t(k) + dt/2, h(k) + k1x * dt/2);

        % Parte 3
        k3x = fh(t(k) + dt, h(k) -k1x*dt+ k2x *2*dt);

        % Update de x e v
        h(k+1) = h(k) + (dt/6) * (k1x + 4*k2x+k3x);
end

%% Método Runge-Kutta 4ª Ordem

% 0   | 
% 1/2 | 1/2 
%1/2  | 0   1/2   
% 1   | 0   0   1
% ----|----------------
%     |  1/6   1/3   1/3  1/6

clc, clear all, close all

% Condições iniciais & finais
t0 = 0 ; t_end = 10;
x0 = 1; v0 = 0;
h = 0.1;

%Constantes 
K = 16 ; m = 1; w = sqrt(K/m) ;

% Funções das derivadas x e v
fv = @(t, x, v) -w^2 * x;  % dv/dt = -w^2 * x
fx = @(t, x, v) v;          % dx/dt = v

% Método Runge-Kutta 4ª ordem
[t, x, v] = runge_kutta_4(fv, fx, t0, x0, v0, h, t_end);

% Cálculo de Energia mecânica
Em = 0.5 * m * v.^2 + 0.5 * K * x.^2;  % Energia total

% Plot das soluções
figure(1)
plot(t, x, '-', t, v, '-');
xlabel('t');
ylabel('x & v');
title('Solução Runge-Kutta 4ª ordem');
legend('x(t)', 'v(t) = dx/dt');
grid on;

figure(2)
plot(t, Em, '-');
xlabel('t');
ylabel('Em');
title('Energia Mecânica');
legend('Energia mecânica');
grid on;

function [t, x, v] = runge_kutta_4(fv, fx, t0, x0, v0, h, t_end)

    % Criar os arrays para armazenar os dados
    t = t0:h:t_end;  % Criar o array do tempo
    N = length(t);   % Número de passos
    x = zeros(1, N); % Inicializar a variável
    v = zeros(1, N); % Inicializar a derivada da variável
    x(1) = x0;       % Condição inicial da variável
    v(1) = v0;       % Condição inicial da derivada da variável

    % Iteração do método Runge-Kutta 2ª ordem
    for k = 1:N-1
        
    % Parte 1
    k1v = fv(t(k), x(k), v(k));
    k1x = fx(t(k), x(k), v(k));

    % Parte 2
    k2v = fv(t(k) + h/2, x(k) + k1x * h/2, v(k) + k1v * h/2);
    k2x = fx(t(k) + h/2, x(k) + k1x * h/2, v(k) + k1v * h/2);

    % Parte 3
    k3v = fv(t(k) + h/2, x(k) + k2x * h/2, v(k) + k2v * h/2);
    k3x = fx(t(k) + h/2, x(k) + k2x * h/2, v(k) + k2v * h/2);

    % Parte 4
    k4v = fv(t(k) + h, x(k) + k3x * h, v(k) + k3v * h);
    k4x = fx(t(k) + h, x(k) + k3x * h, v(k) + k3v * h);

    % Update de x e v
    x(k+1) = x(k) + (h/6) * (k1x + 2*k2x + 2*k3x + k4x);
    v(k+1) = v(k) + (h/6) * (k1v + 2*k2v + 2*k3v + k4v);
    end

end

%% Método de Shooting

clc, clear all, close all

%Constantes
T = 5.0E4 ; w = 1.0E5 ; a = 5.0E-8; L=3; h=0.01;

%Arrays
x = 0:h:L ; N = length(x) ; y = zeros(1,N) ; y(1) = 0 ; y(L) = 0;

%Relacionado a Shoooting
guess(1) = -2; guess(2) = -2.2 ; B=0 ; tol = 1E-4;

for is = 1:50
    dy = zeros(1,N);
    dy(1) = guess(is);

    for k=1:N-1
        dy(k+1) = dy(k) + h*(2*a*T*y(k) + a*w*x(k)*(L-x(k)));
        y(k+1) = y(k) + h*dy(k+1);
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



