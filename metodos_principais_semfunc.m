% Método de euler - ODE de primeira ordem

clc, clear all, close all

% Condições iniciais & finais 
t0 = 0 ; t_end = 10;
x0 = 1;
h = 0.01;

%Constantes k_
a = k1 ; b = k2;

%Definir a função f(t, x)
f = @(t, x) x+t; %se não quiser esta opção, substituir a função no método

% Criar os arrays para armazenar os dados
t = t0:h:t_end; % Criar o array do tempo
x = zeros(1,length(t)); % Inicializar variável
x(1) = x0; % Condicção inicial da variável

% Iteração método de Euler   
for n = 1:length(t)-1
    x(n+1) = x(n) + h * f(t(n), x(n));
end

% Plot da solução
plot(t, x, '-');
xlabel('t');
ylabel('x');
title('Solução pelo método de Euler');
grid on;


%% Método de Euler - sistemas de ODE de primeira ordem
%Euler funciona igual a Euler-cromer mas  não conserva tão bem a energia, não é muito estável 
%para certos tipos de equações diferenciais, especialmente quando há oscilações ou sistemas rígidos.

%v=dx/dt
%funcao=dv/dt

clc, clear all, close all

% Condições iniciais & finais 
t0 = 0 ; t_end = 10;
x0 = 1; v0=0;
h = 0.01;

%Constantes k_
a = k1 ; b = k2;

%Definir a função f(t, x)
f = @(t, x) funcao; %substituir funcao pela pretendida (exemplos no inicio do documento)

% Criar os arrays para armazenar os dados
t = t0:h:t_end;  % Criar o array do tempo
N=length(t);
x = zeros(1,N);  % Inicializar variável
v = zeros(1,N);  % Inciializar derivada da variável
x(1) = x0;  % Condicção inicial da variável
v(1) = v0;  % Condicção inicial da derivada da variável

% Iteração para obter os valores (método de Euler)
for n = 1:N-1     
    x(n+1) = x(n) + h * v(n);  % Update da variável  
    v(n+1) = v(n) + h * f(t(n), x(n), v(n));  % Update da derivada
end

% Plot da solução x
figure(1)
plot(t, x, '-');
xlabel('t');
ylabel('x');
title('Solução pelo método de Euler');
grid on;

% Plot da solução v
figure(2)
plot(t, v, '-');
xlabel('t');
ylabel('v');
title('Solução pelo método de Euler');
grid on;

%% Método de Euler-Cromer

% Método de euler-cromer pode ser aplicado a uma ODE de qualquer ordem, exceto em apenas uma de primeira ordem, pois 
%uma ODE de ordem n pode ser escrita em n ODE’s de primeira ordem.
%Euler-Cromer é um método de primeira ordem. Funciona bem e é mais usado para sistemas físicos conservativos.

%v=dx/dt
%funcao=dv/dt

clc, clear all, close all

% Condições iniciais & finais 
t0 = 0 ; t_end = 10;
x0 = 1; v0=0;
h = 0.01;

%Constantes k_
a = k1 ; b = k2;

%Definir a função f(t, x)
f = @(t, x) funcao; %substituir funcao pela pretendida (exemplos no inicio do documento)

% Criar os arrays para armazenar os dados
t = t0:h:t_end;  % Criar o array do tempo
N=length(t);
x = zeros(1,N);  % Inicializar variável
v = zeros(1,N);  % Inciializar derivada da variável
x(1) = x0;  % Condicção inicial da variável
v(1) = v0;  % Condicção inicial da derivada da variável

% Iteração para obter os valores (método de Euler-Cromer)
for n = 1:N-1 
    v(n+1) = v(n) + h * f(t(n), x(n), v(n));  % Update da derivada
    x(n+1) = x(n) + h * v(n+1);  % Update da variável
end

% Plot da solução x
figure(1)
plot(t, x, '-');
xlabel('t');
ylabel('x');
title('Solução pelo método de Euler');
grid on;

% Plot da solução v
figure(2)
plot(t, v, '-');
xlabel('t');
ylabel('v');
title('Solução pelo método de Euler');
grid on;

%% Método de Kunge-Kutta de 3ª ordem (RK3)
%Mesma utilização que o RK2, mas mais preciso


% Método de Runge-Kutta de 3ªordem


% Exemplo de tabela de Butcher para RK3

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
        k3v = fv(t(k) + 3*h/4, x(k) +(h*0)*r1x + k2x * 3*h/4, v(k) + (h*0)*r1v + k2v * 3*h/4);
        k3x = fx(t(k) + 3*h/4, x(k) + (h*0)*r1x + k2x * 3*h/4, v(k) + (h*0)*r1v + k2v * 3*h/4);

        % Update de x e v
        x(k+1) = x(k) + (h/9) * (2*k1x + 3*k2x + 4*k3x);
        v(k+1) = v(k) + (h/9) * (2*k1v + 3*k2v + 4*k3v);
    end

end
