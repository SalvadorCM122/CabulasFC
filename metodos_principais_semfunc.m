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

