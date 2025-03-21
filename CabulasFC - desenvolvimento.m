
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

% Métodos que faltam: Euler Implicito e formulas de Erros Globais

% PS: O método de Euler Implicito não foi adicionado pois ele depende de
% cada exercicio. Existe sim uma fórmula geral mas a Dr.Sofia Latas não 
% iria gostar.

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


%% Crank-Nicolson - 1ª ordem

clc, clear all, close all

% Condições iniciais & finais
t0 = 0 ; t_end = 10;
x0 = 1;
h = 0.01;

% Definir a função f(t, y)
f = @(t, x) -x;  

% Método Crank-Nicolson
[t, x] = crank_nicolson_first_order(f, t0, x0, h, t_end);

% Plot da solução
figure(1)
plot(t, x, '-');
xlabel('t');
ylabel('y');
title('Solução Crank-Nicolson para ODE de 1ª ordem');
legend('y(t)');
grid on;

function [t, x] = crank_nicolson_first_order(f, t0, x0, h, t_end)
    % Criar os arrays para armazenar os dados
    t = t0:h:t_end;  % Array do tempo
    N = length(t);   % Número de passos
    x = zeros(1, N); % Inicializar a variável
    x(1) = x0;       % Condição inicial da variável

    % Iteração do método Crank-Nicolson
    for n = 1:N-1
        % Tempo e y atuais
        tn = t(n);
        xn = x(n);

        % Definição do sistema linear implícito
        A = 1 + (h/2);  % Matriz dos coeficientes
        B = xn + (h/2) * f(tn, xn);  % Lado direito

        % Resolução do sistema linear
        x(n+1) = B / A;  
    end
end


%% Crank-Nicolson - 2ª ordem

clc, clear all, close all

% Condições iniciais & finais
t0 = 0 ; t_end = 10;
x0 = 1; v0 = 0;
h = 0.01;

%Constantes 
a = 1 ; b = 2; c = 3; m = 1;

% Definir a função f(t, x, v)
f = @(t, x, v) -2*v - x;

%  Método Crank-Nicolson
[t, x, v] = crank_nicolson(f, t0, x0, v0, h, t_end);

%Cálculo de Energia mecânica
Em = 0.5*m*v.^2 + 0.5*1*x.^2;

% Plot das soluções
figure(1)
plot(t, x, '-', t, v, '-');
xlabel('t');
ylabel('x & v');
title('Solução Crank-Nicolson');
legend('x(t)', 'v(t) = dx/dt');
grid on;

figure(2)
plot(t, Em, '-');
xlabel('t');
ylabel('Em');
title('Solução pelo método de Euler');
legend('Energia mecânica');
grid on;

function [t, x, v] = crank_nicolson(f, t0, x0, v0, h, t_end)

    % Criar os arrays para armazenar os dados
    t = t0:h:t_end;  % Criar o array do tempo
    N = length(t);   % Número de passos
    x = zeros(1, N); % Inicializar a variável
    v = zeros(1, N); % Inicializar a derivada da variável
    x(1) = x0;       % Condição incial da variável
    v(1) = v0;       % Condição incial da derivada da variável

    % Crank-Nicolson iteration
    for n = 1:N-1
        % tempo, x, v atuais
        tn = t(n);
        xn = x(n);
        vn = v(n);

        % Definição das matrizes do sistema
        A = [1, -h/2; (h/2), 1];  % Matriz dos coeficientes
        B = [xn + (h/2) * vn; vn + (h/2) * f(tn, xn, vn)];  % Matriz do lado direito

        % Resolução do sistema linear
        Z = linsolve(A,B);

        % Update da variável e a sua derivada
        x(n+1) = Z(1);
        v(n+1) = Z(2);
    end
end

%% Método Runge-Kutta 2ª Ordem

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

%% Método Runge-Kutta 4ª Ordem

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

%% Função islocalmax (achar picos numa função)

Im=find(islocalmax(V)>0); %armazenamento dos índices dos picos locais
Tm=t(Im); armazenamento dos tempos que correspondem a cada um desses ppicos
nI=length(Im); % numero de picos
Texp=diff(Tm); %retorna automaticamente todas as diferenças entre os picos consecutivoss
T=mean(Texp) %faz a média das diferenças retornando o período médio

%% ERROS GLOBAIS

% Erro global corresponde à diferença entre o valor da diferença entre a solução analítica y(tk) e a solução numérica y(k)
Ek=y(tk)-y(k)

%ainda não sei o que fazer aqui
