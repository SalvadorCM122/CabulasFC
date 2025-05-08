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

%% Diferenças finitas (condiçao de neumann e dirichlet)

clc; clear; close all;

% Dados do problema
Q = 2.1e6;           % [W/m^3]
lambda = 0.1;        % [W/(m.K)]
R = 1e-3;            % [m]
T_ext = 20;          % Temperatura na superfície externa [ºC]

% Parâmetros da malha
h = 0.00001;         % Passo radial
r = 0:h:R;
Nx = length(r);      % Número total de pontos
N = Nx - 2;          % Número de pontos internos

% Inicializa matriz A e vetor b
A = zeros(N, N);
b = zeros(N, 1);

% Construção de A e b com diferenças finitas centradas
for i = 1:N
    ri = r(i+1); % Ponto interior correspondente (pula r=0 e r=R)
    A(i,i) = -2 / h^2;
    if i > 1
        A(i,i-1) = 1/h^2 - 1/(2*h*ri);
    end
    if i < N
        A(i,i+1) = 1/h^2 + 1/(2*h*ri);
    end
    b(i) = -Q / lambda; % constante do termo fonte
end

% Condições de fronteira:
% T'(0) = 0 ⇒ T_1 = T_0 → modificar primeira equação
A(1,1) = -2 / h^2;
A(1,2) = 2 / h^2;

% T(R) = 20 ºC ⇒ última equação (i = N) ajusta b(N)
b(N) = b(N) - (1/h^2 + 1/(2*h*r(end-1))) * T_ext;

% Resolve o sistema linear
T_internal = linsolve(A, b);

% Concatena com as condições de fronteira
T = [T_internal(1); T_internal; T_ext];

% Plot do perfil de temperatura
plot(r, T, 'LineWidth', 2)
xlabel('r [m]')
ylabel('Temperatura [ºC]')
title('Perfil de Temperatura na Resistência Elétrica Cilíndrica')
grid on

% Encontra valor máximo de temperatura
[T_max, idx_max] = max(T);
r_max = r(idx_max);
fprintf('Temperatura máxima: %.2f ºC ocorre em r = %.6f m\n', T_max, r_max);

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

%% Crank Nicholson - barra temperatura
clc; close all; clear all

%constantes
k=0.93; c=0.094; p=8.9;
L=50; dx=0.1; dt=1;
x=0:dx:L; t=0:dt:500; Nx=length(x);Nt=length(t);
NM=Nx-2;
eta=k*dt/(c*p*dx^2);
f=2*exp(-(x-L/2).^2);

T=zeros(Nx,Nt);
T(:,1)=linspace(0,20,Nx);
T(1,:)=0;
T(Nx,:)=20;

A=zeros(NM); b=zeros(NM,1);

for i = 1:NM
    A(i,i) = (2+2/eta);
    if i > 1
        A(i,i-1) = -1;
    end
    if i ~= NM
        A(i,i+1) = -1;
    end
end

for n=1:Nt-1
    for i=1:NM
        b(i)=T(i+2,n)+(2/eta-2)*T(i+1,n)+T(i,n)+dt*f(i+1)/eta;
    end
    b(1) = b(1) + T(1,n+1);   
    b(NM) = b(NM) + T(end,n+1); 
    T(2:Nx-1,n+1)=linsolve(A,b);
end

figure(1)
contourf(x,t,T')

figure(2)
mesh(x,t,T')

%% Diferenças finitas centradas (espaço)

f' = (T(i+1,n) - T(i-1,n)) / 2h
f'' = (T(i-1,n) - 2*T(i,n) + T(i+1,n)) / h^2

%% Diferenças finitas avançadas (temporal)

f' = T(i,n+1) - T(i,n) / h
f'' = (T(i,n-1) - 2*T(i,n) + T(i,n+1)) / h^2

%% Formula das diferenças por Crank Nicolson

f'' = (T(i-1,n+1) - 2*T(i,n+1) + T(i+1,n+1) + T(i-1,n) - 2*T(i,n) + T(i+1,n)) / (2*h^2) % Se for centrada

%% Relaxação de Jacobi (exemplo do ex3.4 FR2)
clc; clear all; close all

Ms = [21, 41, 61, 81, 101, 121];  % diferentes valores de M
numIter = zeros(1,length(Ms));       % número de iterações para cada M
tol = 1e-5;                       % tolerância de convergência

% loop para os valores de M
for m = 1:length(Ms)
    M = Ms(m);
    x = linspace(-1, 1, M);
    y = linspace(-1, 1, M);
    h = 2/(M - 1); %espaçamento da malha
    f = zeros(M, M); %matriz fonte

    % Montar f(x,y)
    for i = 1:M
        for j = 1:M
            f(i,j) = -2*(2 - x(j)^2 - y(i)^2);
        end
    end

    % Inicializar T
    Told = zeros(M, M); %solução iteração anterior
    Tnew = Told; %nova solução com zeros

    %loop iterativo de Jacobi
    itmax = 10000;
    for l = 1:itmax
        for i = 2:M-1
            for j = 2:M-1
                Tnew(i,j) = 0.25 * (Told(i-1,j) + Told(i+1,j) + Told(i,j-1) + Told(i,j+1) - h^2 * f(i,j));
            end
        end

        % critério de paragem (erro relativo entre duas iterações consecutivas)
        num = sqrt(sum(sum((Tnew - Told).^2)));
        den = sqrt(sum(sum(Tnew.^2)));
        if num / den < tol
            break
        end
        Told = Tnew;
    end

    numIter(m) = l;  % guardar iterações
end

% Malha para gráfico
[X, Y] = meshgrid(x, y);

% Gráficos
figure(1)
contourf(X, Y, Tnew, 20)
colorbar
title('Contorno de T(x,y)')

figure(2)
mesh(X, Y, Tnew)
xlabel('x'); ylabel('y'); zlabel('T(x,y)')
title('Distribuição de T(x,y)')

% Plot log-log
logM = log(Ms);
logIter = log(numIter);

figure (3)
plot(logM, logIter, 'o-')
xlabel('log(M)')
ylabel('log(nº de iterações)')
title('Convergência vs Tamanho da malha')
grid on

% Regressão linear para declive
p = polyfit(logM, logIter, 1);
slope = p(1);
fprintf('Declive da reta (taxa de convergência): %.2f\n', slope);
