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
