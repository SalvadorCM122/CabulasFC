%% Achar um valor por interseção 

clc, clear all, close all

for i = 1:length(y)-1
    if y(i) > 0
        idx = i; %encontrar o primeiro indice negativo
        break
    end
end

%interp1(y(entre positivo a negativo) , x(entre positivo a negativo), valor
%que quero intersetar)

inter = interp1(y(idx-1:idx),x(idx-1:idx),  0);
disp(['Valor correspondente: ', num2str(inter)]);

%% Achar vários máximos / minimos sem interpolação

clc, clear all, close all

t = 0:0.1:10 ; x = sin(t);

Imax=find(islocalmax(x)); %armazenamento dos índices dos máximos locais
Imin = find(islocalmin(x)); %armazenamento dos índices dos minimos locais

tmax = t(Imax); xmax = x(Imax); %tempo e posição dos máximos locais
tmin = t(Imin); xmin = x(Imin); %tempo e posição dos mínimos locais

figure(2)
plot(t,x,'-')
hold on
plot(tmax, xmax, 'o')
hold on;
plot(tmin, xmin, 'o')
xlabel('t(s)');
ylabel('x(t)');
grid on;
hold off;

%% Periodo e Amplitude por interpolação de Lagrange (bonus de máximos)

clc, clear all, close all

%primeiro fazer metodo de euler-cromer ou crank-nicholson normalmente

% **Determinação da Amplitude e Período**
Im = find(islocalmax(x));  % Índices dos máximos locais
Tm = t(Im);  % Tempos correspondentes aos máximos
nI = length(Im);  

% Inicializar arrays para amplitudes e períodos refinados
Amp = zeros(1, nI - 2);
Texp = zeros(1, nI - 2);

% Aplicar interpolação de Lagrange aos picos
for j = 2:nI-1
    xm = t(Im(j-1:j+1));  % Três tempos vizinhos
    ym = x(Im(j-1:j+1));  % Três valores de x correspondentes
    max_values = lagr(xm, ym);  % Aplicar interpolação
    Amp(j-1) = max_values(2);  % Amplitude refinada
    Texp(j-1) = Tm(j) - Tm(j-1);  % Período entre máximos consecutivos
end

% Calcular amplitude e período médio
Amplitude = mean(Amp);
Periodo = mean(Texp);

% Exibir resultados
fprintf('Amplitude média: %.4f\n', Amplitude);
fprintf('Período médio: %.4f\n', Periodo);

function lagr=lagr(xm,ym)

% determinacao de o maximo de uma funcao discreta
% input: coordenadas de 3 pontos vizinhos de ordenadas maiores
% matrizes xm e ym
% output: coordenadas do ponto máximo (xmax,ymax)
%cálculo coeficientes para a interpolação quadrática

xab=xm(1)-xm(2);
xac=xm(1)-xm(3);
xbc=xm(2)-xm(3);

a=ym(1)/(xab*xac);
b=-ym(2)/(xab*xbc);
c=ym(3)/(xac*xbc);

xml=(b+c)*xm(1)+(a+c)*xm(2)+(a+b)*xm(3);
xmax=0.5*xml/(a+b+c);

xta=xmax-xm(1);
xtb=xmax-xm(2);
xtc=xmax-xm(3);

ymax=a*xtb*xtc+b*xta*xtc+c*xta*xtb;

lagr(1)=xmax;
lagr(2)=ymax;

end

%% ODE45

clc, clear all, close all

% Condições iniciais & finais
t0 = 0 ; t_end = 100;
x0 = 0.8; y0 = 0.3;
h = 0.1;

%Constantes 
a=2; b=0.74; c=0.5; 

options=odeset('RelTol',3E-14,'AbsTol',[1E-13 1E-13]);
[t, sol] = ode45(@(t, sol) func(t, sol, a, b, c), t0:h:t_end, [x0 y0], options);

x=sol(:,1);y=sol(:,2);

plot(t,x)
figure(2)
plot(t,y)

function derivadas=func(t,sol,a,b,c)
        derivadas=zeros(2,1);
        x=sol(1); y=sol(2);
        derivadas(1)=x*(1-x)-a*x*y/(x+y); %MUDAR A EXPRESSÃO DA DERIVADA
        derivadas(2)=-b*x*y/(x+y)-c*y; %MUDAR A EXPRESSÃO DA DERIVADA
end

%% ERROS GLOBAIS

% Erro global corresponde à diferença entre o valor da diferença entre a solução analítica y(tk) e a solução numérica y(k)

hh = [0.1,0.05,0.025,0.0125]; nh = length(hh) ; 
Erro = zeros(1,nh);

va = NaN; %meter a expressão que é dada para v analitico no determinado instante

for i=1:nh
    
    h = hh(i);

    for k=1:N  %LOOP DO MÉTODO USADO
        ...
    end

    Erro(i) = abs(va-v); % V usado no mesmo instante de va
end

lh = log(hh) ; le = log(Erro);

% Calcular a ordem do método
p = polyfit(lh,le,1);
ordem = p(1);

disp(['Ordem do método: ', num2str(ordem)]);

% Verificar se é linear
figure(1)
plot(lh,le, '-');
xlabel('log(hh)');
ylabel('log(erro)');
title('Erros para diferentes valores de h');
grid on;
