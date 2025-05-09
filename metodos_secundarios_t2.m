%% loop para a matriz A e vetor b
%Exemplo - teste315_16


for i = 1:NM %NM - número de popntos internos
    A(i,i) = (2+2/eta);
    if i > 1
        A(i,i-1) = -1;
    end
    if i ~= NM
        A(i,i+1) = -1;
    end
end

for n=1:Nt-1 %Nt - número de pontos domínio temporal
    for i=1:NM
        b(i)=T(i+2,n)+(2/eta-2)*T(i+1,n)+T(i,n)+dt*f(i+1)/eta;
    end
    b(1) = b(1) + T(1,n+1);   
    b(NM) = b(NM) + T(end,n+1); 
    T(2:Nx-1,n+1)=linsolve(A,b);
end

% A = c1*eye(NM) + diag(c2*ones(NM-1,1),1) +diag(c3*ones(NM-1,1),-1)

%------------------------------------------------------------------------------------------------%

%% sol_sist_trid

function y=sol_sist_trid(A,B)

% Encontra as solucoes y para um sistema de equacoes
% do tipo A y = B, onde A e' uma matriz tridiagonal

A_tri=A; 
b_tri=B;

n_tri=numel(b_tri);
d_tri=diag(A_tri);
c_tri=[diag(A_tri,1); 0];
a_tri=[0; diag(A_tri,-1)];

h_tri(1)=c_tri(1)/d_tri(1);
p_tri(1)=b_tri(1)/d_tri(1);

for i=2:n_tri
    h_tri(i)=c_tri(i)/(d_tri(i)-a_tri(i)*h_tri(i-1));
    p_tri(i)=(b_tri(i)-a_tri(i)*p_tri(i-1))/(d_tri(i)-a_tri(i)*h_tri(i-1));
end

y(n_tri)=p_tri(n_tri);
for i=n_tri-1:-1:1
    y(i)=p_tri(i)-h_tri(i)*y(i+1);
end

%------------------------------------------------------------------------------------------------%

%% crack usando LU

clc; clear all; close all

L=50; T0=0; Tf=0; k=0.93; c=0.094; p=8.9;
dt=0.2; tf=500; t=0:dt:tf; Nt=length(t); 
dx=0.5; x=0:dx:L; Nx=length(x);  
eta=k*dt/(c*p*dx^2);
NM=Nx-2;

T=zeros(Nx,Nt);
T(1,:)=T0; T(Nx,:)=Tf;
T(2:Nx-1,1)=100;

D1=(2/eta+2);
A= eye (NM);
A=D1*A; % matriz diagonal, diagonal principal
A(1,2)=-1; % 2º elemento da 1ª linha
for i=2:NM-1
    A(i,i-1)=-1; % diagonal superior
    A(i,i+1)=-1; % diagonal inferior
end
A (NM,NM-1)=-1; % penúltimo elemento da última linha

b=zeros(Nx-2,1);
D2=(2/eta-2);

[L,U,P]=lu(A);
% L-matriz inferior triangular; U-matriz superior triangular; P-matriz permutação

for n=1:Nt-1
    for i=1:Nx-2
    b(i) = T(i,n)+D2*T(i+1,n)+T(i+2,n);
    end

    b(1)=b(1)+T(1,n+1); % é preciso adicionar a CF
    b(NM)=b(NM)+T(Nx,n+1); % CF
    y=L\b; 
    T(2:Nx-1,n+1)=U\y; 
end


figure(1)
contourf(x,t,T')


