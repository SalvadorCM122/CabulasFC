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


