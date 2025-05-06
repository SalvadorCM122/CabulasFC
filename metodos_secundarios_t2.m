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


