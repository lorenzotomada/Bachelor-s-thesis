% esperimento numerico con una matrice random 20x20


close all
clear all

n=20;
A=rand(n,n);
b=rand(n,1);
if norm(b,2)==0
    return
end
beta=norm(b,2);
q1=b/beta;
f='exp';
tol=1e-13; %tolleranza scelta per il processo di Arnoldi, dovuta a
% eventuali cancellazioni numeriche

[Q,H,HK1K]=Arnoldi(A,q1,tol);
[m,~]=size(H); %m è il passo della computazione raggiunto
nR(1)=0; %inizializziamo il residuo generalizzato a 0
ek(1,1)=1; %definiamo il vettore e_k della base canonica
e1(1,1)=1; %definiamo il vettore e_1 della base canonica
err(1)=0; %errore rispetto alla matrice ottenuta usando la funzione "expm"
%di MATLAB

for k=1:m
    Qk=Q(:,1:k);
    Hk=H(1:k,1:k);
    if k~=1
        ek(k,1)=1; %vettore e_k della base canonica
        ek(k-1,1)=0;
        e1(k,1)=0;
    end
    nR(k)=beta*HK1K(k)*abs(ek'*expm(Hk)*e1); %calcolo dell'esponenziale
    %attraverso la funzione funm2, che fa uso dell'algoritmo di
    %Schur-Parlette;
    err(k)=norm(expm(A)*b-Qk*expm(Hk)*e1*beta,2);
end 

semilogy(1:m, nR,'b'), hold on
xlabel('k')
semilogy(1:m, err,'r')
legend('Residuo generalizzato', 'Errore stimato')
