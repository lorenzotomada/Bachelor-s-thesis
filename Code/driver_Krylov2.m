% esperimento numerico con la matrice derivante dalla discretizzazione del Laplaciano nel caso n=1


close all
clear all

%Delta=[0,1]
N=20;
h=1/N;
c=1e-1;
uDeltaN0=ones(N-1,1);

ADeltaN=zeros(N-1,N-1); %inizializziamo ADeltaN a 0
for i=1:N-1
    ADeltaN(i,i)=-2;
    if i~=N-1
        ADeltaN(i+1,i)=1;
        ADeltaN(i,i+1)=1;
    end
end
T=0.1; %istante di tempo in cui fare la valutazione
A=T*c*ADeltaN/(h^2);

b=uDeltaN0;
beta=norm(b,2);
q1=b/beta;
f='exp';
tol=1e-8; %tolleranza scelta per il processo di Arnoldi, dovuta
%a eventuali cancellazioni numeriche

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
    nR(k)=beta*HK1K(k)*abs(ek'*funm2(Hk, f)*e1); %calcolo dell'esponenziale
    %attraverso la funzione funm2, che fa uso dell'algoritmo di
    %Schur-Parlette;
    err(k)=norm(expm(A)*b-Qk*expm(Hk)*e1*beta,2);
end 

semilogy(1:m, nR,'b'), hold on
xlabel('k')
semilogy(1:m, err,'r')
legend('Residuo generalizzato', 'Errore stimato')
