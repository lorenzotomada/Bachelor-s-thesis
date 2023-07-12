close all
clear all

n=14;
% A(1,1)=0;
% for i=1:n
%     A(i,i)=(i+1)/(n+1);
% end
% v(1,1)=1; %definiamo il vettore [1,...,1]
% for i=1:n
% v(i,1)=1;
% end
% b=exp(A)\v; %scegliamo b in maniera tale che f(A)b=v, che conosciamo

f=@(z) exp(z);
tol=0; %tolleranza scelta per il processo di Arnoldi, dovuta a
% eventuali cancellazioni numeriche

beta=norm(b,2);
q1=b/beta;
[Q,H,HK1K]=Arnoldi(A,q1,tol);
[m,~]=size(H); %m è il passo della computazione raggiunto

errore(1)=0; %inizializziamo l'errore a 0
e1(1)=1; %definiamo il vettore e_1 della base canonica
ek(1)=0; %definiamo il vettore e_k della base canonica
resto(1)=0; %inizializziamo il resto r=f-p a 0
stima(1)=0; %inizializziamo la stima a 0

for k=1:m %approssimazioni di Arnoldi
    Hk=H(1:k,1:k);
    spettro=eig(Hk)
    ptilde=polinomio_interpolante(f,spettro);
    %aggiungere controllo
    ptilde=sym2poly(ptilde); %otteniamo un vettore formato dai coefficienti
%     del polinomio che interpola f sullo spettro di H_k
    for i=1:length(ptilde)
        pA=ptilde(length(ptilde)+1-i)*A^(i-1);
        pHk=ptilde(length(ptilde)+1-i)*Hk^(i-1);
    end
    pAb=pA*b; %calcolo di p(A)b 
    errore(k)=norm(f(A)*b-pAb,2);
    rkA=f(A)-pA;
    rkHk=f(Hk)-pHk;
    resto(k)=norm(b,2)*(norm(rkA,2)+norm(rkHk,2));
    if k~=1
        e1(k,1)=0;
    end
    ek(k,1)=1;
    stima(k)=norm(b)*HK1K(k)*abs(ek'*f(Hk)*e1);
end

plot((1:m),resto(1:m),'g')
hold on
plot((1:m),errore(1:m),'r')
plot((1:m), stima(1:m),'b')
legend('Maggiorazione','Errore','Stima')
xlabel('k')