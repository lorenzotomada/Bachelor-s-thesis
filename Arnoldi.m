function [Q,H,HK1K] = Arnoldi(A,b,tol)
%Funzione che implementa l'algoritmo di Arnoldi. Prende in input una
%matrice A, un vettore b, con norm(b,2)=1; restituisce in output le
%matrici Q_k, H_k, come precedentemente definite, e la matrice HK1K,
%i cui elementi sono i coefficienti h_{k+1,k} ottenuti all'iterazione
%k-esima

Q(:,1)=b; %q1=b/norm(b,2); qua norm(b,2)=1
n=length(b); %dimensione del problema

for k=1:n
    qk=Q(:,k);
    z=A*qk;
    for i=1:k
        qi=Q(:,i);
        H(i,k)=dot(qi,z);
        z=z-H(i,k)*qi;
    end
    HK1K(k)=norm(z,2);
    if HK1K(k)<=tol %se h_{k+1,k}=0
        return
    end
    H(k+1,k)=HK1K(k);
    Q(:,k+1)=z/HK1K(k);
end
end

