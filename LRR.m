% uncleaned version, just for reference.

% isbias=0: min_{Qi}  \sum_i  (rank(X*Qi))^k
% isbias=1: min_{Qi}  \sum_i  (rank((X-bi*1')*Qi))^k
function [Q, Qd, obj, b] = LRR(X, c,si, k,isdisc, isbias)
% X: dim*n data matrix
% c: number of clusters
% si: cofficient of singular value for the matrix XQj
% k: parameter, k>=2
% isdisc: 1 if using discrete constraint of Qi and 0 otherwise
% isbias: 1 if using bias and 0 otherwise
% Q: n*c indicator matrix
% Qd: n*c discrete indicator matrix
% obj: objective values during the iteration
% b: c*1 bias vector if isbias=1

% Ref: Feiping Nie, Wei Chang. Robust Subspace Clustering with Low-Rank Structure Constraint. IEEE TKDE 2020.

p=1;
r = 0.000000001;   % regularized  parameter
if p > 0.95
    inMaxIter = 50;
else
    inMaxIter = 50;
end;
[dim, n] = size(X);
% inti Q, rand init
A = rand(n,c);
Q = A./(sum(A,2)*ones(1,c));
if isdisc == 1
    [min_val, qMat] = min(A, [], 2);
    Q = zeros(n,c);
    for t = 1:n     
        Q(t, qMat(t)) = 1;
    end
end;

%Q=1/c*ones(n,c);
%n1=50;n2=50;p1=1; Q=[[p1*ones(n1,1),(1-p1)*ones(n1,1)];[(1-p1)*ones(n2,1),p1*ones(n2,1)]];
%Q=[ones(n,1),zeros(n,1)];
%Q = kron(eye(5),ones(50,1));

obj = zeros(inMaxIter, 1);
if dim > n
    z = zeros(dim-n,1);
else
    z = [];
end;
b = zeros(dim,c);
for iter = 1: inMaxIter 
    % fix Q, update D
    Qtr = 0;
    for g = 1:c
        XQ = (X-b(:,g)*ones(1,n))*diag(Q(:,g));
        [U,S,V] = svd(XQ);
        s = [diag(S); z]; 
        s0 = sqrt(s.^2+r).^p;
        Qtrg = sum(min(si*s0,1));
        Qtr = Qtr + Qtrg^k;
        D{g} = meth2(si,s,k,r,p,U);
    end
    
    % fix D, update Q
    for g = 1:c
        X1 = X-b(:,g)*ones(1,n);
        A(:,g) = diag(X1'*D{g}*X1);
    end
    if isdisc == 1
        [min_val, qMat] = min(A, [], 2);
        Q = zeros(n,c);
        for t = 1:n     
            Q(t, qMat(t)) = 1;
        end
    else
        for t = 1:n
        a = A(t,:);
        Q(t,:) = 1/sum(1./a)*1./a;
        end;
    end;
    
    % update b
    if isbias == 1
        QQ = Q.^2;
        b = (X*QQ)./((ones(dim,1)*sum(QQ))+eps);
    else
        b = zeros(dim,c);
    end;
    
    % calculate obj   
    obj(iter) = Qtr;
end

if isdisc == 1
    Qd = Q;
else
    [min_val, qMat] = max(Q, [], 2);
    Qd = zeros(n,c);
    for t = 1:n     
        Qd(t, qMat(t)) = 1;
    end
end;





