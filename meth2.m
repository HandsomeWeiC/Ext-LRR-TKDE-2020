function D=meth2(si,s,k,r,p,U)

        s0 = sqrt(s.^2+r).^p;
        re=find(si*s0<1);
        [l,m]=size(re);
        [n,m]=size(s0);
        s1=s0(re,1);

       Qtrg =(n-l)+sum(si*s1);

       d = p/2*(sqrt(s.^2+r).^(p-2));
       [m,n]=size(d);
       L=zeros(m,n);
       L(re,1)=d(re,1);
       d=L;
       D = k*Qtrg^(k-1)*U*diag(d)*U';
end