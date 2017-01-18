function J=algo1(Z,Q,V,B,Ng,Nf,Nt)
% Q (:,:,b)
J=zeros(Ng,B,Nt);
C=zeros(2,2,Ng);
for b=1:B
   for ng=1:Ng
      sum=zeros(2,1);
      for nf=1:Nf
          Ztemp=Q(:,:)\Z(:,:,b,ng,nf);
          sum=sum+Ztemp;
      end
      C(:,:,ng)=sum*sum'/Nf;
   end
   
   for th=1:Nt
       Vt=Q(:,:)\V(:,:,b,th);
       P=Vt/(Vt'*Vt)*Vt';
       Pp=eye(2,2)-P;
       for ng=1:Ng
          J(ng,b,th)=-Nf*(log(det(P*C(:,:,ng)*P+Pp))+trace(Pp*C(:,:,ng)));
       end   
   end
end
fprintf('algo 1 termine  \n')

