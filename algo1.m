function J=algo1(Z,Qn,V,B,Ng,Nf,Nt,sourcetype)

J=zeros(Ng,B,Nt);
C=zeros(2,2,Ng);
fprintf('\t algo 1')
for b=1:B
   for ng=1:Ng
      sum=zeros(2,1);
      for nf=1:Nf
          Ztemp=Qn(:,:,b)\Z(:,:,b,ng,nf);
         %Ztemp=Z(:,:,b,ng,nf);
          sum=sum+Ztemp;
        %  sum=sum+Z(:,:,b,ng,nf);
      end
      C(:,:,ng)=sum*sum'/Nf;
   end
   
   for th=1:Nt
      Vt=Qn(:,:,b)\V(:,:,b,th);
       %Vt=V(:,:,b,th);
       P=Vt/(Vt'*Vt)*Vt';        
       Pp=eye(2,2)-P;
       for ng=1:Ng
           if strcmp(sourcetype,'gauss')
%                J(ng,b,th)=-Nf*(log(det(P*C(:,:,ng)*P+Pp))+trace(Pp*C(:,:,ng)));
                J(ng,b,th)=-Nf*log(det(P*C(:,:,ng)*P+Pp*trace(Pp*C(:,:,ng))));
           else
               J(ng,b,th)=-2*Nf*log(.5*trace(Pp*C(:,:,ng)));
           end
        1;
       end   
   end
end
fprintf(' ... criteres J obtenus \n')

