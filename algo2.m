function theta_estimee=algo2(theta_init,J,Q,Ng,B)

theta_estimee=theta_init;
eta=1;
maxit=100;
count=0;
delta_L=Inf;

% a allouer
 gamma_bar=zeros(Q,Ng,B);
 gamma=zeros(Q,Ng,B);

while delta_L >= eta || count<maxit
    count=count+1;
    fprintf('Iteration numero %2d\n',count);
    for q=1:Q
       if count==1
           theta_estimee(q)=theta_init(q);
       else
           %somme
           sum0=zeros(360,1);
           for ng=1:Ng
              sum1=zeros(360,1);
              for b=1:B
                  int=squeeze(J(ng,b,:));
                  sum1=sum1+gamma(q,ng,b)*int; %dimensions a verifier
              end
              sum0=sum0+sum1;
           end
           int=round(real(max(sum0)));
           theta_estimee(q)=mod(int,360)*(int~=0)+360*(int==0);
       end
    end
    
    
    %LogLikelihood computation
    sum0=1;
    for ng=1:Ng
        sum1=0;
        for b=1:B
           sum2=0;
           for q=1:Q
               %theta_estimee(q)
              sum2=sum2+1/Q*exp(J(ng,b,theta_estimee(q))); 
           end
           sum1=sum1+log(sum2);
        end
        sum0=sum0+sum1;
    end  
    if count > 1
        L_theta_chapeau_precedent=L_theta_chapeau;
        L_theta_chapeau=sum0;
        delta_L=(L_theta_chapeau-L_theta_chapeau_precedent)...
            /L_theta_chapeau_precedent;
    else
        L_theta_chapeau=sum0;
    end
    fprintf('\t Delta_L = %3f \n',delta_L);
    theta_precedent=theta_estimee;
    % Expectation step
    for ng=1:Ng
        for b=1:B
            for q=1:Q
                gamma_bar(q,ng,b)=exp(J(ng,b,theta_precedent(q)));
            end
            sum=0;
            for q=1:Q
                sum=sum+gamma_bar(q,ng,b);
            end
            %for q=1:Q
         %   if sum ~=0
                % gamma(q,ng,b)=gamma_bar(q,ng,b)/sum;
                gamma(:,ng,b)=gamma_bar(:,ng,b)/sum;
          %  else
          %      gamma(:,ng,b)=gamma_bar(:,ng,b)*0;
         %   end
            % end
        end
    end
    
end