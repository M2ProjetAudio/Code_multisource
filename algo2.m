function theta_estimee=algo2(theta_init,J,Q,Ng,B)

theta_estimee=theta_init;
eta=1;
maxit=10;
count=0;
delta_L=Inf;

% a allouer
 gamma_bar=zeros(Q,Ng,B);
 gamma=zeros(Q,Ng,B);

while delta_L >= eta || count<maxit
    count=count+1;
    
    for q=1:Q
       if count==1
           theta_estimee(q)=theta_init(q);
       else
           %somme
           sum0=zeros(360,1);
           for ng=1:Ng
              sum1=zeros(360,1);
              for b=1:B
                  sum1=sum1+gamma(q,ng,b)*J(ng,b,:); %dimensions a verifier
              end
              sum0=sum0+sum1;
           end
           theta_estimee(q)=max(sum0);
       end
    end
    
    
    %LogLikelihood computation
    sum0=1;
    for ng=1:Ng
        sum1=0;
        for b=1:B
           sum2=0;
           for q=1:Q
              sum2=sum2+1/Q*exp(J(ng,b,theta_estimee(q))); 
           end
           sum1=sum1+sum2;
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
            for q=1:Q
                gamma(q,ng,b)=gamma_bar(q,ng,b)/sum;
            end
        end
    end
    
end