function theta_estimee=algo2(theta_init,J,Q,Ng,B,thetaArg)

theta_estimee=theta_init;
eta=.02;
maxit=5;
count=0;
delta_L=Inf;

% a allouer
gamma_bar=zeros(Q,Ng,B);
gamma=zeros(Q,Ng,B);
fprintf('\t algo 2\n');
while delta_L >= eta || count<maxit
    count=count+1;
    fprintf('\t\tIteration numero %2d\n',count);
    for q=1:Q
        if count==1
            theta_estimee(q)=theta_init(q);
        else
            %somme
            sum0=zeros(360,1);
            for ng=1:Ng
                sum1=zeros(360,1);
                for b=1:B
                    tmp=squeeze(J(ng,b,:));
                    sum1=sum1+gamma(q,ng,b)*tmp;
                end
                sum0=sum0+sum1;
            end
            %discutable ....
%             if q>1
%                 sum0(theta_estimee(q-1))=[];
%                 [~,theta_estimee(q)]=max(sum0);
%                 while abs(  theta_estimee(q)-theta_estimee(q-1))<5
%                     sum0(theta_estimee(q))=[];
%                     [~,theta_estimee(q)]=max(sum0);
%                 end
%             else
                [~,theta_estimee(q)]=max(sum0);
%             end
            % theta_estimee(q)=thetaArg(idxmax);
            %int=round(real(max(sum0)));
            %theta_estimee(q)=mod(int,360)*(int~=0)+360*(int==0);
            
        end
    end
    
    
    %LogLikelihood computation
    sum0=0;
    for ng=1:Ng
        sum1=0;
        for b=1:B
            sum2=0;
            for q=1:Q
                %    theta_estimee(q)
                
                sum2=sum2+exp(J(ng,b,theta_estimee(q)))/Q;
                %sum2=sum2+J(ng,b,theta_estimee(q))/Q;
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
    fprintf('\t\t\t\t Estimation= %.0f \n',thetaArg(theta_estimee)*180/pi);
    fprintf('\t\t\t\t\t\t\t\t\t\t  Delta_L = %.0f %%\n',delta_L*100);
    theta_precedent=theta_estimee;
    % Expectation step
    for ng=1:Ng
        for b=1:B
%             fprintf('\n\n')
            for q=1:Q
                %fprintf('\t\t poids brut: gamma_bar(q=%d,ng=%d,b=%d)=%.2f\n',q,ng,b,exp(J(ng,b,theta_precedent(q))))
                gamma_bar(q,ng,b)=exp(J(ng,b,theta_precedent(q)));
            end
            sum1=0;
            for q=1:Q
                sum1=sum1+gamma_bar(q,ng,b);
            end
            
            %for q=1:Q
            %             if sum ~=0
            % gamma(q,ng,b)=gamma_bar(q,ng,b)/sum;
            gamma(:,ng,b)=gamma_bar(:,ng,b)/sum1;
            %             else
            %                 gamma(:,ng,b)=0;
            %             end
%             for q=1:Q
%                 fprintf('\t\t poids relatifs: gamma(:,ng=%d,b=%d)=%.2f\n',ng,b,gamma(q,ng,b))
%             end
            %             criteres_testes=squeeze(J(ng,b,theta_precedent));
            %             for q=1:Q
            %                 [~,qmax]=max(criteres_testes);
            %                 criteres_testes(qmax)=-Inf;
            %                 gamma(qmax,ng,b)= (Q+1-q)^2 ;
            %             end
            %             gamma(:,ng,b)=gamma(:,ng,b)/sum(gamma(:,ng,b));
        end
        1;
    end
    1;
end