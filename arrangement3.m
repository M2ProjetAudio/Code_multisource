function Jp=arrangement3(J)
Nth=size(J,3);

Jp=zeros(size(J));

for ng=1:size(J,1)
    for kb=1:size(J,2)
        Jint=squeeze(J(ng,kb,:));
        m=min(min(min(Jint)));
        M=max(max(max(Jint)));
        
        a=1/(M-m);
        b=-a*m;
        Jint=a*Jint+b*ones(size(Jint));
        
        %tout est entre 0 et 1
        %maintenant, je veux equi-repartir les data
        [~,idx]=sort(Jint);
        
        for th=1:Nth
            Jp(ng,kb,idx(th))=th/Nth;
        end
                
    end
end



