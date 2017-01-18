clear variables
close all
clc

Lframe=512;
Ng=10;
Nf=4;
B=128;
Q=2;
fs=44100;
Nbfreq=B;
freqIndexes=round(linspace(1,round(Lframe) ,B));
az=(-180:179);
Ntheta=length(az);
%%
   y{1}=mean(audioread('chasseurs.wav'),2);
   y{2}=mean(audioread('police.wav'),2);
   
   
taille_min=min(length(y{1}),length(y{2}));

signal_tot=[y{1}(1:taille_min),y{2}(1:taille_min)] ;


duree_son=length(signal_tot)/fs;


[hrir,H,P,V]=get_hrtf(Lframe,B,az,Nbfreq,Ntheta,freqIndexes);
%left_ear=zeros(length(signal_tot),1);
%right_ear=zeros(length(signal_tot),1);
%signal_spa=zeros(length(signal_tot),2);
%Positions
Positions=[90,-90];

% spatialisation  35
impulseResponse = hrir.getImpulseResponses(90);
left_ear=conv(signal_tot(:,1),impulseResponse.left);
right_ear=conv(signal_tot(:,1),impulseResponse.right); 
signal_spa{1} =[left_ear(1:length(signal_tot)),right_ear(1:length(signal_tot)) ];
 
 % spatialisation  -85
impulseResponse = hrir.getImpulseResponses(-90);
left_ear=conv(signal_tot(:,2),impulseResponse.left);
right_ear=conv(signal_tot(:,2),impulseResponse.right); 
signal_spa{2} =[left_ear(1:length(signal_tot)),...
            right_ear(1:length(signal_tot)) ];
        
        
signal_spa=(signal_spa{1}+signal_spa{2})/Q;


%%  ajout du bruit au niveau de la reception

Psig=mean(mean(signal_spa,2).^2);
sigma=sqrt(Psig/10)/10;

% ajout d'un bruit Basse Frequence
bruit1=randn(size(signal_tot,1),1);
bruit2=randn(size(signal_tot,1),1);
[b,a]=butter(3,2000/(fs/2),'low'); % je mets fc a 2khz
bruit1=filter(b,a,bruit1);
bruit2=filter(b,a,bruit2);

signal_spa=signal_spa+sigma*[bruit1,bruit2];
signal_spa=rechelonner(signal_spa);

%%

p=audioplayer(signal_spa,fs);
%play(p)


%% Algorithme de localisation
% le signal sonore s'appelle signal_spa
Taille_1_groupe=2.5*Lframe;
Taille_1_algo=Ng*Taille_1_groupe;   % environ 0.3 secondes
Nb_Loca=floor(length(signal_spa)/Taille_1_algo);

%calcul du V


%allocations
lieu=zeros(Nb_Loca,Q);

x1t{Ng,Nf}=zeros(Lframe,1);
x2t{Ng,Nf}=zeros(Lframe,1);
X1{Ng,Nf}=zeros(Lframe,1);
X2{Ng,Nf}=zeros(Lframe,1);
Zint=zeros(Nf*2,1);
Zint2=zeros(Ng*Nf*2,1);
Z=zeros(2,1,B,Ng,Nf);
w=hanning(Lframe);
coef=sqrt(1/Lframe);

for num_exp=1:Nb_Loca
    % Calcul de Qn
    Qn=sigma^2*eye(2,2);
    % Data acquisition
    deb=Taille_1_algo*(num_exp-1);
    x1=signal_spa(deb+1:deb+Taille_1_algo,1); 
    x2=signal_spa(deb+1:deb+Taille_1_algo,2);
    for ng=1:Ng
        deb1=Taille_1_groupe*(ng-1);
        for nf=1:Nf % segments internes
            %fprintf('%d:%d\n',deb+Lframe*(i-1)*.75+1,deb+Lframe*(i-1)*.75+Lframe);
            x1t{ng,nf}(:,1)=x1(deb1+floor(Lframe/2)*(nf-1)+1:deb1+floor(Lframe/2)*(nf-1)+Lframe).*w;
            x2t{ng,nf}(:,1)=x2(deb1+floor(Lframe/2)*(nf-1)+1:deb1+floor(Lframe/2)*(nf-1)+Lframe).*w;
            %         % fft
            X1{ng,nf}=coef*fft(x1t{ng,nf});
            X2{ng,nf}=coef*fft(x2t{ng,nf});
        end
    end
    
    % sss
    for k=1:B
        for ng=1:Ng
            for nf=1:Nf
                Z(:,:,k,ng,nf)=[X1{ng,nf}(freqIndexes(k));...
                    X2{ng,nf}(freqIndexes(k))];
                % zint est un vecteur (8,1)
            end
        end
    end
    % Z est un vecteur (B*Ng*Nf*2,1)
    
    
    % Data conditionning
    J=algo1(Z,Qn,V,B,Ng,Nf,Ntheta);
    theta_init=ones(Q,1);
    % Localization
    theta_estimee=algo2(theta_init,J,Q,Ng,B);
    lieu(num_exp,:)=theta_estimee';
end