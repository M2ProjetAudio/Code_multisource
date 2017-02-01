clear variables
close all
clc

Lframe=1024;
Ng=10;
Nf=4;
B=128;
Q=2;      Positions=[85,-50,-175];
fs=44100;
Nbfreq=B;
freqIndexes=round(linspace(1,round(Lframe) ,B));
az=(-180:179);
% grille des azimuths testes
thetaArg=az*pi/180;

Ntheta=length(az);
%%
%    y{1}=mean(audioread('chasseurs.wav'),2);
%    y{2}=mean(audioread('police.wav'),2);
%    y{3}=mean(audioread('philo.wav'),2);

y{1}=randn(5*44100,1);
band1=[500/(fs/2) 1500/(fs/2)];
[b,a]=butter(4,band1);
y{1}=filter(b,a,y{1});

y{2}=randn(5*44100,1);
band2=[1500/(fs/2) 2500/(fs/2)];
[b,a]=butter(4,band2);
y{2}=filter(b,a,y{2});

y{3}=randn(5*44100,1);
band3=[3500/(fs/2) 4000/(fs/2)];
[b,a]=butter(4,band3);
y{3}=filter(b,a,y{3});

taille_min=min([length(y{1}),length(y{2}),length(y{3})]);

signal_tot=[y{1}(1:taille_min),y{2}(1:taille_min),y{3}(1:taille_min)] ;


duree_son=length(signal_tot)/fs;


[hrir,H,P,V]=get_hrtf(Lframe,B,az,Nbfreq,Ntheta,freqIndexes);
signal_spa{Q}=[];
for q=1:Q
    signal_spa{q}=zeros(taille_min,1);
end
for q=1:Q
    impulseResponse = hrir.getImpulseResponses(Positions(q));
    left_ear=conv(signal_tot(:,q),impulseResponse.left);
    right_ear=conv(signal_tot(:,q),impulseResponse.right);
    signal_spa{q} =[left_ear(1:length(signal_tot)),right_ear(1:length(signal_tot)) ];
    
end
%signal_spa=(signal_spa{1}+signal_spa{2})/Q;
sum=zeros(size(signal_spa{1}));
for q=1:Q
    sum=sum+signal_spa{q};
end
clear 'signal_spa';
signal_spa=sum/Q;
%%  ajout du bruit au niveau de la reception

Psig=mean(mean(signal_spa,2).^2);
sigma=sqrt(Psig/10)/10;

% ajout d'un bruit Basse Frequence
bruit1=sigma*randn(size(signal_tot,1),1);
bruit2=sigma*randn(size(signal_tot,1),1);
[b,a]=butter(3,2000/(fs/2),'low'); % je mets fc a 2khz
bruit1=filter(b,a,bruit1);
bruit2=filter(b,a,bruit2);

signal_spa=signal_spa+sigma*[bruit1,bruit2];
%signal_spa=rechelonner(signal_spa);
%sigma=1;
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
    Qn=chol(sigma^2*eye(2,2));
    % Data acquisition
    deb=Taille_1_algo*(num_exp-1);   %
    % deb=Taille_1_algo*(num_exp+round(Nb_Loca/2)-1);
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
    J=real(algo1(Z,1,V,B,Ng,Nf,Ntheta));
    J=arrangement3(J);
    SQ=round(linspace(1,359,Q+1));
    theta_init=SQ(1:end-1)';
    theta_init=[50,240];
    % Localization
    theta_estimee=algo2(theta_init,J,Q,Ng,B,thetaArg);
    lieu(num_exp,:)=thetaArg(theta_estimee)*180/pi;
end