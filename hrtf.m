clear all
close all
clc

c = 340 ; %célérité du son en mètre par seconde
a = 0.145/2 ; %rayon de la sphère rigide en mètre
r = 3; %distance à la source

threshold = 3;

Hg = zeros(360,512);
Hd = zeros(360,512);

for theta=1:360
    
    for f=1:512
        
        Hg(theta,f) = oreille_gauche(a,r,theta,f,c,threshold);
        
        Hd(theta,f) = oreille_droite(a,r,theta,f,c,threshold);
    end
end

%%%% Calcul de la HRIR à partir de la HRTF
HRIR_1 = ifftshift(Hg(90,:)); % permet d'enlever le surplus de valeur et de faire la transformé inverse

HRIR_gauche = (1/512)*real(HRIR_1); % on ne garde que la partie réel

%%%% Calcul de la HRIR à partir de la HRTF
HRIR_2 = ifftshift(Hd(90,:)); % permet d'enlever le surplus de valeur et de faire la transformé inverse

HRIR_droite = (1/512)*real(HRIR_2); % on ne garde que la partie réel

%%%% Création du son
x = randn(1,44100);

y1 = conv(x,HRIR_gauche); % on créer le son tel qu'il provient à 90 degrès
y_gauche = y1(1,1:size(x,2));

y2 = conv(x,HRIR_droite); % on créer le son tel qu'il provient à 90 degrès
y_droite = y2(1,1:size(x,2));

output_signal = [y_gauche;y_droite];

output_signal = rechelonner(output_signal');

p=audioplayer(output_signal',44100);
p1=audioplayer(x',44100);

%%% max(max(abs(HRIR_gauche-HRIR_droite))) = 2.03


%%% Deuxième méthode pour HRIR
X1 = fft(x,512) .* Hg(90,:);
X2 = fft(x,512) .* Hd(90,:);

m1 = real(ifft(X1,44100));
m2 = real(ifft(X2,44100));

output_signal1 = [m1 ; m2];
output_signal1 = rechelonner(output_signal1');

p2=audioplayer(output_signal1',44100);


