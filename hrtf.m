clear all
close all
clc

c = 340 ; %c�l�rit� du son en m�tre par seconde
a = 0.145/2 ; %rayon de la sph�re rigide en m�tre
r = 3; %distance � la source

threshold = 3;

Hg = zeros(360,512);
Hd = zeros(360,512);

for theta=1:360
    
    for f=1:512
        
        Hg(theta,f) = oreille_gauche(a,r,theta,f,c,threshold);
        
        Hd(theta,f) = oreille_droite(a,r,theta,f,c,threshold);
    end
end

%%%% Calcul de la HRIR � partir de la HRTF
HRIR_1 = ifftshift(Hg(90,:)); % permet d'enlever le surplus de valeur et de faire la transform� inverse

HRIR_gauche = (1/512)*real(HRIR_1); % on ne garde que la partie r�el

%%%% Calcul de la HRIR � partir de la HRTF
HRIR_2 = ifftshift(Hd(90,:)); % permet d'enlever le surplus de valeur et de faire la transform� inverse

HRIR_droite = (1/512)*real(HRIR_2); % on ne garde que la partie r�el

%%%% Cr�ation du son
x = randn(1,44100);

y1 = conv(x,HRIR_gauche); % on cr�er le son tel qu'il provient � 90 degr�s
y_gauche = y1(1,1:size(x,2));

y2 = conv(x,HRIR_droite); % on cr�er le son tel qu'il provient � 90 degr�s
y_droite = y2(1,1:size(x,2));

output_signal = [y_gauche;y_droite];

output_signal = rechelonner(output_signal');

p=audioplayer(output_signal',44100);
p1=audioplayer(x',44100);

%%% max(max(abs(HRIR_gauche-HRIR_droite))) = 2.03


%%% Deuxi�me m�thode pour HRIR
X1 = fft(x,512) .* Hg(90,:);
X2 = fft(x,512) .* Hd(90,:);

m1 = real(ifft(X1,44100));
m2 = real(ifft(X2,44100));

output_signal1 = [m1 ; m2];
output_signal1 = rechelonner(output_signal1');

p2=audioplayer(output_signal1',44100);


