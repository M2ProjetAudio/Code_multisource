clear all
close all
clc

load('y_stereo.txt','ascii');

y_sc=y_stereo/max(max(y_stereo));   %y scaled entre -1 1

Fe=24414;    % frequence echantillonnage
Te=1/Fe;
Ta=length(y_sc)/24414;  %temps d'acquisition

%echelle des temps :
t=0:Te:Ta-Te;





yg=y_sc(:,1);
yd=y_sc(:,2);

s1=yg(1:end/2);
s2=yg(end/2+1:end);
s1b=yd(1:end/2);
s2b=yd(end/2+1:end);

signal=zeros(length(yg),2);
signal(2:2:end,1)=s1;
signal(1:2:end,1)=s1b;
signal(2:2:end,2)=s2;
signal(1:2:end,2)=s2b;

signal(~(signal(:,1)~=0 .* signal(:,2)~=0),:)=[];

p=audioplayer(signal,Fe);

%play(p)
%audiowrite('tele_pierre.wav',y,Fe);

% figure des sons gauche et droite
figure(1)

subplot 211
plot(signal(:,1));grid on
title('micro gauche')
subplot 212; plot(signal(:,2))
title('micro droit')

grid on
