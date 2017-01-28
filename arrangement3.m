function Jp=arrangement3(J)

Jp=zeros(size(J));

for ng=1:size(J,1)
    Jint=squeeze(J(ng,:,:));
    Jp(ng,:,:)=arrangement(Jint,0.75);
end
% 
% 
% m=min(min(min(J)));
% M=max(max(max(J)));
% 
% 
% a=1/(M-m);
% b=-a*m;
% 
% Jp=a*J+b*ones(size(J));