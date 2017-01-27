function Jp=arrangement3(J)

m=min(min(min(J)));
M=max(max(max(J)));


a=1/(M-m);
b=-a*m;

Jp=a*J+b*ones(size(J));