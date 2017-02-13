function Jpos=arrangement(J)


m=min(J(:));
M=max(J(:));

a=1/(M-m);
b=-a*m;
Jpos=a*J+b*ones(size(J));

