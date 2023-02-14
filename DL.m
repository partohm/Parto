function dl=DL(f,a,b,al,l)

   syms w xx
   fw=subs(f,xx,w);
   f=(xx-w)^(1-al)*exp(l*w)*fw;
   intf=lgwtd(f,a,xx);
   df=diff(intf,'xx',2);
dl=exp(-l*xx)/gamma(2-al)*df;
%%%%%%%%%%%%%%%%%



