function dr=DR(f,a,b,al,l)

   syms w xx
   fw=subs(f,xx,w);
  
%%%%%%%%%%%%%%%%%
fr=(w-xx)^(1-al)*exp(-l*w)*fw;
 intfr=lgwtd(fr,xx,b);
   dfr=diff(intfr,'xx',2);
dr=exp(l*xx)/gamma(2-al)*dfr;
%%%%%%%%%%%%%%%%%%%%%%%


