function ak = B(x,l,al)
syms w xx
h=x(2:end)-x(1:end-1);
n = length(x)-1; % number of subintervals

M = zeros(n+1,n+1); % allocate mass matrix
 B1 = zeros(n-1,n-1);B2= zeros(n-1,n-1);
c_a=-1/(2*cos(0.5*al*pi));
 k=-cos(al*pi*0.5)*l^al;
 teta=atan(abs(w)/l);
%%%%%%%%%%%%%%%Fourier transform by heaviside
    for i=2:n
   hs = (heaviside(xx-x(i-1))-heaviside(xx-x(i)))*(xx-x(i-1))/h(i-1)+ (heaviside(xx-x(i))-heaviside(xx-x(i+1)))*(x(i+1)-xx)/h(i);
   uhat(i)=fourier(hs,xx,w);
   vhat(i)=conj(uhat(i));
    end
    hs = (heaviside(xx-x(1))-heaviside(xx-x(2)))*2*(x(2)-xx)/h(1);
    uhat(1)=fourier(hs,xx,w);
   vhat(1)=conj(uhat(1));
   hs = (heaviside(xx-x(n))-heaviside(xx-x(n+1)))*2*(xx-x(n))/h(n);
    uhat(n+1)=fourier(hs,xx,w);
   vhat(n+1)=conj(uhat(n+1));

%%%%%%%%%%%%%%%Fourier Transform Definition 
%   for i=2:n
%    uhat(i)=int((xx-x(i-1))/h(i-1)*exp(-sqrt(-1)*w*xx),xx,x(i-1),x(i))+int((x(i+1)-xx)/h(i)*exp(-sqrt(-1)*w*xx),xx,x(i),x(i+1));
%    vhat(i)=conj(uhat(i));
%     end
%     uhat(1)=int((x(2)-xx)/h(1)*exp(-sqrt(-1)*w*xx),xx,x(1),x(2));
%    vhat(1)=conj(uhat(1)); 
% 
%     uhat(n+1)=int((xx-x(n))/h(n)*exp(-sqrt(-1)*w*xx),xx,x(n),x(n+1));
%    vhat(n+1)=conj(uhat(n+1));
%
 s=w/(1-w^2);
for i = 1:n+1 % loop over subintervals
    for j=1:n+1
     fii=2*(l^2+w^2)^(0.5*al)*cos(al*teta)*uhat(i)*vhat(j);  
     fii=diff(s)/(2*pi)*subs(fii,w,s);
 B1(i,j) = lgwt(fii,-1,1,n);
    end
end
 %%%%%%%%%%%%%%%%%%%
%   s=w/(1-w^2);
%   %s=(2*w-1)/(w*(1-w));
% for d = 1:n-1 % loop over subintervals
%   
%          fii=2*(l^2+w^2)^(0.5*al)*cos(al*teta)*uhat(d)*vhat(1);  
%      fii=diff(s)/(2*pi)*subs(fii,w,s);
%      b1(d) = lgwt(fii,-1,1);
% end
%     for i=1:n-1
%         m=i-1;
%        for j=1:n-1-m
%  B1(j,j+m) =b1(i);
%  B1(j+m,j) =b1(i);
%        end
%     end


for i = 1:n % loop over subintervals
h = x(i+1) - x(i); % interval length
M(i,i) = M(i,i) + h/3; % add h/3 to M(i,i)
M(i,i+1) = M(i,i+1) + h/6;
M(i+1,i) = M(i+1,i) + h/6;
M(i+1,i+1) = M(i+1,i+1) + h/3;
end

ak=-c_a*(B1)+((2*c_a*l^al+k)*M);
 end
