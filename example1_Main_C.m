clear ; close all;
clc;
n=4;
m=4;
l=1;
al=1.1;
x0=0;xf=1;
t0=0;tf=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5555
fid = fopen('results.txt', 'w');
fprintf(fid," *************************************************************\n");
fprintf(fid," *************************************************************\n");
fprintf(fid,"                                                               \n");
fprintf(fid,"                                                               \n");
fprintf(fid,"            Absolute Error       L_2            Elapsed time \n");
fprintf(fid," --------------------------------------------------------------------------\n");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for q=1:2
    clear F f fm c0 C  beta uh u_exact x t Ax Bx E u0 t1 x1 Z
    N(q)=n;
to=(tf-t0)/m;
t=t0:to:tf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hh=(xf-x0)/n;
x=x0:hh:xf;
x=x';
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
h=x(2:end)-x(1:end-1);
[t1,x1]=meshgrid(t0:to:tf,x0:hh:xf);

syms tt xx 

  

  
 u0=-2*l^(6-al)*cos(pi*al*0.5)*gamma(-al)*x.^3.*(1-x).^3;
 c0=u0; 
 %c0=u0*hh;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %u_exact%%%%%%%%%%%%%%%%%%%%%%%%%
 for j=1:m+1
    for i=1:n+1
u_exact(i,j)=(-2*l^(6-al)*cos(pi*al*0.5)*gamma(-al)*exp(-t(j))*x(i)^3*(1-x(i))^3);
    end
 end
 %%%%%beta%%%%%%%%%%%%%%%%%%%%%%%
    for i=2:n
        beta(i) = piecewise(xx<x(i-1), 0,x(i-1)<=xx<=x(i), (xx-x(i-1))/h(i-1) ,x(i)<xx<=x(i+1), (x(i+1)-xx)/h(i) ,xx>x(i+1), 0);
    end
    beta(n+1)=piecewise(xx<x(n), 0,x(n)<=xx<=x(n+1), (xx-x(n))/h(n));
    beta(1)=piecewise(x(1)<=xx<=x(2), (x(2)-xx)/h(1) ,xx>x(2), 0);
    beta=subs(beta,xx,x);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ue=-2*l^(6-al)*cos(pi*al*0.5)*gamma(-al)*exp(-tt)*xx^3*(1-xx)^3;
    ut=diff(ue,tt);
    Dx=-1/(2*cos(0.5*al*pi))*(DR(xx^3*(1-xx)^3,x0,xf,al,l)+DL(xx^3*(1-xx)^3,x0,xf,al,l)-2*l^al*xx^3*(1-xx)^3);
ux=-2*l^(6-al)*cos(pi*al*0.5)*gamma(-al)*exp(-tt)*Dx;
f=ut-ux;
f=subs(f,xx,x(2:n)');
f=[0,f,0];
%f=subs(f,xx,[0.0000001,x(2:n)',0.9999999]);

C = zeros(n+1,m+1);
C(:,1)=c0;
 k=-cos(al*pi*0.5)*l^al;
Ax=A(x);Bx=real(B(x,l,al));Bx=double(Bx);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
for j=2:m+1
    fm=0.5*(subs(exp(-k*t(j))*f,tt,t(j))+subs(exp(-k*t(j-1))*f,tt,t(j-1)));
    fm=double(fm); 
F(:,j)=fm'.*hh;

    C(:,j)=(Ax+0.5*to*Bx)\((Ax-0.5*to*Bx)*C(:,j-1)+to*F(:,j));
end
Elapsedtime1(q)=toc;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

   
for j=1:m+1

    uh(:,j)=exp(k*t(j))*C(:,j)'*beta;
end
uh=double(uh);
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%
E=(uh-u_exact);

absoluteerror1(q)=max(max(abs(E)));
for i=1:m+1
    s(i)=(E(:,i)'*E(:,i));
end
L21(q)=((tf-t0)/m*(xf-x0)/n*(sum(s)))^0.5;

figure
fs=subs(f',tt,t);
fs=double(fs);
mesh(t1, x1,fs);
xlabel('\bf t');ylabel('\bf x');zlabel('\bf  Source Function');
title([' Source Function for n=',num2str(N(q))])


figure
mesh(t1, x1,real(uh));
xlabel('\bf t');ylabel('\bf x');zlabel('\bf u(t,x)');
title('u_{Approximated} for Crank Nicolson')

figure
 mesh(t1, x1,u_exact)
title('u_{exact} for Crank Nicolson');
xlabel('\bf t');ylabel('\bf x');zlabel('\bf u(t,x)');
figure
mesh(t1, x1,(abs(E)))
title('absolute error for Crank Nicolson');
xlabel('\bf t');ylabel('\bf x');zlabel('\bf absolute error');
uh1=uh; 




%%%%%%%%%%%%%
%%%%%%%%%%ETD
clear C  fm F
C = zeros(n+1,m+1);
C(:,1)=c0;
D=(Ax)\Bx;
tic
for j=1:m+1
    fm=subs(exp(-k*t(j))*f,tt,t(j));
    F(:,j)=double(fm');
end
for j=1:m
     C(:,j+1)=(eye(n+1)+1/3*D*to)\(9*C(:,j)+2*to*F(:,j)+to*(F(:,j+1)))+...
        (eye(n+1)+1/4*D*to)\(-8*C(:,j)-1.5*to*F(:,j)-0.5*to*(F(:,j+1)));

end
Elapsedtime3(q)=toc;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
   
for j=1:m+1

    uh(:,j)=exp(k*t(j))*C(:,j)'*beta;
end
uh=double(uh);
%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%
E=(uh-u_exact);

absoluteerror3(q)=max(max(abs(E)));
for i=1:m+1
    s(i)=(E(:,i)'*E(:,i));
end
L23(q)=((tf-t0)/m*(xf-x0)/n*(sum(s)))^0.5;

figure
mesh(t1, x1,real(uh));
xlabel('\bf t');ylabel('\bf x');zlabel('\bf u(t,x)');
title('u_{Approximated} for ETD')


figure
mesh(t1, x1,(abs(E)))
title('absolute error for ETD');
xlabel('\bf t');ylabel('\bf x');zlabel('\bf absolute error');
 uh3=uh;
figure
j=floor(m/2)+1;
plot(x,u_exact(:,j),'-.');
hold on
plot(x,real(uh1(:,j)),'-*')
hold on

plot(x,real(uh3(:,j)),'-+')
xlabel('\bf x'); ylabel('\bf u(x)');
 legend('u_{Exact}','u_{Approximated} for Crank','u_{Approximated} for ETD')
title(['u_{Approximated} and u_{Exact} for t=',num2str(t(j)),' and n=',num2str(N(q))]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,"N=%.0f\n",N(q));
fprintf(fid," --------------------------------------------------------------------------\n");
fprintf(fid,"Crank       %f             %f       %f \n",absoluteerror1(q),L21(q),Elapsedtime1(q));
fprintf(fid,"ETD         %f             %f       %f \n",absoluteerror3(q),L23(q),Elapsedtime3(q));
fprintf(fid," --------------------------------------------------------------------------\n");

Ep1(q)=L21(q);
Ep3(q)=L23(q);
n=2*n;
m=2*m;
end
fprintf(fid," *************************************************************\n");
fprintf(fid," *************************************************************\n");
fprintf(fid,"                                                               \n");
fprintf(fid," N   &  \\\\       p_Crank     &  \\\\        p_ETD  \\\\   hline\\\\ \n");
fprintf(fid," --------------------------------------------------------------------------\n");
for q=1:size(N,2)-1 
   o1(q)=log(Ep1(q)/Ep1(q+1))/log(2);
   o3(q)=log(Ep3(q)/Ep3(q+1))/log(2);
fprintf(fid," %.0f   &           %f        &         %f  \\\\  \n",N(q),o1(q),o3(q));
end
fclose(fid);
%%%%%%%%%%%%%%%%%%%%COMPARISON FIGURES
figure
plot(N,Elapsedtime1,'-*')
hold on
plot(N,Elapsedtime3,'-o')
xlabel('\bf n');ylabel('\bf Elapsedtime');
legend('CRANK','ETD')



figure
plot(N,absoluteerror1,'-*')
hold on
plot(N,absoluteerror3,'-o')
xlabel('\bf n');ylabel('\bf absolute error');
legend('CRANK','ETD')


figure
plot(N,L21,'-*')
hold on
plot(N,L23,'-o')
xlabel('\bf n');ylabel('\bf L2');
legend('CRANK','ETD')

figure
plot(N(1:end-1),o1,'-*')
hold on
plot(N(1:end-1),o3,'-o')
xlabel('\bf n');ylabel('\bf Order of accuracy');
legend('CRANK','ETD')