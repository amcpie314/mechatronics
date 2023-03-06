%Andrew Cunnnigham, Dylan Furtado, Bailey Koestner
%ME/RBE 4322: Modeling and Analysis of Mechatronic Systems
%Project

clc 
clear

s0=[0 0 0 0 0 0];               %Initial conditions
tspan=[0,0.5];                  %Window of time to solve within
[t,x]=ode45(@states,tspan,s0);

figure

    ax1=  subplot(2,3,1);
    plot(x(:,1));
    title(ax1,'lambda1(t)')

    ax2=  subplot(2,3,2);
    plot(x(:,2));
    title(ax2,'lambda2(t)')

    ax3=  subplot(2,3,3);
    plot(x(:,3));
    title(ax3,'x2(t)')

    ax4=  subplot(2,3,4);
    plot(x(:,4));
    title(ax4,'p(t)')

    ax5=  subplot(2,3,5);
    plot(x(:,5));
    title(ax5,'h1(t)')

    ax6=  subplot(2,3,6);
    plot(x(:,6));
    title(ax6,'h2(t)')

function sprime = states(t,x)
%Define Constants
L = 40;
L2 = 40;
J1 = 0.000034;
J2 = 0.000086;
J3 = 0.000127;
MD = 0.2;
R = 140;
R2 = 140;
B = 0.000000128;
B1 = 0.000000128;
BBelt = 0.05;
RB = 0.1;
Fext = 0;
CB = 0.1;
E = 12;
n = 10;
Gi = 0.01255;
Go = 0.00325;
Rc = 0.019;

%Define ODEs
s1 = E-R*x(1)/L-x(5)/(n*J1);
s2 = E-R2*x(2)/L2-x(6)/(n*J3);
s3 = Gi^2/Go*x(5)/J1-x(4)/MD+Rc*x(6)/J2-CB/RB*x(3);
s5 = (x(1)*n/L-x(5)*B/J1+x(5)*BBelt/J1+CB*x(3)*Gi)/(1+1/J1);
s6 = (x(3)*n/L2-Rc*x(4)-B1*n*E)/(1+Rc^2*MD/J2);
s4 = (x(2)*n/L2-s6-B1*x(6)/J2)/Rc+Fext;
sprime=[s1; s2; s3; s4; s5; s6];
end