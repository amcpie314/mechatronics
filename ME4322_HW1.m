%Andrew Cunningham & Dylan Furtado
%Homework 1

clc
clear

%FIRST POSITION

%joint coordinates
A=[1.4 0.485 0];
B=[1.67 0.99 0];
C=[0.255 1.035 0];
D=[0.285 0.055 0];
E=[0.195 2.54 0];
F=[-0.98 2.57 0];
G=[0.05 0.2 0];

%output load
Force=[0 -200 0];
posForce=[-1.67 4.279 0];

%position vectors
posAB=B-A;
posBC=C-B;
posCE=E-C;
posDE=E-D;
posDC=C-D;
posEF=F-E;
posFG=G-F;
posGWeight=posForce-G;

%masses of links (kg) (from solidworks: cast carbon steel data)
massAB=23.18905;
massBC=56.06995;
massDE=97.83505;
massEF=46.69825;
massFG=101.63755;

%mass moments of inertia of links (kg*m^2) (from solidworks: cast carbon steel data)
%format: [Ixx Iyy Izz]
inertiaAB=[0.03935190579 0.73841247893 0.73911596007];
inertiaBC=[0.09415340579 9.79129714774 9.79200062887];
inertiaDE=[0.16376190579 51.53939125957 51.54009474070];
inertiaEF=[0.07853390579 5.69060889107 5.69131237221];
inertiaFG=[0.17009940579 57.76627240925 57.76697589038];

%weights and weight position vectors of links
w1=[0 -9.81*massAB 0];
w2=[0 -9.81*massBC 0];
w3=[0 -9.81*massDE 0];
w4=[0 -9.81*massEF 0];
w5=[0 -9.81*massFG 0];
posw1=posAB/2;
posw2=posBC/2;
posw3=posDE/2;
posw4=posEF/2;
posw5=posFG/2; 

%static force vectors
syms Ax Ay Bx By Cx Cy Dx Dy Ex Ey Fx Fy Gx Gy T
FA=[Ax Ay 0];
FB=[Bx By 0];
FC=[Cx Cy 0];
FD=[Dx Dy 0];
FE=[Ex Ey 0];
FF=[Fx Fy 0];
FG=[Gx Gy 0];
Torque=[0 0 T];

%Eqn for Link AB
eqn1=FA+FB+w1==0;
eqn2=cross(posAB,FB)+cross(posw1,w1)+Torque==0;

%Eqn for Link BC
eqn3=-FB+FC+w2==0;
eqn4=cross(posBC,FC)+cross(posw2,w2)==0;

%Eqn for Link DE
eqn5=-FC+FD+FE+w3==0;
eqn6=cross(posCE,FE)+cross(-posDC,FD)+cross(posw3,w3)==0;

%Eqn for Link EF
eqn7=-FE+FF+w4==0;
eqn8=cross(posEF,FF)+cross(posw4,w4)==0;

%Eqn for Link FG
eqn9=-FF+FG+Force+w5==0;
eqn10=cross(posFG,-FF)+cross(posGWeight,Force)+cross(posw5,w5)==0;

%solve for static forces
solution = solve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10],[Ax,Ay,Bx,By,Cx,Cy,Dx,Dy,Ex,Ey,Fx,Fy,Gx,Gy,T]);
forceAx=double(solution.Ax)
forceAy=double(solution.Ay)
forceBx=double(solution.Bx)
forceBy=double(solution.By)
forceCx=double(solution.Cx)
forceCy=double(solution.Cy)
forceDx=double(solution.Dx)
forceDy=double(solution.Dy)
forceEx=double(solution.Ex)
forceEy=double(solution.Ey)
forceFx=double(solution.Fx)
forceFy=double(solution.Fy)
forceGx=double(solution.Gx)
forceGy=double(solution.Gy)
InputTorque=double(solution.T)

%linkage lengths
AB=norm(B-A);
BC=norm(C-B);
CD=norm(D-C);
CE=norm(E-C);
DE=norm(E-D);
EF=norm(F-E);
GF=norm(F-G);
WG=norm(posForce-G);
WF=norm(posForce-F);

%angular velocity calcs
syms wBC wDE wEF wFG
wAB=[0 0 1];
omegaBC=[0 0 wBC];
omegaDE=[0 0 wDE];
omegaEF=[0 0 wEF];
omegaFG=[0 0 wFG];
eqn1=cross(wAB,posAB)+cross(omegaBC,posBC)+cross(omegaDE,-posDC)==0;
eqn2=cross(omegaDE,posDE)+cross(omegaEF,posEF)+cross(omegaFG,posFG)==0;
solution=solve([eqn1,eqn2],[wBC,wDE,wEF,wFG]);
w_BC=double(solution.wBC)
w_DE=double(solution.wDE)
w_EF=double(solution.wEF)
w_FG=double(solution.wFG)

%linear velocity calcs
v_BA=cross(wAB, posAB)
v_CA=cross([0 0 w_BC], posBC)+v_BA
v_EA=cross([0 0 w_DE], posDE)+v_CA
v_FA=cross([0 0 w_EF], posEF)+v_EA
v_GA=cross([0 0 w_FG], posFG)+v_FA

%angular acceleration calcs
syms aBC aDE aEF aFG
aAB=[0 0 0];
alphaBC=[0 0 aBC];
alphaDE=[0 0 aDE];
alphaEF=[0 0 aEF];
alphaFG=[0 0 aFG];
omegaBC=[0 0 solution.wBC];
omegaDE=[0 0 solution.wDE];
omegaEF=[0 0 solution.wEF];
omegaFG=[0 0 solution.wFG];
eqn1=cross(aAB,posAB)+cross(wAB, cross(wAB,posAB))+cross(alphaBC,posBC)+cross(omegaBC, cross(omegaBC,posBC))+cross(alphaDE,-posDC)+cross(omegaDE, cross(omegaDE,-posDC))==0;
eqn2=cross(alphaEF,posEF)+cross(omegaEF, cross(omegaEF,posEF))+cross(alphaFG,posFG)+cross(omegaFG, cross(omegaFG,posFG))+cross(alphaDE,posDE)+cross(omegaDE, cross(omegaDE,posDE))==0;
solution=solve([eqn1, eqn2],[aBC, aDE, aEF, aFG]);
a_BC=double(solution.aBC)
a_DE=double(solution.aDE)
a_EF=double(solution.aEF)
a_FG=double(solution.aFG)

%linear acceleration calcs
la_BA=cross(aAB, posAB)
la_CA=cross([0 0 a_BC], posBC)+la_BA
la_EA=cross([0 0 a_DE], posDE)+la_CA
la_FA=cross([0 0 a_EF], posEF)+la_EA
la_GA=cross([0 0 a_FG], posFG)+la_FA

%Eqn for Link AB
eqn1=FA+FB+w1==massAB*la_BA;
eqn2=cross(posAB,FB)+cross(posw1,w1)+Torque==0;

%Eqn for Link BC
eqn3=-FB+FC+w2==massBC*la_CA;
eqn4=cross(posBC,FC)+cross(posw2,w2)==norm(inertiaBC)*[0 0 a_BC];

%Eqn for Link DE
eqn5=-FC+FD+FE+w3==massDE*la_EA;
eqn6=cross(posCE,FE)+cross(-posDC,FD)+cross(posw3,w3)==norm(inertiaDE)*[0 0 a_DE];

%Eqn for Link EF
eqn7=-FE+FF+w4==massEF*la_FA;
eqn8=cross(posEF,FF)+cross(posw4,w4)==norm(inertiaEF)*[0 0 a_EF];

%Eqn for Link FG
eqn9=-FF+FG+Force+w5==massFG*la_GA;
eqn10=cross(posFG,-FF)+cross(posGWeight,Force)+cross(posw5,w5)==norm(inertiaFG)*[0 0 a_FG];

%solve for dynamic forces
solution = solve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10],[Ax,Ay,Bx,By,Cx,Cy,Dx,Dy,Ex,Ey,Fx,Fy,Gx,Gy,T]);
dynamicAx=double(solution.Ax)
dynamicAy=double(solution.Ay)
dynamicBx=double(solution.Bx)
dynamicBy=double(solution.By)
dynamicCx=double(solution.Cx)
dynamicCy=double(solution.Cy)
dynamicDx=double(solution.Dx)
dynamicDy=double(solution.Dy)
dynamicEx=double(solution.Ex)
dynamicEy=double(solution.Ey)
dynamicFx=double(solution.Fx)
dynamicFy=double(solution.Fy)
dynamicGx=double(solution.Gx)
dynamicGy=double(solution.Gy)
dynamicTin=double(solution.T)

%ROTATIONAL ANAYSIS

%initial angle 
theta_AB=atan2(B(2)-A(2),B(1)-A(1));
if(theta_AB<0)
    angleAB_horizontal=2*pi+theta_AB; %adjusting the angle to be in the ccw direction from the hozizontal
else 
    angleAB_horizontal = theta_AB;
end

%allocate position storage vectors
newB_x=zeros(360,1);
newB_y=zeros(360,1);
newC_x=zeros(360,1);
newC_y=zeros(360,1);
newE_x=zeros(360,1);
newE_y=zeros(360,1);
newF_x=zeros(360,1);
newF_y=zeros(360,1);
newW_x=zeros(360,1);
newW_y=zeros(360,1);

%allocate static force storage vectors
new_A=zeros(360,1);
new_B=zeros(360,1);
new_C=zeros(360,1);
new_D=zeros(360,1);
new_E=zeros(360,1);
new_F=zeros(360,1);
new_G=zeros(360,1);
new_Input=zeros(360,1);

%allocate rotational velocity storage vectors
new_wBC=zeros(360,1);
new_wDE=zeros(360,1);
new_wEF=zeros(360,1);
new_wFG=zeros(360,1);

%allocate rotational acceleration storage vectors
new_aBC=zeros(360,1);
new_aDE=zeros(360,1);
new_aEF=zeros(360,1);
new_aFG=zeros(360,1);

%allocate linear velocity storage vectors
new_vBA=zeros(360,1);
new_vCA=zeros(360,1);
new_vEA=zeros(360,1);
new_vFA=zeros(360,1);
new_vGA=zeros(360,1);

%allocate linear acceleration storage vectors
new_laBA=zeros(360,1);
new_laCA=zeros(360,1);
new_laEA=zeros(360,1);
new_laFA=zeros(360,1);
new_laGA=zeros(360,1);

%allocate dynamic force/torque storage vectors
new_FA=zeros(360,1);
new_FB=zeros(360,1);
new_FC=zeros(360,1);
new_FD=zeros(360,1);
new_FE=zeros(360,1);
new_FF=zeros(360,1);
new_FG=zeros(360,1);
new_Tin=zeros(360,1);

for theta=0:1:359

    %Calculate B
    B_new = vpa(A+[AB*cos(angleAB_horizontal+deg2rad(theta)) AB*sin(angleAB_horizontal+deg2rad(theta)) 0]);
    
    %Calculate C
    [C_x,C_y]=circcirc(B_new(1),B_new(2),BC,D(1),D(2),CD);
    circIntersect_x = any(isnan(vpa(C_x))); %checking for Not-a-Number
    circIntersect_y = any(isnan(vpa(C_y)));
    C_1=[C_x(1) C_y(1) 0]; %adding the two solutions
    C_2=[C_x(2) C_y(2) 0];
    dist1 = norm(C_1-C);
    dist2 = norm(C_2-C);
    if(dist1<dist2) %checking which new C is closer
      C_new=vpa(C_1);
    else
      C_new=vpa(C_2);
    end

    %Calculating E
    [E_x,E_y]=circcirc(C_new(1),C_new(2),CE,D(1),D(2),DE);
    circIntersect_Ex = any(isnan(vpa(E_x)));
    circIntersect_Ey = any(isnan(vpa(E_y)));
    E_1=[E_x(1) E_y(1) 0];
    E_2=[E_x(2) E_y(2) 0];
    dist3 = norm(E_1-E);
    dist4 = norm(E_2-E);
    if circIntersect_Ex==0 && circIntersect_Ey==0
        if(dist3<dist4) %checking which new E is closer
            E_new=vpa(E_1);
        else
            E_new=vpa(E_2);
        end
    end

    %Calculating F
    [F_x,F_y]=circcirc(E_new(1),E_new(2),EF,G(1),G(2),GF);
    circIntersect_Fx = any(isnan(vpa(F_x)));
    circIntersect_Fy = any(isnan(vpa(F_y)));
    F_1=[F_x(1) F_y(1) 0];
    F_2=[F_x(2) F_y(2) 0];
    dist5 = norm(F_1-F);
    dist6 = norm(F_2-F);
    if circIntersect_Ex==0 && circIntersect_Ey==0
        if(dist5<dist6) %checking which new F is closer
            F_new=vpa(F_1);
        else
            F_new=vpa(F_2);
        end
    end

    %Calculating W
    [W_x,W_y]=circcirc(F_new(1),F_new(2),WF,G(1),G(2),WG);
    circIntersect_Wx = any(isnan(vpa(W_x)));
    circIntersect_Wy = any(isnan(vpa(W_y)));
    W_1=[W_x(1) W_y(1) 0];
    W_2=[W_x(2) W_y(2) 0];
    dist7 = norm(W_1-posForce);
    dist8 = norm(W_2-posForce);
    if circIntersect_Wx==0 && circIntersect_Wy==0
        if(dist7<dist8) %checking which new W is closer
            W_new=vpa(W_1);
        else
            W_new=vpa(W_2);
        end
    end

    %storing values
    newB_x(theta+1)=B_new(1);
    newB_y(theta+1)=B_new(2);
    newC_x(theta+1)=C_new(1);
    newC_y(theta+1)=C_new(2);
    newE_x(theta+1)=E_new(1);
    newE_y(theta+1)=E_new(2);
    newF_x(theta+1)=F_new(1);
    newF_y(theta+1)=F_new(2);
    newW_x(theta+1)=W_new(1);
    newW_y(theta+1)=W_new(2);
    B=B_new;
    C=C_new;
    E=E_new;
    F=F_new;
    posForce=W_new;

    %new position vectors
    posAB=B-A;
    posBC=C-B;
    posCE=E-C;
    posDE=E-D;
    posDC=C-D;
    posEF=F-E;
    posFG=F-G;
    posGWeight=posForce-G;
    posw1=posAB/2;
    posw2=posBC/2;
    posw3=posDE/2;
    posw4=posEF/2;
    posw5=posFG/2; 

    %Eqn for Link AB
    eqn1=FA+FB+w1==0;
    eqn2=cross(posAB,FB)+cross(posw1,w1)+Torque==0;
    
    %Eqn for Link BC
    eqn3=-FB+FC+w2==0;
    eqn4=cross(posBC,FC)+cross(posw2,w2)==0;
    
    %Eqn for Link DE
    eqn5=-FC+FD+FE+w3==0;
    eqn6=cross(posCE,FE)+cross(-posDC,FD)+cross(posw3,w3)==0;
    
    %Eqn for Link EF
    eqn7=-FE+FF+w4==0;
    eqn8=cross(posEF,FF)+cross(posw4,w4)==0;
    
    %Eqn for Link FG
    eqn9=-FF+FG+Force+w5==0;
    eqn10=cross(posFG,-FF)+cross(posGWeight,Force)+cross(posw5,w5)==0;
    
    %solve for static forces
    solution = solve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10],[Ax,Ay,Bx,By,Cx,Cy,Dx,Dy,Ex,Ey,Fx,Fy,Gx,Gy,T]);
    new_A(theta+1)=sqrt((double(solution.Ax))^2+(double(solution.Ay))^2);
    new_B(theta+1)=sqrt((double(solution.Bx))^2+(double(solution.By))^2);
    new_C(theta+1)=sqrt((double(solution.Cx))^2+(double(solution.Cy))^2);
    new_D(theta+1)=sqrt((double(solution.Dx))^2+(double(solution.Dy))^2);
    new_E(theta+1)=sqrt((double(solution.Ex))^2+(double(solution.Ey))^2);
    new_F(theta+1)=sqrt((double(solution.Fx))^2+(double(solution.Fy))^2);
    new_G(theta+1)=sqrt((double(solution.Gx))^2+(double(solution.Gy))^2);
    new_Input(theta+1)=double(solution.T);
    
    %angular velocity calcs
    wAB=[0 0 1];
    omegaBC=[0 0 wBC];
    omegaDE=[0 0 wDE];
    omegaEF=[0 0 wEF];
    omegaFG=[0 0 wFG];
    eqn1=cross(wAB,posAB)+cross(omegaBC,posBC)+cross(omegaDE,-posDC)==0;
    eqn2=cross(omegaDE,posDE)+cross(omegaEF,posEF)+cross(omegaFG,posFG)==0;
    solution=solve([eqn1,eqn2],[wBC,wDE,wEF,wFG]);
    new_wBC(theta+1)=double(solution.wBC);
    new_wDE(theta+1)=double(solution.wDE);
    new_wEF(theta+1)=double(solution.wEF);
    new_wFG(theta+1)=double(solution.wFG);
    
    %linear velocity calcs
    new_vBA(theta+1)=norm(cross(wAB, posAB));
    new_vCA(theta+1)=norm(cross([0 0 double(solution.wBC)], posBC)+cross(wAB, posAB));
    new_vEA(theta+1)=norm(cross([0 0 double(solution.wDE)], posDE)+cross([0 0 double(solution.wBC)], posBC)+cross(wAB, posAB));
    new_vFA(theta+1)=norm(cross([0 0 double(solution.wEF)], posEF)+cross([0 0 double(solution.wDE)], posDE)+cross([0 0 double(solution.wBC)], posBC)+cross(wAB, posAB));
    new_vGA(theta+1)=norm(cross([0 0 double(solution.wFG)], posFG)+cross([0 0 double(solution.wEF)], posEF)+cross([0 0 double(solution.wDE)], posDE)+cross([0 0 double(solution.wBC)], posBC)+cross(wAB, posAB));

    %angular acceleration calcs
    aAB=[0 0 0];
    alphaBC=[0 0 aBC];
    alphaDE=[0 0 aDE];
    alphaEF=[0 0 aEF];
    alphaFG=[0 0 aFG];
    omegaBC=[0 0 solution.wBC];
    omegaDE=[0 0 solution.wDE];
    omegaEF=[0 0 solution.wEF];
    omegaFG=[0 0 solution.wFG];
    eqn1=cross(aAB,posAB)+cross(wAB, cross(wAB,posAB))+cross(alphaBC,posBC)+cross(omegaBC, cross(omegaBC,posBC))+cross(alphaDE,-posDC)+cross(omegaDE, cross(omegaDE,-posDC))==0;
    eqn2=cross(alphaEF,posEF)+cross(omegaEF, cross(omegaEF,posEF))+cross(alphaFG,posFG)+cross(omegaFG, cross(omegaFG,posFG))+cross(alphaDE,posDE)+cross(omegaDE, cross(omegaDE,posDE))==0;
    solution=solve([eqn1, eqn2],[aBC, aDE, aEF, aFG]);
    new_aBC(theta+1)=double(solution.aBC);
    new_aDE(theta+1)=double(solution.aDE);
    new_aEF(theta+1)=double(solution.aEF);
    new_aFG(theta+1)=double(solution.aFG);

    %linear acceleration calcs
    laBA=cross(aAB, posAB);
    new_laBA(theta+1)=norm(laBA);
    laCA=cross([0 0 double(solution.aBC)], posBC)+cross(aAB, posAB);
    new_laCA(theta+1)=norm(laCA);
    laEA=cross([0 0 double(solution.aDE)], posDE)+cross([0 0 double(solution.aBC)], posBC)+cross(aAB, posAB);
    new_laEA(theta+1)=norm(laEA);
    laFA=cross([0 0 double(solution.aEF)], posEF)+cross([0 0 double(solution.aDE)], posDE)+cross([0 0 double(solution.aBC)], posBC)+cross(aAB, posAB);
    new_laFA(theta+1)=norm(laFA);
    laGA=cross([0 0 double(solution.aFG)], posFG)+cross([0 0 double(solution.aEF)], posEF)+cross([0 0 double(solution.aDE)], posDE)+cross([0 0 double(solution.aBC)], posBC)+cross(aAB, posAB);
    new_laGA(theta+1)=norm(laGA);

    %Eqn for Link AB
    eqn1=FA+FB+w1==massAB*laBA;
    eqn2=cross(posAB,FB)+cross(posw1,w1)+Torque==0;
    
    %Eqn for Link BC
    eqn3=-FB+FC+w2==massBC*laCA;
    eqn4=cross(posBC,FC)+cross(posw2,w2)==norm(inertiaBC)*[0 0 solution.aBC];
    
    %Eqn for Link DE
    eqn5=-FC+FD+FE+w3==massDE*laEA;
    eqn6=cross(posCE,FE)+cross(-posDC,FD)+cross(posw3,w3)==norm(inertiaDE)*[0 0 solution.aDE];
    
    %Eqn for Link EF
    eqn7=-FE+FF+w4==massEF*laFA;
    eqn8=cross(posEF,FF)+cross(posw4,w4)==norm(inertiaEF)*[0 0 solution.aEF];
    
    %Eqn for Link FG
    eqn9=-FF+FG+Force+w5==massFG*laGA;
    eqn10=cross(posFG,-FF)+cross(posGWeight,Force)+cross(posw5,w5)==norm(inertiaFG)*[0 0 solution.aFG];

    %solve for dynamic forces
    solution = solve([eqn1,eqn2,eqn3,eqn4,eqn5,eqn6,eqn7,eqn8,eqn9,eqn10],[Ax,Ay,Bx,By,Cx,Cy,Dx,Dy,Ex,Ey,Fx,Fy,Gx,Gy,T]);
    new_FA(theta+1)=sqrt((double(solution.Ax))^2+(double(solution.Ay))^2);
    new_FB(theta+1)=sqrt((double(solution.Bx))^2+(double(solution.By))^2);
    new_FC(theta+1)=sqrt((double(solution.Cx))^2+(double(solution.Cy))^2);
    new_FD(theta+1)=sqrt((double(solution.Dx))^2+(double(solution.Dy))^2);
    new_FE(theta+1)=sqrt((double(solution.Ex))^2+(double(solution.Ey))^2);
    new_FF(theta+1)=sqrt((double(solution.Fx))^2+(double(solution.Fy))^2);
    new_FG(theta+1)=sqrt((double(solution.Gx))^2+(double(solution.Gy))^2);
    new_Tin(theta+1)=double(solution.T);

end

figure
                 
    ax1= subplot(2,3,1);
    plot(newB_x,newB_y);
    title(ax1,'Joint B')
    ax2=  subplot(2,3,2);
    plot(newC_x,newC_y);
    title(ax2,'Joint C')
    ax3= subplot(2,3,3);
    plot(newE_x,newE_y);
    title(ax3,'Joint E')
    ax4=  subplot(2,3,4);
    plot(newF_x,newF_y);
    title(ax4,'Joint F')
    ax13=  subplot(2,3,5);
    plot(newW_x,newW_y);
    title(ax13,'Part')

figure

    ax5=  subplot(2,4,1);
    plot(1:360,new_A);
    title(ax5,'Static Force @ Joint A')
    ax6=  subplot(2,4,2);
    plot(1:360,new_B);
    title(ax6,'Static Force @ Joint B')
    ax7=  subplot(2,4,3);
    plot(1:360,new_C);
    title(ax7,'Static Force @ Joint C')
    ax8=  subplot(2,4,4);
    plot(1:360,new_D);
    title(ax8,'Static Force @ Joint D')
    ax9=  subplot(2,4,5);
    plot(1:360,new_E);
    title(ax9,'Static Force @ Joint E')
    ax10=  subplot(2,4,6);
    plot(1:360,new_F);
    title(ax10,'Static Force @ Joint F')
    ax11=  subplot(2,4,7);
    plot(1:360,new_G);
    title(ax11,'Static Force @ Joint G')
    ax12=  subplot(2,4,8);
    plot(1:360,new_Input);
    title(ax12,'Static Torque')

figure

    ax14=  subplot(2,2,1);
    plot(1:360,new_wBC);
    title(ax14,'Angular Velocity of Link BC')
    ax15=  subplot(2,2,2);
    plot(1:360,new_wDE);
    title(ax15,'Angular Velocity of Link DE')
    ax16=  subplot(2,2,3);
    plot(1:360,new_wEF);
    title(ax16,'Angular Velocity of Link EF')
    ax17=  subplot(2,2,4);
    plot(1:360,new_wFG);
    title(ax17,'Angular Velocity of Link FG')

figure

    ax18=  subplot(2,3,1);
    plot(1:360,new_vBA);
    title(ax18,'Linear Velocity of Joint B wrt A')
    ax19=  subplot(2,3,2);
    plot(1:360,new_vCA);
    title(ax19,'Linear Velocity of Joint C wrt A')
    ax20=  subplot(2,3,3);
    plot(1:360,new_vEA);
    title(ax20,'Linear Velocity of Joint E wrt A')
    ax21=  subplot(2,3,4);
    plot(1:360,new_vFA);
    title(ax21,'Linear Velocity of Joint F wrt A')
    ax38=  subplot(2,3,5);
    plot(1:360,new_vGA);
    title(ax38,'Linear Velocity of Joint G wrt A')

figure

    ax22=  subplot(2,2,1);
    plot(1:360,new_aBC);
    title(ax22,'Angular Acceleration of Link BC')
    ax23=  subplot(2,2,2);
    plot(1:360,new_aDE);
    title(ax23,'Angular Acceleration of Link DE')
    ax24=  subplot(2,2,3);
    plot(1:360,new_aEF);
    title(ax24,'Angular Acceleration of Link EF')
    ax25=  subplot(2,2,4);
    plot(1:360,new_aFG);
    title(ax25,'Angular Acceleration of Link FG')

figure

    ax26=  subplot(2,2,1);
    plot(1:360,new_laBA);
    title(ax26,'Linear Acceleration of Joint B wrt A')
    ax27=  subplot(2,2,2);
    plot(1:360,new_laCA);
    title(ax27,'Linear Acceleration of Joint C wrt A')
    ax28=  subplot(2,2,3);
    plot(1:360,new_laEA);
    title(ax28,'Linear Acceleration of Joint E wrt A')
    ax29=  subplot(2,2,4);
    plot(1:360,new_laFA);
    title(ax29,'Linear Acceleration of Joint F wrt A')
    ax39=  subplot(2,2,4);
    plot(1:360,new_laGA);
    title(ax39,'Linear Acceleration of Joint G wrt A')

figure

    ax30=  subplot(2,4,1);
    plot(1:360,new_FA);
    title(ax30,'Dynamic Force @ Joint A')
    ax31=  subplot(2,4,2);
    plot(1:360,new_FB);
    title(ax31,'Dynamic Force @ Joint B')
    ax32=  subplot(2,4,3);
    plot(1:360,new_FC);
    title(ax32,'Dynamic Force @ Joint C')
    ax33=  subplot(2,4,4);
    plot(1:360,new_FD);
    title(ax33,'Dynamic Force @ Joint D')
    ax34=  subplot(2,4,5);
    plot(1:360,new_FE);
    title(ax34,'Dynamic Force @ Joint E')
    ax35=  subplot(2,4,6);
    plot(1:360,new_FF);
    title(ax35,'Dynamic Force @ Joint F')
    ax36=  subplot(2,4,7);
    plot(1:360,new_FG);
    title(ax36,'Dynamic Force @ Joint G')
    ax37=  subplot(2,4,8);
    plot(1:360,new_Tin);
    title(ax37,'Dynamic Input Torque')
