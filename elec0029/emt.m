global wN Re Le Ldd Ldf Lff Lqq Lqq1 Lq1q1 Loo Ra Rf Rq1 thetar0 ;

j=sqrt(-1) ;

% --------------- DATA ---------------
% network and operating point
% P and Q are in per unit on a single-phase base power (three-phase model)

Xe=0.20 ;
Re=0.01 ;
P=0.5 ;
Q=0.1 ;
V=1 ;

% machine

Ra=0.005 ;
Xd=2.4 ;
Xq=2.4 ;
Xl=0.2 ;
Xpd=0.40 ;
Xpq=0.25 ;
Loo=0.1 ;
fN=50. ;
Tpdo_s=7. ;
Tpqo_s=0.3 ;

% time step size and simulation interval

h= 100e-06 ;
Tsim= 3 ;

% -------------- END OF DATA --------------

% derive parameters from data

wN= 2*pi*fN ;
TB= 1/(2*pi*fN) ;
Le= Xe ;
Ldd= Xd ;
Ldf= Xd-Xl ;
Lff= Ldf^2/(Xd-Xpd) ;
Tpdo= Tpdo_s/TB ;
Rf= Lff/Tpdo ;
Lqq= Xq ;
Lqq1= Xq-Xl ;
Lq1q1= Lqq1^2/(Xq-Xpq) ;
Tpqo= Tpqo_s/TB ;
Rq1= Lq1q1/Tpqo ;

% initial value of state variables

phasI= (P-j*Q)/V ;             % network voltages and currents
va= sqrt(2)*V ;
vb= sqrt(2)*V*cos(-2*pi/3) ;
vc= sqrt(2)*V*cos( 2*pi/3) ;
ia= sqrt(2)*abs(phasI)*cos(angle(phasI)) ;
ib= sqrt(2)*abs(phasI)*cos(angle(phasI)-2*pi/3) ;
ic= sqrt(2)*abs(phasI)*cos(angle(phasI)+2*pi/3) ;
phasE= V-(Re+j*Xe)*phasI ;

phi= angle(V+Ra*phasI+j*Xq*phasI);     % machine
thetar0= phi+pi/2 ;

id= sqrt(2/3)*(cos(thetar0)*ia+cos(thetar0-2*pi/3)*ib+cos(thetar0+2*pi/3)*ic) ;
iq= sqrt(2/3)*(sin(thetar0)*ia+sin(thetar0-2*pi/3)*ib+sin(thetar0+2*pi/3)*ic) ;
io= 0 ;
vd= sqrt(2/3)*(cos(thetar0)*va+cos(thetar0-2*pi/3)*vb+cos(thetar0+2*pi/3)*vc) ;
vq= sqrt(2/3)*(sin(thetar0)*va+sin(thetar0-2*pi/3)*vb+sin(thetar0+2*pi/3)*vc) ;
vo= 0 ;

psid=  vq+Ra*iq ;
psiq= -vd-Ra*id ;
psio= 0 ;
ifd= (psid-Ldd*id)/Ldf ;
psif= Lff*ifd+Ldf*id ;
Vf= Rf*ifd ;
iq1= 0 ;
psiq1= Lqq1*iq ;

% identify differential and algebraic states
% note: the i-th differential equation is  d x(i)/dt = ...

dif= [1 2 3 15 16 17 18 19] ;
alg= [4 5 6 7 8 9 10 11 12 13 14] ;

% allocate space to store trajectory

nbsteps= ceil(Tsim/h) ;
time= zeros(1,nbsteps+1) ;
xhist= zeros(19,nbsteps+1) ;

% store initial values of time and state vector

time(1)= 0 ;
x= [ia ib ic va vb vc vd vq vo id iq io ifd iq1 psid psiq psif psiq1 psio]' ; 
xhist(1:19,1)= x ;

% compute initial value of u and of Ax+u

u= zeros(19,1) ; 
u(1) = -sqrt(2)*abs(phasE)*cos(angle(phasE))/Le ;
u(2) = -sqrt(2)*abs(phasE)*cos(angle(phasE)-2*pi/3)/Le ;
u(3) = -sqrt(2)*abs(phasE)*cos(angle(phasE)+2*pi/3)/Le ; 
u(17)= Vf ;

[A]= matrixA(0) ;
c= A*x+u ;
 
% integration by Trapezoidal maethod

for k=1:nbsteps
    
    t=k*h ;
    time(k+1)=t ;
    
    if t<= 0.05
        ea= sqrt(2)*abs(phasE)*cos(wN*t+angle(phasE)) ;
        eb= sqrt(2)*abs(phasE)*cos(wN*t+angle(phasE)-2*pi/3) ;
        ec= sqrt(2)*abs(phasE)*cos(wN*t+angle(phasE)+2*pi/3) ;
    else
        ea=0. ;
        eb=0. ;
        ec=0. ;
    end
    u(1)= -ea/Le ;
    u(2)= -eb/Le ;
    u(3)= -ec/Le ;

    [A]= matrixA(t) ;
    M= A;
    M(dif,dif)= M(dif,dif) - (2/h/wN)*eye(size(dif,2)) ;
    
    b= -u ;              
    b(dif)= b(dif)- (2/h/wN)*x(dif) - c(dif) ; 
    
    x= M\b ;
    
    xhist(1:19,k+1)=x ;
    
    c= A*x+u ;
end
