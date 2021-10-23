global wN Re Le Ldd Ldf Lff Lqq Lqq1 Lq1q1 Loo Ra Rf Rq1 theta0 phasE Vf ;

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
Tsim= 3.0 ;

% -------------- END OF DATA --------------

% derived parameters

wN=2*pi*fN ;
TB=1/(2*pi*fN) ;
Le=Xe ;
Ldd=Xd ;
Ldf=Xd-Xl ;
Lff=Ldf^2/(Xd-Xpd) ;
Tpdo=Tpdo_s/TB ;
Rf=Lff/Tpdo ;
Lqq=Xq ;
Lqq1=Xq-Xl ;
Lq1q1=Lqq1^2/(Xq-Xpq) ;
Tpqo=Tpqo_s/TB ;
Rq1=Lq1q1/Tpqo ;

% allocate space to store trajectory

nbsteps=ceil(Tsim/h) ;
time=zeros(1,nbsteps+1) ;
xhist=zeros(17,nbsteps+1) ;

% initialization of state variables

phasI=(P-j*Q)/V ;             % network voltages and currents
va=sqrt(2)*V ;
vb=sqrt(2)*V*cos(-2*pi/3) ;
vc=sqrt(2)*V*cos( 2*pi/3) ;
ia=sqrt(2)*abs(phasI)*cos(angle(phasI)) ;
ib=sqrt(2)*abs(phasI)*cos(angle(phasI)-2*pi/3) ;
ic=sqrt(2)*abs(phasI)*cos(angle(phasI)+2*pi/3) ;
phasE=V-(Re+j*Xe)*phasI ;

phi=angle(V+Ra*phasI+j*Xq*phasI);     % machine
theta0=phi+pi/2 ;

id=(sqrt(2)/3)*(cos(theta0)*ia+cos(theta0-2*pi/3)*ib+cos(theta0+2*pi/3)*ic) ;
iq=(sqrt(2)/3)*(sin(theta0)*ia+sin(theta0-2*pi/3)*ib+sin(theta0+2*pi/3)*ic) ;
io=0 ;
vd=(sqrt(2)/3)*(cos(theta0)*va+cos(theta0-2*pi/3)*vb+cos(theta0+2*pi/3)*vc) ;
vq=(sqrt(2)/3)*(sin(theta0)*va+sin(theta0-2*pi/3)*vb+sin(theta0+2*pi/3)*vc) ;
vo=0 ;

psid= vq+Ra*iq ;
psiq=-vd-Ra*id ;
psio=0 ;
ifd=(psid-Ldd*id)/Ldf ;
psif=Lff*ifd+Ldf*id ;
Vf=Rf*ifd ;
iq1=0 ;
psiq1=Lqq1*iq ;

% identify differential and algebraic states
% note: the i-th differential equation is  d x(i)/dt = ...

dif= [1 2 3 15 16 17 18 19] ;
alg= [4 5 6 7 8 9 10 11 12 13 14] ;

% initialisation of integration

time(1)=0 ;
x=[ia ib ic va vb vc vd vq vo id iq io ifd iq1 psid psiq psif psiq1 psio]' ; 
xhist(1:19,1)=x ;

[A,u]=linmodel(0) ;

% integration

for k=1:nbsteps
    c= A*x + u ;  
    
    t=k*h ;
    time(k+1)=t ;
    
    [A,u]=linmodel(t) ;
    
    b= -u ;
    b(dif)= b(dif)- (2/h/wN)*x(dif) - c(dif) ;
    
    M= A;
    M(dif,dif)= M(dif,dif) - (2/h/wN)*eye(size(dif,2)) ;
    x= M\b ;
    xhist(1:19,k+1)=x ;
    
end
