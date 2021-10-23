% ONE-MACHINE INFINITE-BUS SYSTEM WITHOUT PSS

% ----- DATA -----

% network 

Xe=0.45 ;

% machine

Xd=2.2 ;
Xq=2.2 ;
Xl=0.2 ;
Xpd=0.3 ;
Xpq=0.25 ;
Tpdo_s=7.0 ;
Tpqo_s=0.4 ;
H=4.0 ;
G=70. ;
T=0.4 ;
fN=50. ;

% operating point

V=1.0 ;
P=0.7 ;
Q=0.15 ;

% ----- DERIVED PARAMETERS AND OPERATING POINT -----

% derived parameters

wN=2*pi*fN ;
TB=1/(2*pi*fN) ;

Ldd=Xd ;    % direct axis
Ldf=Xd-Xl ;
Tpdo=Tpdo_s/TB ;
Lff=(Ldf^2)/(Xd-Xpd) ;
Rf=Lff/Tpdo ;
K=Rf/Ldf ;

Lqq=Xq ;    % quadrature axis
Lqq1=Xq-Xl ;
Tpqo=Tpqo_s/TB ;
Lq1q1=(Lqq1^2)/(Xq-Xpq) ;
Rq1=Lq1q1/Tpqo ;

% initial values of differential and algebraic states

phasVinf=V-j*Xe*(P-j*Q)/V ;      % angul ref = machine
Vinf=abs(phasVinf) ;
theta=-angle(phasVinf) ;         % angul ref = inf bus
phasV=V*exp(j*theta) ;           
vx=real(phasV) ;
vy=imag(phasV) ;
phasI=(phasV-Vinf)/(j*Xe) ;
ix=real(phasI) ;
iy=imag(phasI) ;
phasE=phasV+j*Xq*phasI;          % this voltage falls in the q axis
delta=angle(phasE);
id=-ix*sin(delta)+iy*cos(delta) ;
iq= ix*cos(delta)+iy*sin(delta) ;
vd=-vx*sin(delta)+vy*cos(delta) ;
vq= vx*cos(delta)+vy*sin(delta) ;
psid=vq ;
psiq=-vd ;
ifd=(psid-Ldd*id)/Ldf ;
psif=Lff*ifd+Ldf*id ;
vf=Rf*ifd/K ;
iq1=0. ;
psiq1=Lqq1*iq ;

% ----- LINEARIZATION -----

% full (unreduced) Jacobian
% u = setpoint of AVR          z = omega

%   1    2    3    4    5    6    7    8    9    10  11   12   13    14   15   16   17   18
%  vx   vy   ix   iy   vd   vq   id   iq   ifd  iq1 psid psiq psif psiq1 delt omeg  vf    u

J = [ ...
   -1    0    0   -Xe   0    0    0    0    0    0    0    0    0     0    0    0    0    0   ;  %  1 net equ x 
    0   -1   Xe    0    0    0    0    0    0    0    0    0    0     0    0    0    0    0   ;  %  2 net equ y 
    0    0    0    0    0    0   Ldd   0   Ldf   0   -1    0    0     0    0    0    0    0   ;  %  3 psid
    0    0    0    0    0    0    0   Lqq   0   Lqq1  0   -1    0     0    0    0    0    0   ;  %  4 psiq
    0    0    0    0    0    0   Ldf   0   Lff   0    0    0   -1     0    0    0    0    0   ;  %  5 psif
    0    0    0    0    0    0    0  Lqq1   0  Lq1q1  0    0    0    -1    0    0    0    0   ;  %  6 psiq1
    0    0    0    0    0    0    0    0  -Rf*wN 0    0    0    0     0    0    0  K*wN   0   ;  %  7 d psif/dt
    0    0    0    0    0    0    0    0    0 -Rq1*wN 0    0    0     0    0    0    0    0   ;  %  8 d psiq1/dt
    0    0    0    0   -1    0    0    0    0    0    0   -1    0     0    0    0    0    0   ;  %  9 Park d axis
    0    0    0    0    0   -1    0    0    0    0    1    0    0     0    0    0    0    0   ;  % 10 Park q axis
    0    0    0    0    0    0 psiq/(2*H) -psid/(2*H) 0 0 -iq/(2*H) id/(2*H) 0 0 0 0 0    0   ;  % 11 d omega/dt
    0    0    0    0    0    0    0    0    0    0    0    0    0     0    0   wN    0    0   ;  % 12 d delta/dt
 -G*vx/(T*V) -G*vy/(T*V) 0 0 0  0 0    0    0    0    0    0    0     0    0    0  -1/T  G/T  ;  % 13 d vf/dt
 -sin(delta) cos(delta) 0 0 -1  0 0    0    0    0    0    0    0     0   -vq   0    0    0   ;  % 14 vd(vx,vy,delta)
  cos(delta) sin(delta) 0 0  0 -1 0    0    0    0    0    0    0     0    vd   0    0    0   ;  % 15 vq(vx,vy,delta)
 0 0 -sin(delta) cos(delta) 0 0  -1    0    0    0    0    0    0     0   -iq   0    0    0   ;  % 16 id(vx,vy,delta)
 0 0  cos(delta) sin(delta) 0 0   0   -1    0    0    0    0    0     0    id   0    0    0   ;  % 17 iq(vx,vy,delta) 
    0    0    0    0    0    0    0    0    0    0    0    0    0     0    0    1    0    0 ] ;  % 18 z = h(x,y,u)

% linearized model : matrix A

difeq= [12 11  7  8 13] ;
difst= [15 16 13 14 17] ;
algeq= [1 2 3 4 5 6 9 10 14 15 16 17] ;
algst= [1 2 3 4 5 6 7  8  9 10 11 12];
dfdx=J(difeq,difst);
dfdy=J(difeq,algst);
dgdx=J(algeq,difst);
dgdy=J(algeq,algst);
A=dfdx-dfdy*inv(dgdy)*dgdx;

% linearized model : matrices B, C and D

dfdu=J(difeq,18);
dgdu=J(algeq,18);
dhdx=J(18,difst);
dhdy=J(18,algst);
dhdu=J(18,18);

B=dfdu-dfdy*inv(dgdy)*dgdu;
C=dhdx-dhdy*inv(dgdy)*dgdx;
D=dhdu-dhdy*inv(dgdy)*dgdu;
