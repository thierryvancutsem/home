function [A,u] = linmodel(t) ;

global wN Re Le Ldd Ldf Lff Lqq Lqq1 Lq1q1 Loo Ra Rf Rq1 theta0 phasE Vf ;

A=zeros(19,19) ;
theta=theta0+wN*t ;

A(1,1)  = -Re/Le ;   % d ia/dt
A(1,4)  = 1/Le ;
A(2,2)  = -Re/Le ;   % d ib/dt
A(2,5)  = 1/Le ;
A(3,3)  = -Re/Le ;   % d ic/dt
A(3,6)  = 1/Le ;
A(4,4)  = (sqrt(2)/3)*cos(theta) ;         % vd 
A(4,5)  = (sqrt(2)/3)*cos(theta-2*pi/3) ;
A(4,6)  = (sqrt(2)/3)*cos(theta+2*pi/3) ;
A(4,7)  = -1 ;
A(5,4)  = (sqrt(2)/3)*sin(theta) ;         % vq 
A(5,5)  = (sqrt(2)/3)*sin(theta-2*pi/3) ;
A(5,6)  = (sqrt(2)/3)*sin(theta+2*pi/3) ;
A(5,8)  = -1 ; 
A(6,4)  = 1/3 ;                            % vo
A(6,5)  = 1/3 ; 
A(6,6)  = 1/3 ; 
A(6,9)  = -1 ;
A(7,1)  = (sqrt(2)/3)*cos(theta) ;         % id 
A(7,2)  = (sqrt(2)/3)*cos(theta-2*pi/3) ;
A(7,3)  = (sqrt(2)/3)*cos(theta+2*pi/3) ;
A(7,10) = -1 ;
A(8,1)  = (sqrt(2)/3)*sin(theta) ;         % iq 
A(8,2)  = (sqrt(2)/3)*sin(theta-2*pi/3) ;
A(8,3)  = (sqrt(2)/3)*sin(theta+2*pi/3) ;
A(8,11) = -1 ; 
A(9,1)  = 1/3 ;     % io
A(9,2)  = 1/3 ; 
A(9,3)  = 1/3 ;
A(9,12) = -1 ;
A(10,10)= Ldd ;     % psid
A(10,13)= Ldf ;
A(10,15)= -1 ;
A(11,11)= Lqq ;     % psiq
A(11,14)= Lqq1 ;
A(11,16)= -1 ; 
A(12,13)= Lff ;     % psif
A(12,10)= Ldf ;
A(12,17)= -1 ;
A(13,14)= Lq1q1 ;   % psiq1
A(13,11)= Lqq1 ;
A(13,18)= -1 ;
A(14,12)= Loo ;     % psio
A(14,19)= -1 ;
A(15,16)= -1 ;      % d psid/dt
A(15,7) = -1 ;
A(15,10)= -Ra ;
A(16,15)=  1 ;      % d psiq/dt
A(16,8) = -1 ;
A(16,11)= -Ra ;
A(17,13)= -Rf ;     % d psif/dt
A(18,14)= -Rq1 ;    % d psiq1/dt
A(19,12)= -Ra ;     % d psio/dt
A(19,9) = -1 ;

if t<= 0.05  || t>= 99
    ea=sqrt(2)*abs(phasE)*cos(wN*t+angle(phasE)) ;
    eb=sqrt(2)*abs(phasE)*cos(wN*t+angle(phasE)-2*pi/3) ;
    ec=sqrt(2)*abs(phasE)*cos(wN*t+angle(phasE)+2*pi/3) ;
else
    ea=0. ;
    eb=0. ;
    ec=0. ;
end

u= [ -ea/Le -eb/Le -ec/Le 0 0 0 0 0 0 0 0 0 0 0 0 0 Vf 0 0]' ;

end

