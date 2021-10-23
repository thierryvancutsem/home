function [A] = matA() ;

global Re Le Ldd Ldf Lff Lqq Lqq1 Lq1q1 Ra Rf Rq1 theta0 ;

A=zeros(14,14) ;

A(1,3)=-Rf ;          % d psif/dt
A(2,4)=-Rq1;          % d psiq1/dt
A(3,7)=1 ;            % net equ. x comp
A(3,11)=-Re ;
A(3,12)=Le ;
A(4,8)=1 ;            % net equ. y comp
A(4,12)=-Re ;
A(4,11)=-Le ;
A(5,5)=1 ;            % vd(vx,vy)
A(5,7)=-cos(theta0) ;
A(5,8)=-sin(theta0) ;
A(6,6)=1 ;            % vq(vx,vy)
A(6,7)=-sin(theta0) ;
A(6,8)=cos(theta0) ;
A(7,9)=1 ;            % id(vx,vy)
A(7,11)=-cos(theta0) ;
A(7,12)=-sin(theta0) ;
A(8,10)=1 ;           % iq(vx,vy)
A(8,11)=-sin(theta0) ;
A(8,12)=cos(theta0) ;
A(9,5)=1 ;            % Park - d axis
A(9,9)=Ra ;
A(9,14)=1 ;
A(10,6)=1 ;           % Park - q axis
A(10,10)=Ra ;
A(10,13)=-1 ;
A(11,13)=1 ;          % psid(id,if)
A(11,9)=-Ldd ; 
A(11,3)=-Ldf ;
A(12,14)=1 ;          % psiq(iq)
A(12,10)=-Lqq ;
A(12,4)=-Lqq1 ;
A(13,1)=1 ;           % psif(id,idf)
A(13,3)=-Lff ;
A(13,9)=-Ldf ;
A(14,2)=1 ;           % psiq1(iq,iq1)
A(14,4)=-Lq1q1 ;
A(14,10)=-Lqq1 ;
