Sb= 10e6 ;           % single-phase base power (VA)

Vb(1:3)= 36e3/sqrt(3) ; % vector of phase-neutral base voltages (V)
Vb(7:9)= Vb(1) ;
Vb(4:6)= 150e3/sqrt(3) ;
Vb(10:12)= 6e3/sqrt(3) ;
Vb(13:15)= 15e3/sqrt(3) ;

Zb=Vb.^2/Sb ;    % vector of base impedances (ohm)

YcableBC= Yline(0.909/Zb(1), 1.659/Zb(1), 645.1e-06*Zb(1), 7.870/Zb(1), ...
    3.470/Zb(1), 645.1e-06*Zb(1)) ;

k= 3*Sb/27e6 ;
YtrfoBA= Ytrfo(0.005*k, 0.11*k, 0, 0.95, pi/6, 0, 1, 0.005*k, 0.175*k, ...
    1e6, 1e6, 0, 0) ;

k= 3*Sb/20e6 ;
YtrfoCD= Ytrfo(0.006*k, 0.10*k, 0, 1.03, 0, 0, 2, 0.006*k, 0.15*k, ...
    1e6, 1e6, 1e6, 1e6) ;

k= 3*Sb/10e6 ;
YtrfoCE= Ytrfo(0.006*k, 0.126*k, 0, 0.97, pi/6, 0, 1, 0.006*k, ...
    0.136*k, 1e6, 1e6, 0, 0) ;

p= 15e6/(3*Sb) ; q= 7e6/(3*Sb) ;
YloadD= Yload(0, 1e6, 1e6, p, q, 1.0011372) ;

pcc= 3000e6/(3*Sb) ; p= 9.21e6/(3*Sb) ; q= 5.75e6/(3*Sb) ;
[Ynet,Inet]= Norton_equ(0., 1/pcc, p, q, 1, 0) ;

k= 3*Sb/10e6 ; p= 6e6/(3*Sb) ; q= 2e6/(3*Sb) ;
[Ygen,Igen]= Norton_gen(0.005*k, 0.13*k, 0.010*k, 0.13*k, 0.005*k, ...
    0.07*k, 1e6, 1e6, p, q, 1.0093312, 0.02972656) ;

Y=zeros(15,15) ;
Y(1:6,1:6)= Y(1:6,1:6)+YtrfoBA ;
ind= [1 2 3 7 8 9] ;
Y(ind,ind)= Y(ind,ind)+YcableBC ;
ind= [7 8 9 10 11 12] ;
Y(ind,ind)= Y(ind,ind)+YtrfoCD ;
ind= [7 8 9 13 14 15] ;
Y(ind,ind)= Y(ind,ind)+YtrfoCE ;
ind= [4 5 6] ;
Y(ind,ind)= Y(ind,ind)+Ynet ;
ind= [10 11 12] ;
Y(ind,ind)= Y(ind,ind)+YloadD ;
ind= [13 14 15] ;
Y(ind,ind)= Y(ind,ind)+Ygen ;

I= zeros(15,1) ;
I(4:6)= Inet ;
I(13:15)= Igen ;

V=Y\I;


