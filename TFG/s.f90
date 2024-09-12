program main
    implicit none

    integer i1,i2,L1,L2,n1,n2,j1,j2,C1,C2,MMAX,Cd,cont1,cont2,i,m1,m2,m3,m4,j3,j4,kron
    real*16 VH,IntAng,P,M,Inte,VAng,Norm,delta,Ig,CG,Z


    write(*,*) 'Introduce la dimension de la base: '
	read(*,*) MMAX

    write(*,*) 'Introduce la carga del núcleo '
	read(*,*) Z

    i1=5
    i2=5
    C1=0
    C2=0
    cont1=0
    cont2=0
    M=0

    m1=0
    m2=0
    m3=0
    m4=0

    j1=0
    j2=0
    j3=0
    j4=0



    open(1, file='data.dat', status= 'unknown')
    open(2, file='norm.dat', status= 'unknown')

    cont1=0


    do while (cont1.lt.MMAX)
        do L1=0, (i1-5)/2, 1
        cont1=cont1+1
        print*, cont1
        cont2=0
        i2=5
        C2=0
        
                    do while (cont2.lt.MMAX)
                                    do L2=0, (i2-5)/2, 1
                                    cont2=cont2+1
                                        !SACAMOS VALORES DEL HAMILTONIANO
                                                P=VH(i1,i2,L1,L2,j1,j2,j3,j4,m1,m2,m3,m4,Z)   
                                                write(1,*) P
                                        if(M.gt.P) then
                                        M=P
                                        end if
                                    end do
                        i2=i2+2
                   end do
        end do


        i1=i1+2
    end do


    close(1)
    close(2)

    print *, cont1, M,  (i1-2)*0.5



STOP
end program


! -----------------------------------------------------------------------------------------------
!   Función que calcula el número de estados con energía k  
! -----------------------------------------------------------------------------------------------
function Cd(k)

    implicit none

    integer k,L,Cd,n,l1,kron
    
    Cd=0

    do L=0,(k-5)/2,1
    do n=0,L/2,1
                Cd=Cd+1
    end do
    end do


    
    RETURN
end function








! -----------------------------------------------------------------------------------------------
!   Función que calcula int[cos^a*sen^b]
! -----------------------------------------------------------------------------------------------
function IntAng (a,b)

    implicit none

    integer a,b
    real*16 IntAng

    IntAng=gamma(0.5*(a+1))*gamma(0.5*(b+1))/(2*gamma(0.5*(a+b)+1))

    RETURN
end function



! -----------------------------------------------------------------------------------------------
!   Función que calcula el valor esperado de cos^a*sen^b 
! -----------------------------------------------------------------------------------------------
function VAng (c,b,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)

    implicit none

    integer c,b,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4,i,j,n1,n2
    real*16 VAng,AP,IntAng,S1,S2
    
    S1=0
    S2=0

    if ((v1.eq.v3).and.(v2.eq.v4).and.(m1.eq.m3).and.(m2.eq.m4)) then

    n1=(L1-v1-v2)/2
    n2=(L2-v3-v4)/2

    do i=0, n1,1
        do j=0, n2,1
                if (mod(i+j,2).eq.0) then 
                    S1=S1+AP(n1,v1,v2,i)*AP(n2,v3,v4,j)*IntAng(2+v2+v4+2*(n1-i)+2*(n2-j)+c, 2+v1+v3+2*i+2*j+b)
                else
                    S2=S2+AP(n1,v1,v2,i)*AP(n2,v3,v4,j)*IntAng(2+v2+v4+2*(n1-i)+2*(n2-j)+c, 2+v1+v3+2*i+2*j+b)
                end if
        end do
    end do

    VAng=S1-S2

    else

    VAng=0.e0

    end if

    RETURN
end function






! -----------------------------------------------------------------------------------------------
!   Función que calcula las constantes A(n,l1,l2,s) 
! -----------------------------------------------------------------------------------------------
function AP (n,l1,l2,s) 

    implicit none

    integer n,l1,l2,s,i
    real*16 B,AP,Comb2,z1,z2


    z1=1.e0*n+l1+0.5
    z2=1.e0*n+l2+0.5    

    AP=Comb2(z1,n-s)*Comb2(z2,s)

    B=2*(2.0*n+l1+l2+2)*Gamma(n+l1+l2+2.0)*Gamma(n+1.0)/(Gamma(n+l1+1.5)*Gamma(n+l2+1.5))

    AP=AP*sqrt(B)

    RETURN
end function




! -----------------------------------------------------------------------------------------------
!   Función que calcula las integrales I(mu;L1,L2,k1,k2)
! -----------------------------------------------------------------------------------------------
function Inte (mu,L1,L2,k1,k2,Z)

    implicit none

    integer mu,L1,L2,k1,k2,p1,p2,q1,q2,i,j
    real*16 gamma1,gamma2,c,Inte,S1,S2,alpha1,alpha2,Norm,Comb,fact1,fact2,f1,f2,Z



        gamma1=4.0*Z/k1
        gamma2=4.0*Z/k2
        c=0.5*(gamma1+gamma2)
        S1=0
        S2=0

        p1=(k1-5)/2-L1
        p2=(k2-5)/2-L2

        q1=2*L1+4
        q2=2*L2+4

        fact1=1
        fact2=1

        do i=2,mu,1
            fact1=fact1*i
        end do



        f1=1.0
        alpha1=1


        do i=0,p1,1
        fact2=fact1
        f2=1.0
        alpha2=alpha1




        do j=0,p2,1

            if(mod(i+j,2).eq.0) then

                S1=S1+Comb(p1+q1,p1-i)*Comb(p2+q2,p2-j)*fact2*alpha2/(f1*f2)

             else

                S2=S2+Comb(p1+q1,p1-i)*Comb(p2+q2,p2-j)*fact2*alpha2/(f1*f2)

            end if

        f2=f2*(j+1)
        fact2=fact2*(mu+i+j+1)
        alpha2=alpha2*gamma2/c
        end do





        f1=f1*(i+1)
        fact1=fact1*(mu+i+1)
        alpha1=alpha1*gamma1/c
        end do

        Inte=(S1-S2)/(c**(mu+1))
        Inte=Inte*(gamma1**L1)*(gamma2**L2)




    RETURN  
end function







! -----------------------------------------------------------------------------------------------
!   Función que calcula las combinaciones de a y b
! -----------------------------------------------------------------------------------------------
function Comb(a,b)

    implicit none

    integer a,b,i
    real*16 Comb

    Comb=1.0


    if (b.eq.0) then

        Comb=1.e0
    
    else if(a.eq.0) then

        Comb=1.e0

    else if (a/2.lt.b) then

        do i=0,b-1,1
        
            Comb=Comb*(a-i)/(b-i)
        end do

    else

        do i=a,b+1,-1
            Comb=Comb*(i)/(a-i+1)
        end do

    end if

    RETURN
end function



! -----------------------------------------------------------------------------------------------    
!   Función que calcula los productos cruzados <i|H|j> 
! -----------------------------------------------------------------------------------------------
function VH(k1,k2,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4,Z) 

    implicit none

    integer k1,k2,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4,i,kron,n1,n2,l,r1,r2
    real*16 aux,aux2,gamma1,gamma2,VH,VAng,Norm,Inte,pi,Ig2,CG,V,Z

    pi=3.14159265359

    gamma1=4.0*Z/k1
    gamma2=4.0*Z/k2


    
    VH=0.e0

    !Interacción con el núcleo
    aux=Inte(4+L1+L2,L1,L2,k1,k2,Z)
    aux2=VAng (-1,0,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)+VAng (0,-1,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)
    VH=VH-aux*(2.0*aux2-Z*kron(L1-L2))*kron(v1-v3)*kron(v2-v4)


    VH=VH*Norm(k1,L1,Z)*Norm(k2,L2,Z)

    !Término Ho
    if((k1.eq.k2).and.(L1.eq.L2).and.(v1.eq.v3).and.(v2.eq.v4).and.(m1.eq.m3).and.(m2.eq.m4)) then
        VH=VH-2.0*Z*Z/(k1*k1)
    end if
    

    RETURN
end function



! -----------------------------------------------------------------------------------------------
!   Función que calcula las constantes de normalización N(2k,L) 
! -----------------------------------------------------------------------------------------------
function Norm(k,L,Z)   

    implicit none

    integer k,L,i
    real*16 aux,Norm,gamma,inte,Z

    Norm=1/sqrt(Inte(5+2*L,L,L,k,k,Z))

    RETURN
end function




! -----------------------------------------------------------------------------------------------
!   Función que calcula las combinaciones (z n) con z semientero
! -----------------------------------------------------------------------------------------------
function Comb2(z,n)   

    implicit none

    integer n,i
    real*16 z,Comb2

    
    Comb2=Gamma(z+1)/(Gamma(n+1.0)*Gamma(z-n+1))

    RETURN
end function




! -----------------------------------------------------------------------------------------------
!   Función delta de kronecker
! -----------------------------------------------------------------------------------------------
function kron(a)

    implicit none

    integer a,kron


    if(a.eq.0) then

    kron=1

    else

    kron=0

    end if
    


    RETURN
end function



! -----------------------------------------------------------------------------------------------
!   Función que calcula las integrales angulares entre 0 y pi/4
! -----------------------------------------------------------------------------------------------
recursive function Ig(a,b) result(ang)

    implicit none

    integer a,b,k
    real*16 pi,Comb,ang,aux
    pi=3.14159265359

    ang=0

    if ((a.eq.0).and.(b.eq.0)) then

        ang=pi/4

    else if ((a.eq.0).and.(b.eq.1)) then

        ang=1-1/sqrt(2.0)

    else if ((a.eq.1).and.(b.eq.0)) then

        ang=1/sqrt(2.0)

    else if ((a.eq.1).and.(b.eq.1)) then

        ang=0.25

    else if ((a.gt.1).and.(b.eq.0)) then

        do k=1,a/2,1
            ang=ang+Comb(a/2,k)*Gamma(0.5*k+0.5)/Gamma(0.5*k+1)
        end do

        ang=ang*Gamma(0.5)+pi
        ang=ang/(2**(a/2+2))

    else if ((a.gt.1).and.(b.eq.1)) then

        ang=(sqrt(2.0))**(a+1)
        ang=(1-1/ang)/(b+1)

    else if ((a.eq.0).and.(b.gt.1)) then

        do k=1,b/2,1
            ang=ang+2*(0.5-mod(k,2))*Comb(b/2,k)*Gamma(0.5*k+0.5)/Gamma(0.5*k+1)
        end do

        ang=ang*Gamma(0.5)+pi
        ang=ang/(2**(b/2+2))


    else if ((a.eq.1).and.(b.gt.1)) then

        ang=(sqrt(2.0))**(b+1)
        ang=1/(ang*(b+1))

    else if ((a.gt.1).and.(b.gt.1)) then

        ang=(sqrt(2.0))**(b+a)
        ang=1/(ang*(b+a))

        aux=Ig(a-2,b)
        ang=ang+(a-1)*aux/(a+b)

    else

        write(*,*) 'Falta algo', a, b


    end if
    
    


    RETURN
end function





! -----------------------------------------------------------------------------------------------
!   Función que calcula los coeficientes de Clebsch–Gordan
! -----------------------------------------------------------------------------------------------
function CG(j,m,j1,m1,j2,m2)

    implicit none

    integer j,m,j1,m1,j2,m2,k,km,ko
    real*16 CG,aux,aux1,aux2,fact

    km=1000000

    if ((j1+j2-j).lt.km) then
        km=j1+j2-j
    end if

    if ((j1-m1).lt.km) then
        km=j1-m1
    end if

    if ((j2+m2).lt.km) then
        km=j2+m2
    end if

    ko=0

    if ((j2-j-m1).gt.ko) then
        ko=j2-j-m1
    end if

    if ((j1-j+m2).gt.ko) then
        ko=j1+m2-j
    end if


    CG=0

    if((m1+m2).eq.m) then

        aux1=((2*j+1.e0)*fact(j+j1-j2)*fact(j-j1+j2)/fact(j1+j2+j+1))*fact(j1+j2-j)

        aux2=fact(j+m)*fact(j-m)*fact(j1+m1)*fact(j1-m1)*fact(j2+m2)*fact(j2-m2)

        do k=ko,km,1
            aux=fact(k)*fact(j1+j2-j-k)*fact(j1-m1-k)*fact(j2+m2-k)*fact(j-j2+m1+k)*fact(j-j1-m2+k)
                CG=2*(0.5-mod(k,2))/aux
        end do

        CG=CG*sqrt(aux1)*sqrt(aux2)

    else

        CG=0

    end if

    if((j.lt.abs(j1-j2).or.(j.gt.(j1+j2)))) then
        CG=0
    end if

    


    RETURN
end function


! -----------------------------------------------------------------------------------------------
!   Función factorial
! -----------------------------------------------------------------------------------------------
function fact(a)

    implicit none

    integer a,i
    real*16 fact

    fact=1.e0

    do i=2,a,1
        fact=fact*i
    end do
    

    RETURN
end function






! -----------------------------------------------------------------------------------------------
!   Función que calcula las integrales angulares del r12
! -----------------------------------------------------------------------------------------------
function IG2(L1,L2,j1,j2,j3,j4,j)

    implicit none

    integer L1,L2,j1,j2,j3,j4,n1,n2,k,i,j
    real*16 S1,S2,aux,AP,Ig,Ig2,aux1
    
    S1=0
    S2=0


    n1=(L1-j1-j2)/2
    n2=(L2-j3-j4)/2

    do i=0, n1,1
        do k=0, n2,1
        aux1=Ig(2+j1+j3+2*i+2*k-(j+1),2+j2+j4+2*(n1-i)+2*(n2-k)+j)
        aux=Ig(2+j2+j4+2*(n1-i)+2*(n2-k)-(j+1), 2+j1+j3+2*i+2*k+j)
        aux=aux+aux1
                if (mod(i+k,2).eq.0) then 
                    S1=S1+AP(n1,j1,j2,i)*AP(n2,j3,j4,k)*aux
                else
                    S2=S2+AP(n1,j1,j2,i)*AP(n2,j3,j4,k)*aux
                end if
        end do
    end do

    Ig2=S1-S2



    RETURN
end function


















