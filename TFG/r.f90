program main
    implicit none

    integer i1,i2,L1,L2,n1,n2,j1,j2,C1,C2,MMAX,Cd,cont1,cont2,i,m1,m2,m3,m4,j3,j4,kron,lmax
    real*8 VH,IntAng,P,M,Inte,VAng,Norm,delta,Ig,CG,Ig2


    write(*,*) 'Introduce la dimension de la base: '
	read(*,*) MMAX

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

    lmax=0



    open(1, file='data.dat', status= 'unknown')
    open(2, file='norm.dat', status= 'unknown')


    i1=5
    C1=0

    do while (cont1.lt.MMAX)
        do L1=0, (i1-5)/2, 1
        do j1=0,lmax,1
        do j2=0,lmax,1
        do m1=-j1,j1,1
        do m2=-j2,j2,1
        n1=(L1-j1-j2)/2 
        if ((n1.ge.0).and.(mod((L1-j1-j2),2).eq.0)) then
        print *, cont1
        cont1=cont1+1
        cont2=0
        i2=5
        C2=0
        
                    do while (cont2.lt.MMAX)
                                    do L2=0, (i2-5)/2, 1
                                    do j3=0,lmax,1
                                    do j4=0,lmax,1
                                    do m3=-j3,j3,1
                                    do m4=-j4,j4,1
                                    n2=(L2-j3-j4)/2
                                    if ((n2.ge.0).and.(mod((L2-j3-j4),2).eq.0)) then
                                        cont2=cont2+1
                                        !SACAMOS VALORES DEL HAMILTONIANO
                                                P=VH(i1,i2,L1,L2,j1,j2,j3,j4,m1,m2,m3,m4)   
                                                write(1,*) P
                                        if(M.gt.P) then
                                        M=P
                                        end if
                                    end if
                                        
                                    end do
                                    end do
                                    end do
                                    end do
                                    end do
                        i2=i2+2
                   end do
        end if
        end do
        end do
        end do
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
!   Función que calcula int[cos^a*sen^b]
! -----------------------------------------------------------------------------------------------
function IntAng (a,b)

    implicit none

    integer a,b
    real*8 IntAng

    IntAng=gamma(0.5*(a+1))*gamma(0.5*(b+1))/(2*gamma(0.5*(a+b)+1))

    RETURN
end function



! -----------------------------------------------------------------------------------------------
!   Función que calcula el valor esperado de cos^a*sen^b 
! -----------------------------------------------------------------------------------------------
function VAng (c,b,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)

    implicit none

    integer c,b,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4,i,j,n1,n2
    real*8 VAng,AP,IntAng,S1,S2
    
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
    real*8 B,AP,Comb2,z1,z2


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
function Inte (mu,L1,L2,k1,k2)

    implicit none

    integer mu,L1,L2,k1,k2,p1,p2,q1,q2,i,j
    real*8 gamma1,gamma2,c,Inte,S1,S2,alpha1,alpha2,Norm,Comb,fact



        gamma1=8.e0/k1
        gamma2=8.e0/k2
        c=0.5*(gamma1+gamma2)
        S1=0
        S2=0

        p1=(k1-5)/2-L1
        p2=(k2-5)/2-L2

        q1=2*L1+4
        q2=2*L2+4


        alpha1=1


        do i=0,p1,1
        alpha2=alpha1

        do j=0,p2,1

            if(mod(i+j,2).eq.0) then

                S1=S1+Comb(p1+q1,p1-i)*Comb(p2+q2,p2-j)*fact(i+j+mu)*alpha2/(fact(i)*fact(j))

             else

                S2=S2+Comb(p1+q1,p1-i)*Comb(p2+q2,p2-j)*fact(i+j+mu)*alpha2/(fact(i)*fact(j))

            end if
        alpha2=alpha2*gamma2/c
        end do

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
    real*8 Comb

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
function VH(k1,k2,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4) 

    implicit none

    integer k1,k2,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4,i,kron,n1,n2,l,r1,r2,m
    real*8 aux,aux2,gamma1,gamma2,VH,VAng,Norm,Inte,pi,Ig2,CG,V,IJ

    pi=3.14159265359

    gamma1=8.0/k1
    gamma2=8.0/k2


    
    VH=0.e0

    !Interacción con el núcleo
    aux=Inte(4+L1+L2,L1,L2,k1,k2)
    aux2=VAng (-1,0,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)+VAng (0,-1,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)
    VH=VH-2.e0*aux*(aux2-1.e0*kron(L1-L2))*kron(v1-v3)*kron(v2-v4)*kron(m1-m3)*kron(m2-m4)

    !Repulsión entre los electrones
    V=0.e0
    
    r1=abs(v1-v3)
    r2=v1+v3

    if (abs(v2-v4).gt.r1) then
        r1=abs(v2-v4)
    end if

    if ((v2+v4).lt.r2) then
        r2=(v2+v4)
    end if

    do l=r1,r2,1
        aux=1.e0*(2*l+1)*(2*l+1)*(2*v1+1)*(2*v2+1)*(2*v3+1)*(2*v4+1)
        aux=sqrt(aux)*IJ(l,m1+m3,v1,m1,v3,m3)*IJ(l,0,v1,0,v3,0)*2*(0.5-mod(m1+m3,2))
        aux=aux*IJ(l,-m1-m3,v2,m2,v4,m4)*IJ(l,0,v2,0,v4,0)
        V=V+Ig2(L1,L2,v1,v2,v3,v4,l)*aux/(2*l+1)
    end do

    VH=VH+V*Inte(4+L1+L2,L1,L2,k1,k2)
    VH=VH*Norm(k1,L1)*Norm(k2,L2)

    !Término Ho
    if((k1.eq.k2).and.(L1.eq.L2).and.(v1.eq.v3).and.(v2.eq.v4).and.(m1.eq.m3).and.(m2.eq.m4)) then
        VH=VH-8.0/(k1*k1)
    end if

    

    RETURN
end function


function integrate_cos_sin(n, m, a, b) result(integral)
  implicit none
  integer, intent(in) :: n, m
  real*8, intent(in) :: a, b
  real*8 :: integral

  ! Método de Simpson para la integración numérica
  real*8 :: h, x, sum
  integer :: i, n_intervals

  n_intervals = 1000  ! Número de intervalos para la integración

  h = (b - a) / n_intervals
  sum = function_to_integrate(a, n, m) + function_to_integrate(b, n, m)

  do i = 1, n_intervals-1, 2
    x = a + i * h
    sum = sum + 4.0d0 * function_to_integrate(x, n, m)
  end do

  do i = 2, n_intervals-2, 2
    x = a + i * h
    sum = sum + 2.0d0 * function_to_integrate(x, n, m)
  end do

  integral = h / 3.0d0 * sum

contains

  ! Función a integrar: (cos(x))^n * (sin(x))^m
  function function_to_integrate(x, n, m) result(y)
    real*8, intent(in) :: x
    integer, intent(in) :: n, m
    real*8 :: y

    y = (cos(x))**n * (sin(x))**m
  end function function_to_integrate

end function integrate_cos_sin





! -----------------------------------------------------------------------------------------------
!   Función que calcula las constantes de normalización N(2k,L) 
! -----------------------------------------------------------------------------------------------
function Norm(k,L)   

    implicit none

    integer k,L,i
    real*8 aux,Norm,gamma,inte

    Norm=1/sqrt(Inte(5+2*L,L,L,k,k))

    RETURN
end function



! -----------------------------------------------------------------------------------------------
!   Función que calcula los 3j
! -----------------------------------------------------------------------------------------------
function IJ(j,m,j1,m1,j2,m2)   

    implicit none

    integer j,m,j1,j2,m2,m1
    real*8 IJ,CG

    IJ=CG(j,-m,j1,m1,j2,m2)*2*(0.5-mod(j1-j2-m,2))/sqrt(1.d0*(2*j+1))

    RETURN
end function




! -----------------------------------------------------------------------------------------------
!   Función que calcula las combinaciones (z n) con z semientero
! -----------------------------------------------------------------------------------------------
function Comb2(z,n)   

    implicit none

    integer n,i
    real*8 z,Comb2

    
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
function Ig(a,b)

    implicit none

    integer a,b,k
    real*8 pi,Comb,ang,aux,c,d,integrate_cos_sin,Ig
    pi=3.14159265359

    c=0.d0
    d=pi/4.d0

    Ig=integrate_cos_sin(a, b, c, d)
    


    RETURN
end function





! -----------------------------------------------------------------------------------------------
!   Función que calcula los coeficientes de Clebsch–Gordan
! -----------------------------------------------------------------------------------------------
function CG(j,m,j1,m1,j2,m2)

    implicit none

    integer j,m,j1,m1,j2,m2,k,km,ko
    real*8 CG,aux,aux1,aux2,fact

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
    real*8 fact

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
    real*8 S1,S2,aux,AP,Ig,Ig2,aux1
    
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


















