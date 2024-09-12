program main
    implicit none

    integer i1,i2,L1,L2,n1,n2,j1,j2,C1,C2,MMAX,Cd,cont1,cont2,i,m1,m2,m3,m4,j3,j4,kron
    real*8 VH,IntAng,P,M,Inte,VAng,Norm,delta


    write(*,*) 'Introduce la dimension de la base: '
	read(*,*) MMAX

    i1=5
    i2=5
    C1=0
    C2=0
    cont1=0
    cont2=0
    M=0



    open(1, file='data.dat', status= 'unknown')
    open(2, file='norm.dat', status= 'unknown')

    do while (C1.lt.MMAX)
        C1=C1+Cd(i1)
        do L1=0, (i1-5)/2, 1
        do n1=0,L1/2,1
        do j1=0,L1-2*n1,1
        j2=L1-2*n1-j1
        do m1=-j1,j1,1
        do m2=-j2,j2,1
        if ((mod(n1,2).eq.1).and.(j1.eq.j2).and.(m1.eq.m2)) then
        else
        cont1=cont1+1
        cont2=0
        i2=5
        C2=0
        
                    do while (C2.lt.MMAX)
                        C2=C2+Cd(i2)
                                    do L2=0, (i2-5)/2, 1
                                    do n2=0,L2/2,1
                                    do j3=0,L2-2*n2,1
                                    j4=L2-2*n2-j3
                                    do m3=-j3,j3,1
                                    do m4=-j4,j4,1
                                    if ((mod(n2,2).eq.1).and.(j3.eq.j4).and.(m3.eq.m4)) then
                                    else
                                    cont2=cont2+1
                                        !SACAMOS VALORES DEL HAMILTONIANO
                                        if (cont1.le.cont2) then
                                                P=VH(i1,i2,L1,L2,j1,j2,j3,j4,m1,m2,m3,m4)
                                                
                                                write(1,*) P
                                        !SACAMOS PRODUCTOS ESCALARES
                                        else    
                            
                                            P=0.5*kron(i1-i2)*kron(L1-L2)*delta(L1,L2,j1,j2,j3,j4,m1,m2,m3,m4)
                                            write(2,*) p
                                        end if
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
!   Función que calcula el número de estados con energía k  
! -----------------------------------------------------------------------------------------------
function Cd(k)

    implicit none

    integer k,L,Cd,n,l1,kron
    
    Cd=0

    do L=0,(k-5)/2,1
    do n=0,L/2,1
        if (mod(n,2).eq.1) then 
            do l1=0,(L-2*n),1
                Cd=Cd+(2*l1+1)*(2*(L-2*n-l1)+1-kron(L-2*n-2*l1))
            end do
        else
            do l1=0,(L-2*n),1
                Cd=Cd+(2*l1+1)*(2*(L-2*n-l1)+1)
            end do
        end if
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
    real*8 gamma1,gamma2,c,Inte,S1,S2,alpha1,alpha2,Norm,Comb,fact1,fact2,f1,f2



        gamma1=8.0/k1
        gamma2=8.0/k2
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
!   Función que calcula las integrales de dr
! -----------------------------------------------------------------------------------------------
function Partr (mu,L1,L2,k1,k2)

    implicit none

    integer mu,L1,L2,k1,k2,p1,p2,q1,q2,i,j
    real*8 gamma1,gamma2,c,Partr,S1,S2,alpha1,alpha2,Norm,Comb,fact1,fact2,f1,f2



        gamma1=8.0/k1
        gamma2=8.0/k2
        c=0.5*(gamma1+gamma2)
        S1=0
        S2=0

        p1=(k1-5)/2-L1
        p2=(k2-5)/2-L2

        q1=2*L1+4
        q2=2*L2+4



        f1=1.0
        alpha1=1


        do i=0,p1,1
        f2=1.0
        alpha2=alpha1




        do j=0,p2,1

            if(mod(i+j,2).eq.0) then

                S1=S1+Comb(p1+q1,p1-i)*Comb(p2+q2,p2-j)*alpha2*( (L1+i)*Gamma(1.e0*mu+i+j)*c-0.5*Gamma(1.e0*mu+i+j+1) )/(f1*f2)

             else

                S2=S2+Comb(p1+q1,p1-i)*Comb(p2+q2,p2-j)*alpha2*((L1+i)*Gamma(1.e0*mu+i+j)*c-0.5*Gamma(1.e0*mu+i+j+1))/(f1*f2)

            end if

        f2=f2*(j+1)
        alpha2=alpha2*gamma2/c
        end do





        f1=f1*(i+1)
        fact1=fact1*(mu+i+1)
        alpha1=alpha1*gamma1/c
        end do

        Partr=(S1-S2)/(c**(mu+1))
        Partr=Partr*(gamma1**(L1+1))*(gamma2**L2)





    RETURN  
end function






! -----------------------------------------------------------------------------------------------
!   Función que calcula las integrales de dr²
! -----------------------------------------------------------------------------------------------
function Partr2 (mu,L1,L2,k1,k2)

    implicit none

    integer mu,L1,L2,k1,k2,p1,p2,q1,q2,i,j
    real*8 gamma1,gamma2,c,S1,S2,alpha1,alpha2,Norm,Comb,fact1,fact2,f1,f2,Partr2,aux



        gamma1=8.0/k1
        gamma2=8.0/k2
        c=0.5*(gamma1+gamma2)
        S1=0
        S2=0

        p1=(k1-5)/2-L1
        p2=(k2-5)/2-L2

        q1=2*L1+4
        q2=2*L2+4



        f1=1.0
        alpha1=1


        do i=0,p1,1
        f2=1.0
        alpha2=alpha1




        do j=0,p2,1

            if(mod(i+j,2).eq.0) then

                aux=(L1+i)*(L1+i-1)*Gamma(1.e0*mu+i+j-1)*c*c+0.25*Gamma(1.e0*mu+i+j+1)-(L1+i)*Gamma(1.e0*mu+i+j)*c
                S1=S1+Comb(p1+q1,p1-i)*Comb(p2+q2,p2-j)*alpha2*aux/(f1*f2)

             else

                aux=(L1+i)*(L1+i-1)*Gamma(1.e0*mu+i+j-1)*c*c+0.25*Gamma(1.e0*mu+i+j+1)-(L1+i)*Gamma(1.e0*mu+i+j)*c
                S2=S2+Comb(p1+q1,p1-i)*Comb(p2+q2,p2-j)*alpha2*aux/(f1*f2)

            end if

        f2=f2*(j+1)
        alpha2=alpha2*gamma2/c
        end do





        f1=f1*(i+1)
        fact1=fact1*(mu+i+1)
        alpha1=alpha1*gamma1/c
        end do

        Partr2=(S1-S2)/(c**(mu+1))
        Partr2=Partr2*(gamma1**(L1+2))*(gamma2**L2)





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

    integer k1,k2,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4,i,kron,n1,n2
    real*8 aux,aux2,gamma1,gamma2,VH,VAng,Norm,Inte,Partr,Partr2,delta

    gamma1=8.0/k1
    gamma2=8.0/k2


    
   VH=0.e0

    
    
    !(1)x(1)
    if ((v1.eq.v2).and.(v3.eq.v4).and.(m1.eq.m2).and.(m3.eq.m4)) then

        if ((k1.eq.k2).and.(L1.eq.L2)) then

            !Término cinético
            VH=-8.0/(k1*k1)

            !Interacción culombiana
            aux=Inte(L1+L2+4,L1,L2,k1,k2)
            VH=VH-2*aux*(VAng (-1,0,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)+VAng (0,-1,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)-0.5)
    
            !Normalización
            VH=VH*Norm(k1,L1)*Norm(k2,L2)


            else if ((v1.eq.v3).and.(v2.eq.v4).and.(m1.eq.m3).and.(m2.eq.m4)) then

            !Interacción culombiana
            aux=Inte(L1+L2+4,L1,L2,k1,k2)
            VH=VH-2*aux*(VAng (-1,0,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)+VAng (0,-1,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)-0.5*kron(L1-L2))

            !Normalización
            VH=VH*Norm(k1,L1)*Norm(k2,L2)
        

        else 

               !Interacción culombiana
               aux=Inte(L1+L2+4,L1,L2,k1,k2)
               VH=-2*aux*(VAng (-1,0,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)+VAng (0,-1,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)-0.5*kron(L1-L2))

               !Normalización
               VH=VH*Norm(k1,L1)*Norm(k2,L2)


        end if


    !(1)x(2)=0

    
    !(2)x(2)
    !Solo se consideran estos casos porque si l1!=l2, por ejemplo, entonces l1!=l3 o l2!=l4
    else if (((v1.ne.v2).and.(v3.ne.v4)).or.((m1.ne.m2).and.(m3.ne.m4))) then

        if ((k1.eq.k2).and.(L1.eq.L2)) then

            VH=-4.0/(k1*k1)*delta(L1,L2,v1,v2,v3,v4,m1,m2,m3,m4) 

        end if

        n1=(L1-v1-v2)/2
        n2=(L2-v3-v4)/2

        
        aux=Inte(L1+L2+4,L1,L2,k1,k2)
        
        aux2= VAng (-1,0,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)+VAng (0,-1,L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)
        aux2=aux2+((VAng (-1,0,L1,L2,v2,v1,v3,v4,m2,m1,m3,m4)+VAng (0,-1,L1,L2,v2,v1,v3,v4,m2,m1,m3,m4))*2*(0.5-mod(n1,2)))
        aux2=aux2+((VAng (-1,0,L1,L2,v1,v2,v4,v3,m1,m2,m4,m3)+VAng (0,-1,L1,L2,v1,v2,v4,v3,m1,m2,m4,m3))*2*(0.5-mod(n2,2)))
        aux2=aux2+((VAng (-1,0,L1,L2,v2,v1,v4,v3,m2,m1,m4,m3)+VAng (0,-1,L1,L2,v2,v1,v4,v3,m2,m1,m4,m3))*2*(0.5-mod(n1+n2,2)))

        VH=VH+aux*(kron(L1-L2)*delta(L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)-aux2)

        VH=VH*Norm(k1,L1)*Norm(k2,L2)


    end if

    RETURN
end function



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
!   Función auxiliar de VH
! -----------------------------------------------------------------------------------------------
function delta(L1,L2,v1,v2,v3,v4,m1,m2,m3,m4)

    implicit none

    integer L1,L2,v1,v2,v3,v4,m1,m2,m3,m4,n1,n2
    real*8 delta

    n1=(L1-v1-v2)/2
    n2=(L2-v3-v4)/2

    if ((v1.eq.v3).and.(v2.eq.v4).and.(m1.eq.m3).and.(m2.eq.m4)) then

        delta=1+2*(0.5-mod(n1+n2,2))

    else if ((v1.eq.v4).and.(v2.eq.v3).and.(m1.eq.m4).and.(m2.eq.m3)) then

        delta=2*(0.5-mod(n1,2))+2*(0.5-mod(n2,2))

    else


    delta=0

    end if
    


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




















