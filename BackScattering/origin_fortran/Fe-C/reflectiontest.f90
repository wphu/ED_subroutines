! This code is suitable for the material with the mass number smaller than 92
program reflectiontest
implicit none

	
    real :: rn(91,500),reenergy(91,500)
    integer :: nz1, m1,nz2(1), nw(1), ne
    integer :: i, j
    real*8  energy, theta,ee,re
	real*8, external :: rnion,reion


	!	nz2	= array for atomic numbers of the constituents.
	!	nw	= array for relative numbers of the constituents.
    !   theta 	= angle of incidence ions with the normal to the plasma facing surface,  in degree.
    open(10,file="reflection_rate.txt")
    open(20,file="energy_rate.txt")
	
	nz1 = 26
    m1 = 56
	
    ne =1
    nz2(1) = 6
    nw(1) = 1
	
	
	
    do i = 1, 91
        do j =1, 500
            theta = 1.0*(i-1)
            energy = 0.1*j
            Rn(i,j)= rnion( theta,energy,nz1,m1,ne,nz2,nw ) !  number backscattering coefficient of light ions.
            reenergy(i,j)=reion( theta,energy,nz1,m1,ne,nz2,nw )
            if(Rn(i,j)>= 1) Rn(i,j) = 1
        end do
        write(10,"(500d)") rn(i,:)
        write(20,"(500d)") reenergy(i,:)
    end do
	
	
	
    theta=65.53
    energy=0.1 !0.8565
    ee=reion( theta,energy,nz1,m1,ne,nz2,nw )
    re= rnion( theta,energy,nz1,m1,ne,nz2,nw )
    write(*,*) energy,ee,re
!    theta = 0.
!    do i = 1, 500
!            energy = 0.1*dble(i)
!            Rn(i,1)= rnion( theta,energy,nz1,m1,ne,nz2,nw ) !  number backscattering coefficient of light ions.
!            if(Rn(i,1)>= 1) Rn(i,1) = 1
!            write(10,"(2d)") energy, rn(i,1)
!    end do



    energy = 100.
    theta = 0.
    nz1 = 1
    m1 = 2
    ne = 1
    nz2(1)= 74
    nw(1) = 1


	Rn(1,1)= rnion( theta,energy,nz1,m1,ne,nz2,nw ) !  number backscattering coefficient of light ions.
	write(*,*) "backscattering coefficient" , rn(1,1)
	stop


print *, "end of reflection test"
end program reflectiontest
