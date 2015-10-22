!===============================================================================
!          velocities
!  monomer1 = M1;  monomer2 = M2, CM1 CM2 center of mass
!  CM2 is on the origin
!  b_impact and beta are defined as:
!  the impact parameter in the rest frame of M1: the distance between M1 and projectile of M2
!  the angle between x-axis and projected velocity of M2 on xy plane in the rest frame of M1.


    subroutine calvel(Natom,atomsymbol,velCM1CM2,b_impact,beta,Temperature,tranVect,velo1,velo2,LDebug)
    implicit none
    
    real(kind=8) :: Mass, MostProbVel,atommass(Natom), velCM1CM2, Temperature
    integer :: Natom,atomsymbol(Natom),i
    logical :: LDebug
    real(kind=8) :: b_impact,beta, delta, tranVect(3),len_tranVect
    real(kind=8) :: direction(3), relative_vel(3)
    real(kind=8) :: velo1(Natom,3),velo2(Natom,3)
!delta is the angle between velocity of monomer2 and z-axis in the rest frame of mon1
    
    Mass = 0.0
	do i=1, Natom, 1
	    if (atomsymbol(i)==1) then
	        atommass(i) = 1.0078
	    elseif ( atomsymbol(i) == 6 ) then
	        atommass(i) = 12
	    else
	        write(*,*) 'Only H and C are supported at present'
	        stop
        endif
		Mass = Mass+atommass(i)
	end do


! magnitude of relative velocity
	MostProbVel = sqrt(2*Temperature*3.167D-6*0.263/Mass)
	velCM1CM2 = velCM1CM2*MostProbVel
	if (LDebug .eqv. .True.) then
		write(*,*) '***********'
		write(*,*) 'subroutine calvel'
		write(*,*) 'N, Mass of monomer, Most probable velocity, velCM1CM2'
		write(*,*) Natom, Mass, MostProbVel, velCM1CM2
	end if
	
! components of relative velocity in the rest frame of 	monomer1
    len_tranVect = sqrt(dot_product(tranVect,tranVect))
    delta = asin(b_impact/len_tranVect)
	direction(3)= cos(delta)
	direction(1)= sin(delta)*cos(beta)
	direction(2)= sin(delta)*sin(beta)
	relative_vel = velCM1CM2*direction
	
! velocities of M1 and M2 in the frame of CM, which is stationary relative to laboratory frame
! velocity_M1 - relative_vel/2; velocity_M2 + relative_vel/2
    do i=1, Natom, 1
        velo1(i,:) = velo1(i,:) - relative_vel/2.0
        velo2(i,:) = velo2(i,:) + relative_vel/2.0
    end do
	
	end subroutine calvel
	
!============================================================================
!============================================================================
! velocity debug

    subroutine velocity_debug
    implicit none
    
    write(*,*)	'****************'
    write(*,*)  'subroutine velocity_debug'
    write(*,*)  'in developing'
    
    end subroutine velocity_debug
	
	
	

