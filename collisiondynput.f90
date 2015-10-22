! read the initial conditions of two monomers in collision from    
! an input file named 'ColDynGen.ini', and generate the *.dynput 
! 131012 distance of CM and relative velocity of CM, only

	subroutine collision_dynput_generate(arg)
	implicit none

	logical :: LColDynGenDebug
	character(len=256) :: ColDynGen_flag
	real(kind=8) :: ColDynGen_temp, Temperature ! , disCM1CM2, xcoorCM2, ycoorCM2
!	real(kind=8) :: theta, phi, alpha
    type :: typeEulerAngle
        real(kind = 8) :: phi, theta, psi !notation in Goldstein's book
    end type typeEulerAngle
    type(typeEulerAngle) :: EulAng1, EulAng2
	real(kind=8) :: velCM1CM2, b_impact, beta
!	real(kind=8) :: Mass, MostProbVel
	integer :: Natom, Nsnap, i
	real(kind=8),allocatable :: coord1(:,:), coord2(:,:), velo1(:,:), velo2(:,:)
	integer,allocatable :: atomsymbol(:)
!	real(kind=8),allocatable :: atommass(:)
	real(kind=8) :: tranVect(3)
	integer :: arg
	
	LColDynGenDebug = .False.
	tranVect(1) = 0
	tranVect(2) = 0
	tranVect(3) = 12
	
!	read ColDynGen.ini
!	call read_ColDynGen_ini_old(ColDynGen_flag, ColDynGen_temp,disCM1CM2, &
!	     xcoorCM2, ycoorCM2,theta, phi, alpha,velCM1CM2, b_impact, beta, &
!	     LColDynGenDebug)

	call read_ColDynGen_ini(ColDynGen_flag, ColDynGen_temp,EulAng1, &
	  &   EulAng2 ,velCM1CM2, b_impact, beta, tranVect(3),arg,LColDynGenDebug)
!	Extract information form relax.dynput, 
!	and write first 15 lines in collision.dynput
	call Extract_inform_dynput(Natom, Nsnap, Temperature,arg,LColDynGenDebug)

!	Extract the coordinates and velocities of 
!	the Nsnap-th frame and the (Nsnap-10)-th frame
	allocate(coord1(Natom,3),coord2(Natom,3))	
	allocate(velo1(Natom,3),velo2(Natom,3))
	allocate(atomsymbol(Natom))
!	allocate(atommass(Natom))
	call Extract_coord_vel(Natom,Nsnap,coord1,coord2,velo1,velo2,atomsymbol, &
		 LColDynGenDebug)

!!========OLD VERSION==============================================
!!      MAIN PART
!!	coord1 and velo1 are unchanged.
!!	modify coord2 and velo2 to arbitrary CM-position, direction and CM-velocity
!!=======================================================
!!=======================================================
!!==========NEW VERSION==================================
!!          MAIN PART
!! coord1 and coord2 rotate EulAng1 and EulAng2, respectively
!! coord1 translate by tranVect, which is set to be (0,0,12)
!!=======================================================
!	most probable speed \sqrt{ 2 k_b T / M 
	
    call system('mkdir ./Debug/')
!   debug the codes on rigid body rotation	
    if( LColDynGenDebug .eqv. .True. ) then
        call rotation_debug(Natom, coord1, atomsymbol)
        call coordinate_Debug(Natom,coord1,EulAng1,tranVect,coord2,EulAng2,atomsymbol)
    end if
    
    call rotation_rigid_body(Natom,coord1,EulAng1)
    call translation_rigid_body(Natom,coord1,tranVect)
    call rotation_rigid_body(Natom,coord2,EulAng2) !CM of monomer2 at origin
    
!    if( LColDynGenDebug .eqv. .True. ) then
    if( .True. ) then
        call visualize(Natom,coord1,coord2,atomsymbol,arg)
    end if    

!!=====================================================================
!!==   velocities. 
    call calvel(Natom,atomsymbol,velCM1CM2,b_impact,beta,Temperature,tranVect,velo1,velo2,LColDynGenDebug)

    if( LColDynGenDebug .eqv. .True. ) then
        call velocity_debug
    end if

    call writedynput(Natom,coord1,coord2,velo1,velo2,atomsymbol,arg)

    write(*,"('collision_',I2,'.dynput has been generated')") arg
    
	end subroutine
	
	
!====================================================================	
!====================================================================
!====================================================================
!read ColDynGen.ini file
	subroutine read_ColDynGen_ini_old(ColDynGen_flag, ColDynGen_temp,disCM1CM2, &
	     xcoorCM2, ycoorCM2,theta, phi, alpha,velCM1CM2, b_impact, beta, LDebug)
	implicit none

	logical :: LDebug
	character(len=256) :: ColDynGen_flag
	real(kind=8) :: ColDynGen_temp, disCM1CM2, xcoorCM2, ycoorCM2
	real(kind=8) :: theta, phi, alpha
	real(kind=8) :: velCM1CM2, b_impact, beta

	open(unit=444,file='ColDynGen.ini')
	read(444,*) ColDynGen_flag
	read(444,*) ColDynGen_temp
	read(444,*) disCM1CM2
	read(444,*) xcoorCM2
	read(444,*) ycoorCM2
	read(444,*) theta
	read(444,*) phi
	read(444,*) alpha
	read(444,*) velCM1CM2
	read(444,*) b_impact
	read(444,*) beta
	close(444)

	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine read_ColDynGen_ini'
		write(*,*) trim(adjustl(ColDynGen_flag)), ColDynGen_temp
		write(*,*) disCM1CM2, xcoorCM2, ycoorCM2
		write(*,*) theta, phi, alpha
		write(*,*) velCM1CM2, b_impact, beta
	endif

	end subroutine read_ColDynGen_ini_old


!====================================================================
!read ColDynGen.ini file
	subroutine read_ColDynGen_ini(ColDynGen_flag, ColDynGen_temp,EulAng1, &
	     EulAng2 ,velCM1CM2, b_impact, beta, tranVectz,arg,LDebug)
	implicit none

    type :: typeEulerAngle
        real(kind = 8) :: phi, theta, psi !notation in Goldstein's book
    end type typeEulerAngle
    type(typeEulerAngle) :: EulAng1, EulAng2
	logical :: LDebug
	character(len=256) :: ColDynGen_flag,filename,arg_c
	real(kind=8) :: ColDynGen_temp !, disCM1CM2, xcoorCM2, ycoorCM2
!	real(kind=8) :: theta, phi, alpha
	real(kind=8) :: velCM1CM2, b_impact, beta, tranVectz
	integer :: arg

    write(arg_c,'(I2)') arg
    filename = 'ColDynGen_'//trim(adjustl(arg_c))//'.ini'
    
	open(unit=444,file=trim(adjustl(filename)))
	read(444,*) ColDynGen_flag
	read(444,*) ColDynGen_temp
	read(444,*) EulAng1%phi
	read(444,*) EulAng1%theta
	read(444,*) EulAng1%psi
	read(444,*) EulAng2%phi
	read(444,*) EulAng2%theta
	read(444,*) EulAng2%psi
	read(444,*) velCM1CM2
	read(444,*) b_impact
	read(444,*) beta
	read(444,*) tranVectz
	close(444)

	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine read_ColDynGen_ini'
		write(*,*) trim(adjustl(filename))
		write(*,*) trim(adjustl(ColDynGen_flag)), ColDynGen_temp
		write(*,*) EulAng1
		write(*,*) EulAng2
		write(*,*) velCM1CM2, b_impact, beta
		write(*,*) tranVectz
	endif

	end subroutine read_ColDynGen_ini


!=============================================================
!=============================================================
! obtain # of atoms and # of snap in relax.dynput
! and write first 15 lines in collision.dynput
! the fisrt line in relax.dynput should be '*_relax*', there must be '_relax' in the first line

	subroutine Extract_inform_dynput(Natom,Nsnap,Temperature,arg,LDebug)
	implicit none

	integer :: Natom, Nsnap, n_, nl, i,arg
	real(kind=8) :: junk1,junk2,Temperature
	character(len=256) :: symbol,arg_c,filename
	logical :: LDebug

    write(arg_c,'(I2)') arg
    filename = 'collision_'//trim(adjustl(arg_c))//'.dynput'

	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine Extract_inform_dynput'
	endif

	open(unit=445,file='relax.dynput')
	open(unit=446,file=trim(adjustl(filename)))
	
! the 1st line
	read(445,'(A)') symbol
	n_ = scan(trim(adjustl(symbol)),'_')
	nl = len(trim(adjustl(symbol)))
	do i= nl,n_+6,-1
		symbol(i+4:i+4) = symbol(i:i)
	end do
	symbol(n_+1:n_+9) = 'collision'
	write(446,'(A)') trim(adjustl(symbol))//'_'//trim(adjustl(arg_c))

	If ( LDebug .eqv. .True. ) then
		write(*,*) 'lenght of 1st line', nl, 'position of prefix', n_
		write(*,*) trim(adjustl(symbol))
	endif
	
! the 2nd-8th line
	do i=2,8,1
		read(445,'(A)') symbol
		write(446,'(A)') trim(adjustl(symbol))
	end do

! the 9th line
	read(445,*) Natom, n_, nl
	write(446,'(T3,I4,T8,I4,T13,I4,T41,"! N,NIMAGE,NRECT")') Natom*2, n_, nl

! the 10th line
	read(445,*) nl, junk1, junk2, Temperature, n_
	write(446,'(T5,I1,T9,F8.6,T17,F8.6,T27,F6.1,T33,I6,T41,"! ICF,VMAX,ZNOISE,TREF,IDENRC")') &
	      nl,junk1,junk2,Temperature,n_
	If ( LDebug .eqv. .True. ) then
		write(*,*) Temperature
	endif

! the 11th line
	read(445,*) junk1, nl, Nsnap
	write(446,'(T4,F9.6,T11,I5,T20,I6,T41,"! DT,NTIME,NSNAP")') junk1, nl, Nsnap

! the 12th-15th line
	do i=12,15,1
		read(445,'(A)') symbol
		write(446,'(A)') trim(adjustl(symbol))
	end do

	close(445)
	close(446)
	


	end subroutine Extract_inform_dynput

!=============================================================
!=============================================================
!	Extract the coordinates and velocities of 
!	the Nsnap-th frame and the (Nsnap-10)-th frame
!	from relax.coord and relax.veloc

	subroutine Extract_coord_vel(Natom,Nsnap,coord1,coord2,velo1,velo2,  &
			   atomsymbol, LDebug)
	implicit none

	logical :: LDebug
	integer :: Natom, Nsnap, i
	real(kind=8) :: coord1(Natom,3), coord2(Natom,3), velo1(Natom,3), velo2(Natom,3)
	character(256) :: junk
	integer :: atomsymbol(Natom)
	character :: junk1

!	coordinates

	open(unit=447,file='relax.coord')

	do i=1,(Natom+2)*(Nsnap-10),1
		read(447,*)
	end do

	read(447,*) junk
	read(447,'(A)') junk

	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine Extract_coord_vel'
		write(*,*) trim(adjustl(junk))
	endif

	do i=1,Natom,1
		read(447,*) junk1, coord1(i,:), atomsymbol(i)
	end do

	if ( LDebug .eqv. .True.) then
		write(*,*) coord1(2,:), atomsymbol(2)
		write(*,*) 'coordinates of monomer 1: read'
	end if

	do i=1, (Natom+2)*8, 1
		read(447,*)
	end do

	read(447,*) junk
	read(447,'(A)') junk

	If ( LDebug .eqv. .True. ) then
		write(*,*) trim(adjustl(junk))
	endif

	do i=1,Natom,1
		read(447,*) junk1,coord2(i,:)
	end do

	if ( LDebug .eqv. .True.) then
		write(*,*) coord2(2,:), atomsymbol(2)
		write(*,*) 'coordinates of monomer 2: read'
	end if

	close(447)

!	velocities
	
	open(unit=448,file='relax.veloc')

	do i=1,(Natom+1)*(Nsnap-10),1
		read(448,*)
	end do

	do i=1,Natom,1
		read(448,*) velo1(i,:)
	end do

	if ( LDebug .eqv. .True.) then
		write(*,*) velo1(2,:), atomsymbol(2)
		write(*,*) 'velocities of monomer 1: read'
	end if

	read(448,*)

	do i=1, (Natom+1)*8, 1
		read(448,*)
	end do

	do i=1,Natom,1
		read(448,*) velo2(i,:)
	end do

	if ( LDebug .eqv. .True.) then
		write(*,*) velo2(2,:), atomsymbol(2)
		write(*,*) 'velocities of monomer 2: read'
	end if
	
	close(448)

	end subroutine Extract_coord_vel
!=============================================================
!=============================================================
!  write the atom specification in collision.dynput

    subroutine writedynput(Natom,coord1,coord2,velo1,velo2,atomsymbol,arg)
    implicit none
    
    integer :: Natom, atomsymbol(Natom), i, icol, jfix, arg
    real(kind=8) :: coord1(Natom,3), coord2(Natom,3), velo1(Natom,3), velo2(Natom,3)
    character(256) :: arg_c, filename
    
    write(arg_c,'(I2)') arg
    filename = 'collision_'//trim(adjustl(arg_c))//'.dynput'    
    
    icol = 3
    jfix = 0
    
    open(unit=446, file=trim(adjustl(filename)), position='APPEND')
    
    do i=1, Natom, 1
        write(446,'(T4,F10.6,T17,F10.6,T30,F10.6,T43,I3,T49,I3,T55,I3)') coord1(i,:), atomsymbol(i), icol, jfix
    end do
    
    do i=1, Natom, 1
        write(446,'(T4,F10.6,T17,F10.6,T30,F10.6,T43,I3,T49,I3,T55,I3)') coord2(i,:), atomsymbol(i), icol, jfix
    end do
    
    do i=1, Natom, 1
        write(446,'(T4,F10.6,T17,F10.6,T30,F10.6)') velo1(i,:) 
    end do
    
    do i=1, Natom, 1
        write(446,'(T4,F10.6,T17,F10.6,T30,F10.6)') velo2(i,:) 
    end do
    
    close(446)
    
    end subroutine writedynput
    
    
    

