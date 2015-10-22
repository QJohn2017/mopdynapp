! dynamics of center of mass of each monomer.
! the coordinates, velocities, accelerations of CM. & bond orders
! the work done by the force between two monomers. !!NO, it is impossible to cal the work!!

    subroutine MassCenterDyn(arg)
    implicit none
    
    integer :: Natom,Nsnap,arg
    logical :: LMasCenDynDebug
    integer, allocatable :: atomsymbol(:)
    character, allocatable :: atomsymbol_c(:)
    real(kind=8),allocatable :: atommass(:),CMcoord(:,:,:),CMvel(:,:,:),CMaccel(:,:,:)
    real(kind=8) :: Mass, Delta
    
    LMasCenDynDebug = .True.
    
    call get_monomer_number(Natom,Nsnap,Delta,arg,LMasCenDynDebug)
    
   	allocate(atomsymbol(Natom),atomsymbol_c(Natom),atommass(Natom))
   	allocate(CMcoord(Nsnap,2,3),CMvel(Nsnap,2,3),CMaccel(Nsnap,2,3))
   	
    call get_atom_symbol(Natom,atomsymbol,atomsymbol_c,atommass,Mass,arg,LMasCenDynDebug) 

    call system('mkdir ./result/')
    
    call monomer_CM_coordinates(Natom,atommass,Mass,Nsnap,CMcoord,arg,LMasCenDynDebug,Delta)
    call monomer_CM_velocities(Natom,atommass,Mass,Nsnap,CMvel,arg,LMasCenDynDebug)
    call monomer_CM_accelerations(Natom,atommass,Mass,Nsnap,CMaccel,arg,LMasCenDynDebug)
    
    call output_CM(Natom,Nsnap,Delta,CMcoord,CMvel,CMaccel,arg,LMasCenDynDebug)
    
	call bondorder(Natom,Nsnap,Delta,atomsymbol,atomsymbol_c,arg,LMasCenDynDebug)
    		   
    end subroutine MassCenterDyn 
    
!================================================================================
!=================================================================================
! get # of atoms of a monomer, and # of snap, Delta = time_{snap i+1} - time_{snap i}
    
    subroutine get_monomer_number(Natom,Nsnap,Delta,arg,LDebug)
    implicit none
    
    integer :: Natom,Nsnap,nlines,arg
    logical :: LDebug
    real(kind=8) :: Delta, t1, t2
    character(256) :: arg_c, filename
    
    write(arg_c,'(I2)') arg
    filename = 'collision_'//trim(adjustl(arg_c))//'.coord'    
    
    open(unit=444,file=trim(adjustl(filename)))
    read(444,*) Natom
    Natom = Natom/2
    close(444)
    
    call system('wc -l '//trim(adjustl(filename))//' > temp.dat')
    open(unit=445,file='temp.dat')
    read(445,*) nlines
    Nsnap = nlines / (Natom*2 + 2)
    close(445)
    call system('rm temp.dat')
    
    call system('grep -i " 1TIME" '//trim(adjustl(filename))//' > temp.dat')
    call system('grep -i " 2TIME" '//trim(adjustl(filename))//' >> temp.dat')
    call system("sed -i 's/FRAME= 1TIME= //g' temp.dat")
    call system("sed -i 's/FRAME= 2TIME= //g' temp.dat")
    call system("sed -i 's/.FSEC//g' temp.dat")
    open(unit=445,file='temp.dat')
    read(445,*) t1
    read(445,*) t2
    Delta = t2 - t1
    close(445)
    call system('rm temp.dat')    
    
   	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine get_monomer_number'
		write(*,*) Natom, nlines, Nsnap, Delta
	endif    		
    
    end subroutine get_monomer_number
!=================================================================================
!=================================================================================
! the type of atoms

    subroutine get_atom_symbol(Natom,atomsymbol,atomsymbol_c,atommass,Mass,arg,LDebug)
    implicit none
    
    integer :: Natom,i, atomsymbol(Natom),arg
    real(kind=8) :: xjunk, yjunk,zjunk,Mass, atommass(Natom)
    character :: atomsymbol_c(Natom)
    character(len=256) :: junk,arg_c,filename
    logical :: LDebug
    
    write(arg_c,'(I2)') arg
    filename = 'collision_'//trim(adjustl(arg_c))//'.coord'    
    
    open(unit=444, file=trim(adjustl(filename)))
    
    read(unit=444,fmt=*) i
    
    read(444,'(A)') junk
    
	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine get_atom_symbol'
		write(*,*) i, trim(adjustl(junk))
	endif    
	
	do i=1, Natom, 1
	    read(444,*) atomsymbol_c(i), xjunk, yjunk, zjunk, atomsymbol(i)
	end do
	
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
	
	close(444)
	
	If ( LDebug .eqv. .True. ) then
		write(*,*) Mass, atomsymbol_c(6),atomsymbol(6)
	endif    		
	
	end subroutine get_atom_symbol

!=============================================================================
!=============================================================================
! the coordinates of center of mass of  monomers in dimerization

    subroutine monomer_CM_coordinates(Natom,atommass,Mass,Nsnap,CMcoord,arg,LDebug,Delta)
    implicit none
    
    integer :: Natom,Nsnap,i,j,arg
    real(kind=8) :: atommass(Natom), Mass, CMcoord(Nsnap,2,3),d, Delta
    logical :: LDebug
    real(kind=8) :: coord1(Natom,3), coord2(Natom,3)
    character(256) :: junk,arg_c,filename
    character :: junkc
    
    write(arg_c,'(I2)') arg
    filename = 'collision_'//trim(adjustl(arg_c))//'.coord'    

    open(unit=444,file=trim(adjustl(filename)))
    
    do i=1, Nsnap, 1
    
        read(444,*) junk
        read(444,*) junk
    
        do j=1, Natom, 1
            read(444,*) junkc, coord1(j,:)
        end do
        call CM_Weighted(Natom,atommass,Mass,coord1,CMcoord(i,1,:))
        
        do j=1, Natom, 1
            read(444,*) junkc, coord2(j,:)
        end do
        call CM_Weighted(Natom,atommass,Mass,coord2,CMcoord(i,2,:)) 
        
    end do
    
    close(444)
    
    call maxdistance(Natom,coord1,coord2,d)
    call outputmaxdistance(d,arg,Nsnap,Delta)
    
    end subroutine monomer_CM_coordinates
!=============================================================================
!=============================================================================
!
    subroutine outputmaxdistance(d,arg,Nsnap,Delta)
    implicit none
    
    real(kind=8) :: d,Delta
    integer :: arg, Nsnap
    character(256) :: arg_c, filename
    
    write(arg_c,'(I2)') arg
    filename = 'collision_'//trim(adjustl(arg_c))//'.dob' ! dimerize or bounce off    
    
    open(unit=451,file='./result/'//trim(adjustl(filename)))
    
    if ( d .gt. 32.0) then
        write(451,*) 'N', d, '! bounce off'
    else 
        write(451,*) 'Y', d, '! NOT bounce off in', Nsnap*Delta, ' fs'
    end if
    
    close(451)
    
    end subroutine outputmaxdistance
!=============================================================================
!=============================================================================
!
    subroutine distance(x,y,d)
    implicit none  
    
    real(kind=8) :: x(3), y(3), d
    
    d = sqrt( (x(1)-y(1))*(x(1)-y(1)) + (x(2)-y(2))*(x(2)-y(2)) + (x(3)-y(3))*(x(3)-y(3)) )
    
    end subroutine distance
!=============================================================================
!=============================================================================
!
    subroutine maxdistance(Natom,coord1,coord2,d)
    implicit none
    
    integer :: Natom,i,j
    real(kind=8) :: coord1(Natom,3), coord2(Natom,3), d, dtemp
    
    d=0.0
    
    do i=1, Natom, 1
        do j=1, Natom, 1
            call distance(coord1(i,:),coord2(j,:),dtemp)
            if ( dtemp .gt. d ) d = dtemp
        end do
    end do
    
    end subroutine maxdistance
!=============================================================================
!=============================================================================
! the velocities of center of mass of monomers in dimerization

    subroutine monomer_CM_velocities(Natom,atommass,Mass,Nsnap,CMvel,arg,LDebug)
    implicit none
    
    integer :: Natom,Nsnap,i,j,arg
    real(kind=8) :: atommass(Natom), Mass, CMvel(Nsnap,2,3)
    logical :: LDebug
    real(kind=8) :: veloc1(Natom,3), veloc2(Natom,3)
    character(256) :: arg_c, filename
!    character(256) :: junk
!    character :: junkc

    write(arg_c,'(I2)') arg
    filename = 'collision_'//trim(adjustl(arg_c))//'.veloc'    
    
    open(unit=446,file=trim(adjustl(filename)))
    
    do i=1, Nsnap, 1
    
        do j=1, Natom, 1
            read(446,*) veloc1(j,:)
        end do
        call CM_Weighted(Natom,atommass,Mass,veloc1,CMvel(i,1,:))
        
        do j=1, Natom, 1
            read(446,*) veloc2(j,:)
        end do
        call CM_Weighted(Natom,atommass,Mass,veloc2,CMvel(i,2,:)) 
        
        read(446,*)
        
    end do
    
    close(446)
    
    end subroutine monomer_CM_velocities

!=============================================================================
!=============================================================================
! the accelerations of center of mass of monomers in dimerization

    subroutine monomer_CM_accelerations(Natom,atommass,Mass,Nsnap,CMaccel,arg,LDebug)
    implicit none
    
    integer :: Natom,Nsnap,i,j,arg
    real(kind=8) :: atommass(Natom), Mass, CMaccel(Nsnap,2,3)
    logical :: LDebug
    real(kind=8) :: accel1(Natom,3), accel2(Natom,3)
    character(256) :: arg_c, filename
!    character(256) :: junk
!    character :: junkc

    write(arg_c,'(I2)') arg
    filename = 'collision_'//trim(adjustl(arg_c))//'.accel'    
    
    open(unit=447,file=trim(adjustl(filename)))
    
    do i=1, Nsnap, 1
    
        do j=1, Natom, 1
            read(447,*) accel1(j,:)
        end do
        call CM_Weighted(Natom,atommass,Mass,accel1,CMaccel(i,1,:))
        
        do j=1, Natom, 1
            read(447,*) accel2(j,:)
        end do
        call CM_Weighted(Natom,atommass,Mass,accel2,CMaccel(i,2,:)) 
        
        read(447,*)
        
    end do
    
    close(447)
    
    end subroutine monomer_CM_accelerations
    
!===========================================================================
!===========================================================================
! CM_Weighted;  \sum_i x_i*m_i / \sum_i m_i, where x_i can be 
! velocity, coordinate, acceleration and so on.
    
    subroutine CM_Weighted(Natom,atommass,Mass,vectorsystem,vector)
    implicit none
    
    integer :: Natom,i
    real(kind=8) :: atommass(Natom), Mass, vectorsystem(Natom,3), vector(3)
    
    vector = 0.0
    
    do i=1, Natom, 1
        vector(:) = vector(:) + vectorsystem(i,:)*atommass(i)
    end do
    
    vector = vector / Mass
    
    end subroutine CM_Weighted
    
!============================================================================
!============================================================================    
! output coordinates, velocities and accelerations of CM1 and CM2    
! i, [time(fs)]  [x,y,z of CM1(\AA)] [x,y,z of CM2(\AA)] [vx,vy,vz,|v| of CM1(\AA/fs)] 
! [vx,vy,vz,|v| of CM2 (\AA/fs)] [ax,ay,az,|a| of CM1(\AA/fs^2)] [ax,ay,az,|a| of CM2(\AA/fs^2)]
! 24 column in total
    
    subroutine output_CM(Natom,Nsnap,Delta,CMcoord,CMvel,CMaccel,arg,LDebug)
    implicit none
    
    integer :: Natom, Nsnap, i,arg
    real(kind=8) :: Delta, CMcoord(Nsnap,2,3), CMvel(Nsnap,2,3), CMaccel(Nsnap,2,3)
    logical :: LDebug
    character(256) :: arg_c,filename
    
    write(arg_c,'(I2)') arg
    filename = 'collision_'//trim(adjustl(arg_c))//'.cm'    
    
    open(unit=448,file='./result/'//trim(adjustl(filename)))
    
    do i=1, Nsnap, 1
        write(448,*) i, i*Delta, CMcoord(i,1,:), CMcoord(i,2,:), CMvel(i,1,:), &
     &   sqrt(dot_product(CMvel(i,1,:),CMvel(i,1,:))), CMvel(i,2,:),  &
     &   sqrt(dot_product(CMvel(i,2,:),CMvel(i,2,:))), CMaccel(i,1,:),  &
     &   sqrt(dot_product(CMaccel(i,1,:),CMaccel(i,1,:))), CMaccel(i,2,:), &
     &   sqrt(dot_product(CMaccel(i,2,:),CMaccel(i,2,:)))
    end do
    
    close(448)
    
    end subroutine output_CM

    
    
