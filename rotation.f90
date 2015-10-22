!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
! rotation and translation of monomers.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! rotate a point about origin point by Eular angle EAin
    subroutine rotation_point(x,y,z,EAin)
    implicit none
    

    type :: typeEulerAngle
        real(kind = 8) :: phi, theta, psi !notation in Goldstein's book
    end type typeEulerAngle
    type(typeEulerAngle) :: EAin, EA
    real(kind=8) :: x,y,z, xnew, ynew, znew
    real(kind=8) :: R(3,3)
    
    EA%phi = EAin%phi / 180 * 3.1415926
    EA%theta = EAin%theta / 180 * 3.1415926
    EA%psi = EAin%psi / 180 * 3.1415926
    
!    R(1,1) = cos(EA%psi)*cos(EA%phi)-cos(EA%theta)*sin(EA%phi)*sin(EA%psi)
!    R(1,2) = cos(EA%psi)*sin(EA%phi)+cos(EA%theta)*cos(EA%phi)*sin(EA%psi)
!    R(1,3) = sin(EA%psi)*sin(EA%theta)
!    R(2,1) = -sin(EA%psi)*cos(EA%phi)-cos(EA%theta)*sin(EA%phi)*cos(EA%psi)
!    R(2,2) = -sin(EA%psi)*sin(EA%phi)+cos(EA%theta)*cos(EA%phi)*cos(EA%psi)
!    R(2,3) = cos(EA%psi)*sin(EA%theta)
!    R(3,1) = sin(EA%theta)*sin(EA%phi)
!    R(3,2) = -sin(EA%theta)*cos(EA%phi)
!    R(3,3) = cos(EA%theta)

    R(1,1) = cos(EA%psi)*cos(EA%phi)-cos(EA%theta)*sin(EA%phi)*sin(EA%psi)
    R(1,2) = -sin(EA%psi)*cos(EA%phi)-cos(EA%theta)*sin(EA%phi)*cos(EA%psi)
    R(1,3) = sin(EA%theta)*sin(EA%phi)
    R(2,1) = cos(EA%psi)*sin(EA%phi)+cos(EA%theta)*cos(EA%phi)*sin(EA%psi)
    R(2,2) = -sin(EA%psi)*sin(EA%phi)+cos(EA%theta)*cos(EA%phi)*cos(EA%psi)
    R(2,3) = -sin(EA%theta)*cos(EA%phi)
    R(3,1) = sin(EA%theta)*sin(EA%psi)
    R(3,2) = sin(EA%theta)*cos(EA%psi)
    R(3,3) = cos(EA%theta)

    xnew = R(1,1)*x + R(1,2)*y + R(1,3)*z
    ynew = R(2,1)*x + R(2,2)*y + R(2,3)*z
    znew = R(3,1)*x + R(3,2)*y + R(3,3)*z
    
    x = xnew
    y = ynew
    z = znew
    
    end subroutine rotation_point

!==============================================================================
!==============================================================================
! rotation a rigid-body about origin by Euler angle EA

    subroutine rotation_rigid_body(N,coord,EAin)
    implicit none
    
    type :: typeEulerAngle
        real(kind = 8) :: phi, theta, psi !notation in Goldstein's book
    end type typeEulerAngle
    type(typeEulerAngle) :: EAin
    integer :: N, i 
    real(kind=8) :: coord(N,3)
!    logical :: LDebug
    
    do i=1, N, 1
        call rotation_point(coord(i,1),coord(i,2),coord(i,3),EAin)
    end do    
    
    end subroutine rotation_rigid_body
 
!==============================================================================
!==============================================================================
! rotation matrix of a given Euler angle

    subroutine RotMatrix(EA,R)    
    implicit none
    
    type :: typeEulerAngle
        real(kind = 8) :: phi, theta, psi !notation in Goldstein's book
    end type typeEulerAngle
    type(typeEulerAngle) :: EA    
    real(kind=8) :: R(3,3)
    
    R(1,1) = cos(EA%psi)*cos(EA%phi)-cos(EA%theta)*sin(EA%phi)*sin(EA%psi)
    R(1,2) = -sin(EA%psi)*cos(EA%phi)-cos(EA%theta)*sin(EA%phi)*cos(EA%psi)
    R(1,3) = sin(EA%theta)*sin(EA%phi)
    R(2,1) = cos(EA%psi)*sin(EA%phi)+cos(EA%theta)*cos(EA%phi)*sin(EA%psi)
    R(2,2) = -sin(EA%psi)*sin(EA%phi)+cos(EA%theta)*cos(EA%phi)*cos(EA%psi)
    R(2,3) = -sin(EA%theta)*cos(EA%phi)
    R(3,1) = sin(EA%theta)*sin(EA%psi)
    R(3,2) = sin(EA%theta)*cos(EA%psi)
    R(3,3) = cos(EA%theta)
    
    end subroutine RotMatrix
    
!==============================================================================
!==============================================================================
! product of two Euler angle. EA3 = EA1 * EA2. Pay attention to the order of operation!
! rotates EA2 firstly, and then rotates EA1

    subroutine EulAngProduct(EA3,EA1,EA2)
    implicit none
    
    type :: typeEulerAngle
        real(kind = 8) :: phi, theta, psi !notation in Goldstein's book
    end type typeEulerAngle
    type(typeEulerAngle) :: EA3,EA1,EA2
    real(kind=8) :: R3(3,3), R1(3,3), R2(3,3)
    
    call RotMatrix(EA1,R1)
    call RotMatrix(EA2,R2)
    
    R3 = MATMUL(R1,R2)
    
    end subroutine EulAngProduct
 !=============================================================================
 !=============================================================================
 ! for debug purpose. generate a movie of continuous rotation.
 
    subroutine rotation_debug(N,coord,atomsymbol)
    implicit none
    
    integer :: N, i, j
    real(kind=8) :: coord(N,3), coord_temp(N,3)
    integer :: atomsymbol(N)
    character :: atomsymbol_c(N)
    type :: typeEulerAngle
        real(kind = 8) :: phi, theta, psi !notation in Goldstein's book
    end type typeEulerAngle
    type(typeEulerAngle) :: EAin
    
    write(*,*) '*******************'
    write(*,*) 'subroutine rotation_debug'
    write(*,*) 'Use JMOL to visualize phi.xyz, theta.xyz and psi.xyz'
    
    do i=1, N, 1
        if( atomsymbol(i) == 1 ) then
            atomsymbol_c(i) = 'H'
        elseif (atomsymbol(i)==6) then
            atomsymbol_c(i) = 'C'
        else
            stop
            write(*,*) 'only C and H are supported at present'
        end if
    end do
    
    open(unit=449,file='./Debug/phi.xyz')
    do i=1, 11, 1
        EAin%phi=(i-1)*18.0
        EAin%theta=0.0
        EAin%psi=0.0
        coord_temp = coord
        call rotation_rigid_body(N,coord_temp,EAin)
        write(449,'(I2)') N
        write(449,'("phi=",F5.1,"theta=",F3.1,"psi=",F3.1)') EAin%phi, EAin%theta, EAin%psi
        do j=1,N,1
            write(449,*) atomsymbol_c(j), coord_temp(j,:)
        end do
    end do
    close(449)
    
    open(unit=450, file='./Debug/theta.xyz')
    do i=1, 11, 1
        EAin%phi=45.0
        EAin%theta=(i-1)*9.0        
        EAin%psi=0.0
        coord_temp = coord
        call rotation_rigid_body(N,coord_temp,EAin)
        write(450,'(I2)') N
        write(450,'("phi=",F4.1,"theta=",F5.1,"psi=",F3.1)') EAin%phi, EAin%theta, EAin%psi
        do j=1,N,1
            write(450,*) atomsymbol_c(j), coord_temp(j,:)
        end do
    end do
    close(450)    
    
    open(unit=451, file='./Debug/psi.xyz')
    do i=1, 11, 1
        EAin%phi=45.0
        EAin%theta=45.0        
        EAin%psi=(i-1)*18.0
        coord_temp = coord
        call rotation_rigid_body(N,coord_temp,EAin)
        write(451,'(I2)') N
        write(451,'("phi=",F4.1,"theta=",F4.1,"psi=",F5.1)') EAin%phi, EAin%theta, EAin%psi
        do j=1,N,1
            write(451,*) atomsymbol_c(j), coord_temp(j,:)
        end do
    end do
    close(451)    
    
    end subroutine rotation_debug
    
!=============================================================================
!=============================================================================
! translate a point

    subroutine translation_point(x,y,z,tranVect)
    implicit none
    
    real(kind=8) :: x, y, z, tranVect(3), xnew, ynew, znew
    
    xnew = x+tranVect(1)
    ynew = y+tranVect(2)
    znew = z+tranVect(3)
    
    x=xnew
    y=ynew
    z=znew
    
    end subroutine translation_point
   
!=============================================================================
!=============================================================================
! translate a rigid-body

    subroutine translation_rigid_body(N,coord,tranVect)
    implicit none
    
    integer :: N,i
    real(kind=8) :: coord(N,3), tranVect(3)
    
    do i=1, N, 1
       call translation_point(coord(i,1),coord(i,2),coord(i,3),tranVect)
    end do
    
    end subroutine translation_rigid_body
    
!=============================================================================
!=============================================================================
! Debug: rotation + translation.

    subroutine coordinate_Debug(N,coord1,EulAng1,tranVect,coord2,EulAng2,atomsymbol)    
    implicit none
    
    integer :: N, i, j, atomsymbol(N)
    real(kind=8) :: coord1(N,3), coord2(N,3), coord_temp(N,3), tranVect(3), tV(3),coordin(N,3)
    character :: atomsymbol_c(N)
    type :: typeEulerAngle
        real(kind = 8) :: phi, theta, psi !notation in Goldstein's book
    end type typeEulerAngle
    type(typeEulerAngle) :: EulAng1, EulAng2, EAin
    
    write(*,*) "*********************"
    write(*,*) "subroutine coordinate_Debug"
    write(*,*) "Use JMOL to visualize coordDebug.xyz"
    
    do i=1, N, 1
        if( atomsymbol(i) == 1 ) then
            atomsymbol_c(i) = 'H'
        elseif (atomsymbol(i)==6) then
            atomsymbol_c(i) = 'C'
        else
            stop
            write(*,*) 'only C and H are supported at present'
        end if
    end do    
    
    open(unit=452,file='./Debug/coordDebug.xyz')
    do i=0, 10, 1
        EAin%phi = EulAng1%phi*i/10.0
        EAin%theta = 0.0
        EAin%psi = 0.0
        coord_temp = coord1
        call rotation_rigid_body(N,coord_temp,EAin)
        write(452,'(I2)') N
        write(452,'("phi=",F5.1,"theta=",F3.1,"psi=",F3.1)') EAin%phi, EAin%theta, EAin%psi
        do j=1,N,1
            write(452,*) atomsymbol_c(j), coord_temp(j,:)
        end do
    end do
    
    do i=0, 10, 1
        EAin%phi = EulAng1%phi
        EAin%theta = EulAng1%theta*i/10.0
        EAin%psi = 0.0
        coord_temp = coord1
        call rotation_rigid_body(N,coord_temp,EAin)
        write(452,'(I2)') N
        write(452,'("phi=",F5.1,"theta=",F5.1,"psi=",F3.1)') EAin%phi, EAin%theta, EAin%psi
        do j=1,N,1
            write(452,*) atomsymbol_c(j), coord_temp(j,:)
        end do
    end do
    
    do i=0, 10, 1
        EAin%phi = EulAng1%phi
        EAin%theta = EulAng1%theta
        EAin%psi = EulAng1%psi*i/10.0
        coord_temp = coord1
        call rotation_rigid_body(N,coord_temp,EAin)
        write(452,'(I2)') N
        write(452,'("phi=",F5.1,"theta=",F5.1,"psi=",F5.1)') EAin%phi, EAin%theta, EAin%psi
        do j=1,N,1
            write(452,*) atomsymbol_c(j), coord_temp(j,:)
        end do
    end do    
    
    coordin =coord1
    call rotation_rigid_body(N,coordin,EulAng1)
    
    do i=0, 10, 1
        tV = tranVect*i/10.0
        coord_temp = coordin
        call translation_rigid_body(N,coord_temp,tV)
        write(452,'(I2)') N
        write(452,'(F4.1,F4.1,F5.1)') tV
        do j=1,N,1
            write(452,*) atomsymbol_c(j), coord_temp(j,:)
        end do
    end do
    
    call translation_rigid_body(N,coordin, tranVect)
    
    
    do i=0, 10, 1
        EAin%phi = EulAng2%phi*i/10.0
        EAin%theta = 0.0
        EAin%psi = 0.0
        coord_temp = coord2
        call rotation_rigid_body(N,coord_temp,EAin)
        write(452,'(I2)') N*2
        write(452,'("phi=",F5.1,"theta=",F3.1,"psi=",F3.1)') EAin%phi, EAin%theta, EAin%psi
        do j=1,N,1
            write(452,*) atomsymbol_c(j), coordin(j,:)
        end do
        do j=1,N,1
            write(452,*) atomsymbol_c(j), coord_temp(j,:)
        end do
    end do
    
    do i=0, 10, 1
        EAin%phi = EulAng2%phi
        EAin%theta = EulAng2%theta*i/10.0
        EAin%psi = 0.0
        coord_temp = coord2
        call rotation_rigid_body(N,coord_temp,EAin)
        write(452,'(I2)') N*2
        write(452,'("phi=",F5.1,"theta=",F5.1,"psi=",F3.1)') EAin%phi, EAin%theta, EAin%psi
        do j=1,N,1
            write(452,*) atomsymbol_c(j), coordin(j,:)
        end do
        do j=1,N,1
            write(452,*) atomsymbol_c(j), coord_temp(j,:)
        end do
    end do
    
    do i=0, 10, 1
        EAin%phi = EulAng2%phi
        EAin%theta = EulAng2%theta
        EAin%psi = EulAng2%psi*i/10.0
        coord_temp = coord2
        call rotation_rigid_body(N,coord_temp,EAin)
        write(452,'(I2)') N*2
        write(452,'("phi=",F5.1,"theta=",F5.1,"psi=",F5.1)') EAin%phi, EAin%theta, EAin%psi
        do j=1,N,1
            write(452,*) atomsymbol_c(j), coordin(j,:)
        end do
        do j=1,N,1
            write(452,*) atomsymbol_c(j), coord_temp(j,:)
        end do
    end do    
   
    close(452)
   
   end subroutine coordinate_Debug    
    
!============================================================================
!============================================================================
! visualization.

    subroutine visualize(N,coord1,coord2,atomsymbol,arg)
    implicit none
    integer :: N, i, atomsymbol(N),arg
    real(kind=8) :: coord1(N,3), coord2(N,3)
    character :: atomsymbol_c(N)
    character(256) :: arg_c, filename
    
    write(arg_c,'(I2)') arg
    filename = 'collision_initial_'//trim(adjustl(arg_c))//'.xyz'
    
    write(*,*) "*********************"
    write(*,*) "subroutine visualize"
    write(*,*) "Use JMOL to visualize collision.xyz"
    
    do i=1, N, 1
        if( atomsymbol(i) == 1 ) then
            atomsymbol_c(i) = 'H'
        elseif (atomsymbol(i)==6) then
            atomsymbol_c(i) = 'C'
        else
            stop
            write(*,*) 'only C and H are supported at present'
        end if
    end do    
    
    open(unit=453,file='./Debug/'//trim(adjustl(filename)))    
    write(453,'(I2)') N*2
    write(453,*) 'initial positions of collision'
    do i=1,N,1
        write(453,*) atomsymbol_c(i), coord1(i,:)
    end do
    do i=1,N,1
        write(453,*) atomsymbol_c(i), coord2(i,:)
    end do   
    close(453)
    
    end subroutine visualize
    
    

