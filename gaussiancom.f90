! gaussiancom log2com *.log : *.log of a monomer -> *.com of a dimer
! gaussiancom xyz2com  *.xyz Nthframe: the Nth frame in *.xyz -> *.com
! gaussiancom log2dynput *.log :  *.log -> *.dynput, useful to generate dynput for thermalization
!
! There should be ColDynGen_0.ini in the working directory



    program gaussiancom
    implicit none
    
	integer :: argue_count, argue_int, Nthframe
	character(len=256) :: argue_type, argue_filename, Nthframe_c
	logical :: LGauComDebug
	character(len=256) :: ColDynGen_flag
	real(kind=8) :: ColDynGen_temp
    type :: typeEulerAngle
        real(kind = 8) :: phi, theta, psi !notation in Goldstein's book
    end type typeEulerAngle
    type(typeEulerAngle) :: EulAng1, EulAng2
	real(kind=8) :: velCM1CM2, b_impact, beta
	integer :: Natom, Nsnap, i
	real(kind=8),allocatable :: coord1(:,:), coord2(:,:), velo1(:,:), velo2(:,:)
	integer,allocatable :: atomsymbol(:)
	character,allocatable :: atomsymbol_c(:)
!	real(kind=8),allocatable :: atommass(:)
	real(kind=8) :: tranVect(3)
	
	LGauComDebug = .True.
	tranVect(1) = 0
	tranVect(2) = 0
	tranVect(3) = 12

! read the command
	argue_count = command_argument_count()
	if ( (argue_count .ne. 2) .and. (argue_count .ne. 3) ) then
		write(*,*) 'Number of command argument should be two or three'
		stop
	end if
    write(*,*) argue_count
    	
	call get_command_argument(1,argue_type)
	call get_command_argument(2,argue_filename)
	
!	write(*,*) trim(adjustl(argue_filename))

    
    if( argue_type == 'log2com' ) then ! *.log of a monomer -> *.com of a dimer

        call read_ColDynGen_ini(ColDynGen_flag, ColDynGen_temp,EulAng1, &
    	  &   EulAng2 ,velCM1CM2, b_impact, beta, tranVect(3),0,LGauComDebug)
    
        call get_monomer_number_log(Natom,trim(adjustl(argue_filename)),  &
                len(trim(adjustl(argue_filename))),LGauComDebug)
    
        allocate(atomsymbol_c(Natom),atomsymbol(Natom),coord1(Natom,3), coord2(Natom,3))
    
        call readmonomerlog(Natom,coord1,atomsymbol,trim(adjustl(argue_filename)),  &
               len(trim(adjustl(argue_filename))),LGauComDebug)
               
        call get_atomsymbol_c(Natom,atomsymbol,atomsymbol_c,LGauComDebug)
        
        coord2 = coord1
        
        call rotation_rigid_body(Natom,coord1,EulAng1)
        call translation_rigid_body(Natom,coord1,tranVect)
        call rotation_rigid_body(Natom,coord2,EulAng2) !CM of monomer2 at origin
        
        call writedimercom(Natom,atomsymbol_c,coord1,coord2,trim(adjustl(argue_filename)),  &
                len(trim(adjustl(argue_filename))), LGauComDebug)


    elseif (argue_type == 'xyz2com') then ! the Nth frame in *.xyz -> *.com
    
        call get_command_argument(3,Nthframe_c)
        
        read(Nthframe_c,*) Nthframe
        
!        write(*,*) Nthframe_c, Nthframe
    
        call get_monomer_number_xyz(Natom,trim(adjustl(argue_filename)),  &
                len(trim(adjustl(argue_filename))),LGauComDebug)
                
        allocate(atomsymbol_c(Natom),coord1(Natom,3))
        
        call readxyzNthframe(Nthframe,Natom,atomsymbol_c,coord1,trim(adjustl(argue_filename)),  &
                len(trim(adjustl(argue_filename))),LGauComDebug)
                
        call writexyz2com(Nthframe,Natom,atomsymbol_c,coord1,trim(adjustl(argue_filename)),  &
                len(trim(adjustl(argue_filename))), LGauComDebug)


    elseif (argue_type == 'log2dynput') then ! *.log -> *.dynput, useful to generate dynput for thermalization

        call get_monomer_number_log(Natom,trim(adjustl(argue_filename)),  &
                len(trim(adjustl(argue_filename))),LGauComDebug)
    
        allocate(atomsymbol(Natom),coord1(Natom,3))
    
        call readmonomerlog(Natom,coord1,atomsymbol,trim(adjustl(argue_filename)),  &
               len(trim(adjustl(argue_filename))),LGauComDebug)
               
        call writelog2dynput(Natom,atomsymbol,coord1,trim(adjustl(argue_filename)),  &
                len(trim(adjustl(argue_filename))), LGauComDebug)
               
                       
    else
        stop 'only "log2com", "xyz2com" and "log2dynput" are allowed'
    end if
    
    
    end program
    
!==============================================================================
!==============================================================================
! get the number of atoms in a monomer from *.log file

    subroutine get_monomer_number_log(Natom,filename,length,LDebug)
    implicit none    
    
    integer :: Natom, length, junk
    character(length) :: filename
    logical :: LDebug
    character(256) :: aline
    real(kind=8) :: junkr
    
    Natom = 0
    
	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine get_monomer_number_log'
	endif
    
    open(unit=444,file=trim(adjustl(filename)))
    
    do while (.True.)
        read(444,'(A)') aline
        if ( trim(adjustl(aline)) == 'Standard orientation:') then
            write(*,*) trim(adjustl(aline))
            read(444,*)
            read(444,*)
            read(444,*)
            read(444,*)
!            read(444,'(A)') aline
            do while ( aline(2:2) .ne. '-' )
                read(444,'(A)') aline
              	If ( LDebug .eqv. .True. ) then
!                    write(*,*) trim(adjustl(aline))
                    write(*,*) trim(aline)
               	endif
                if (aline(2:2) .ne. '-') then
                    read(aline,'(i7,i11,18x,3f12.6)') Natom,junk,junkr,junkr,junkr
!                    write(*,*) Natom,junk,junkr                    
                endif
            end do        
            exit
         end if
    end do
    
    close(444)
    
	If ( LDebug .eqv. .True. ) then
		write(*,*) Natom
	endif
    
    end subroutine get_monomer_number_log
!===============================================================================
!===============================================================================
!   read information from *.log, which is the output file of Gaussian09

    subroutine readmonomerlog(Natom,coord1,atomsymbol,filename,length,LDebug)
    implicit none

    integer :: Natom, i,length,junk, Nstandard
    integer ::  atomsymbol(Natom)
    real(kind=8) :: coord1(Natom,3), junkr
!    character :: atomsymbol_c(Natom)
    character(length) :: filename
    logical :: LDebug
    character(256) :: aline

	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine readmonomerlog'
	endif    
    
    Nstandard = 0
    
    open(unit=444,file=trim(adjustl(filename)))
    
    do while (.True.)
        read(444,'(A)',end=50) aline
        if ( trim(adjustl(aline)) == 'Standard orientation:') then
            Nstandard = Nstandard +1
            i = 0
!            write(*,*) trim(adjustl(aline))
            read(444,*)
            read(444,*)
            read(444,*)
            read(444,*)
!            read(444,'(A)') aline
            do while ( i < Natom  )
                i = i +1
                read(444,'(A)') aline
!              	If ( LDebug .eqv. .True. ) then
!                    write(*,*) trim(adjustl(aline))
!                    write(*,*) trim(aline)
!               endif
                read(aline,'(i7,i11,18x,3f12.6)') junk,atomsymbol(i),coord1(i,1),coord1(i,2),coord1(i,3)
            end do        
         end if
    end do
    
50  close(444)

	If ( LDebug .eqv. .True. ) then
		do i = 1, Natom,1
    		write(*,'(i7,3f12.6)') atomsymbol(i), coord1(i,1),coord1(i,2),coord1(i,3)
    	end do
    	write(*,*) Nstandard
	endif

    end subroutine readmonomerlog

!=====================================================================================
!=====================================================================================
! atomsymbol -> atomsymbol_c

    subroutine get_atomsymbol_c(Natom,atomsymbol,atomsymbol_c,LDebug)
    implicit none
    
    integer :: Natom, atomsymbol(Natom),i
    character :: atomsymbol_c(Natom)
    logical :: LDebug
    

    do i=1, Natom, 1
        if( atomsymbol(i) == 1 ) then
            atomsymbol_c(i) = 'H'
        elseif (atomsymbol(i)==6) then
            atomsymbol_c(i) = 'C'
        else
            stop
            write(*,*) 'only C and H are supported at present'
        end if
    end do

    end subroutine get_atomsymbol_c
    
!==================================================================================
!===================================================================================
! write *.com of dimer

    subroutine writedimercom(Natom,atomsymbol_c,coord1,coord2,filein,length,LDebug)
    implicit none
    
    integer :: Natom, i, pos, length
    character :: atomsymbol_c(Natom)
    real(kind=8) :: coord1(Natom,3), coord2(Natom,3)
    logical :: LDebug
    character(length) :: filein
    character(length+1) :: fileout
    
    
    pos=index(filein,'.')
    fileout(1:pos-1) = filein(1:pos-1)
    fileout(pos:pos) = 'd'
    fileout(pos+1:length+1) = '.com'
    
	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine writedimercom'
	endif    
    
    open(unit=445,file=fileout)
    write(445,*) '%nprocshared=12'
    write(445,*) '%mem=32GB'
    write(445,*) '# opt upm6'
    write(445,*) ''
    write(445,*) 'r51 dimer with Euler angles and CM distance given in ColDynGen.ini'
    write(445,*) ''
    write(445,*) '0 1'
    
    do i=1, Natom, 1
!        write(445,'(2A,3F11.9)') atomsymbol_c(i), coord1(i,:)
        write(445,*) atomsymbol_c(i), coord1(i,:)
    end do
    
    do i=1, Natom, 1
!        write(445,'(2A,3F11.9)') atomsymbol_c(i), coord2(i,:)
        write(445,*) atomsymbol_c(i), coord2(i,:)
    end do
    
    write(445,*) ''
    write(445,*) ''
    
    close(445)
    
    end subroutine writedimercom
    
!======================================================================================
!======================================================================================
! write *.dynput

    subroutine writelog2dynput(Natom,atomsymbol,coord1,filein,length, LDebug)    
    implicit none
    
    integer :: Natom, i, pos, length
    integer :: atomsymbol(Natom)
    real(kind=8) :: coord1(Natom,3)
    logical :: LDebug
    character(length) :: filein
    character(length+3) :: fileout
    
    
    pos=index(filein,'.')
    fileout(1:pos) = filein(1:pos)
    fileout(pos+1:length+3) = 'dynput'

    
	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine writelog2dynput'
	endif        
    
    open(unit=446,file=fileout)
    
    write(446,*) fileout(1:pos-1)//'m_thermal_5ps_1600K'
    write(446,*) 'MOPDYN FOR RATES'
    write(446,*) 'RUN  /home/youxq/you/mopac/MOPAC2012.exe'
    write(446,*) 'PARAMETERIZATION PM6'
    write(446,*) 'uhf bonds 1SCF SCFCRT=1.D-10 AUX(PRECISION=9)'
    write(446,*) '    11  F   T   T   F                   ! IPOT,LENER,LMOMEN,LXROT,LBOUND'
    write(446,*) '        F   F   F   T                   ! LTEMP,LTEST,RIMAGE,LCHARG'
    write(446,*) '        F   F   F   T                   ! LNOISE,LTRAMP,LQUENCH,LTSTAT'
    write(446,*) Natom,'  0    0                         ! N,NIMAGE,NRECT'
    write(446,*) '    1    .00001   .00000   1600.0   1   ! ICF,VMAX,ZNOISE,TREF,IDENRC'
    write(446,*) '    0.5     5       2000                ! DT,NTIME,NSNAP'
    write(446,*) '   0.0 0.0 0.0                          ! SX,SY,SZ'
    write(446,*) '     .000000         F                  ! TRAMP,RTSTAT'
    write(446,*) '     1000.00     100.000000             ! TCONST,ZTHERM'
    write(446,*) '       0.544  20.0                      ! SCALE DBOUND'  

    do i=1, Natom, 1
        write(446,*) coord1(i,:), atomsymbol(i), '3   0'
    end do
    
    write(446,*) ''
    write(446,*) ''
    
    close(446)
    
    end subroutine writelog2dynput
    
!==================================================================================
!==================================================================================
! 
    subroutine get_monomer_number_xyz(Natom,filename,length,LDebug)
    implicit none
    
    integer :: Natom, length
    character(length) :: filename
    logical :: LDebug
    
    open(unit=447,file=filename)
    read(447,*) Natom
    close(447)
                
	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine get_monomer_number_xyz'
		write(*,*) Natom
	endif                    
                
    end subroutine get_monomer_number_xyz            
                
!====================================================================================
!====================================================================================
! read the Nth frame in *.xyz

    subroutine readxyzNthframe(Nthframe,Natom,atomsymbol_c,coord1,filename, length,LDebug)
    implicit none
    
    integer :: Nthframe, Natom,i,length    
    real(kind=8) :: coord1(Natom,3)
    character :: atomsymbol_c(Natom)
    character(length) :: filename
    logical :: LDebug
    
    open(unit=447,file=filename)
    
    do i=1, (Nthframe-1)*(Natom+2), 1
        read(447,*)
    end do
    
    read(447,*)
    read(447,*)
    do i=1, Natom, 1
        read(447,*) atomsymbol_c(i), coord1(i,:)
    end do
    
    close(447)
                
    end subroutine readxyzNthframe
    
!===================================================================================
!===================================================================================
! write *.com

    subroutine writexyz2com(Nthframe,Natom,atomsymbol_c,coord1,filein,length, LDebug)
    implicit none
    
    integer :: Natom, i, length, Nthframe,pos
    character :: atomsymbol_c(Natom)
    real(kind=8) :: coord1(Natom,3)
    character(length) :: filein
    character(length+6) :: fileout
    character(5) :: num
    logical :: LDebug
                
    pos=index(filein,'.')
    fileout(1:pos-1) = filein(1:pos-1)
    write(num,'(i5)') Nthframe
    write(*,*) num
    fileout(pos:length+1+len(trim(adjustl(num)))) = '_'//trim(adjustl(num))//'.com'
    fileout(length+2+len(trim(adjustl(num))):length+6) =''
    
	If ( LDebug .eqv. .True. ) then
		write(*,*) '***********'
		write(*,*) 'subroutine writexyz2com'
		write(*,*) trim(adjustl(fileout))
	endif    
    
    open(unit=448,file=trim(adjustl(fileout)))
    write(448,*) '%nprocshared=12'
    write(448,*) '%mem=32GB'
    write(448,*) '# opt upm6'
    write(448,*) ''
    write(448,*) 'the '//num//'th frame in '//filein
    write(448,*) ''
    write(448,*) '0 1'
    
    do i=1, Natom, 1
!        write(445,'(2A,3F11.9)') atomsymbol_c(i), coord1(i,:)
        write(448,*) atomsymbol_c(i), coord1(i,:)
    end do
    
    
    write(448,*) ''
    write(448,*) ''
    
    close(448)     
    
    end subroutine writexyz2com           
                
                
                
                    
    
    
    
