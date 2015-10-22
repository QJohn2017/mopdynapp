! auto generate ColDynGen_xx.ini from *.sample file.
! there are two kinds of *.sample.
! the first line of *.sample can be 'Debug', 'Param' or 'Latin'

!! if 'Latin', the program will read in a matrix, in which the # of row is # of sample, 
!!   and # of column is # of parameters
!! Latin
!! 11
!! {1st set}
!! {2nd set}
!! ... 
!! {num_ini th set}
!! NOTICE: the second line should be 11. and # of {set} = 9, where Temperature and tranVectz
!!         are fixed in Latin. Modified subroutine readsamplelatin if needed.

!! if 'Debug', parameters vary independently.
!!   there will be an 'F' or an 'S' at the beginning of a line. 'F' means this parameter
!!   will be fixed, 'S' means this parameter will be scanned.
!!   F p
!!   S p_start step_length data_number
!!   UGLY!!! I don't know how to deal with variable loop number. There should be pn loops, where pn is a variable.

!! if 'Param', parameters vary synchronously.
!!  parameterized sample file is used
!! there will be an 'F' or an 'S' at the beginning of a line. 'F' means this parameter will be fixed,
!! 'S' means this parameter will be scanned.
!!  F 15
!!  S 3 7
!!  p 0 4 10
!!  F 15 means this parameter is fixed, and its value is 15
!!  p 0 4 10 means p_start, p_step_length, p_step_num. In this example,
!!     the parameter p \in {0, 4, 8, 12, 16, 20, 24, 28, 32, 36}
!!  S 3 7 means this parameter will be scanned, and its value is 3+7*p

    program autoini
    implicit none
    
    character(5) :: atype
    character(256) :: filename
    integer :: argue_count, n, pn, i,j,stepn, num_ini, scannum ! parameter_number, # of S
    integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,i10,i11
    real(kind=8), allocatable :: parameters(:), steplength(:), paramwrite(:,:),coef(:)
    integer,allocatable :: stepnum(:), iscan(:) 
    real(kind=8) :: startp, stepl, par
    character,allocatable :: options(:)
    
	argue_count = command_argument_count()
	if ( argue_count .ne. 1 ) then
		write(*,*) 'usage: autoini *.sample'
		stop
	end if
	filename = ''
	call get_command_argument(1,filename)
	
	n=scan(filename,'.')
	if ( filename(n+1:n+6) .ne. 'sample' ) stop 'usage: autoini *.sample'
	
	write(*,*) trim(adjustl(filename))
	
	
	open(unit=444,file=trim(adjustl(filename)))
	
	read(444,*) atype
!    write(*,*) atype
	if ( atype .eq. 'Debug' ) then
	    write(*,*) 'sample mode: Debug'
	    call sampledebuglength(trim(adjustl(filename)),pn)
        if ( pn .ne. 11 ) write(*,*) 'Warning! there should be 11 parameters'
	    allocate(parameters(pn), steplength(pn), stepnum(pn))
	    call readsampledebug(trim(adjustl(filename)), pn, parameters, steplength, stepnum)
!        write(*,*) stepnum(8)
        do i=1, pn, 1
            write(*,*) parameters(i), steplength(i), stepnum(i)
        end do
        num_ini = 1
        do i=1, pn, 1
            num_ini = num_ini*stepnum(i)
        end do
        
!        scannum = 0
!        do i=1, pn, 1
!            if( stepnum(i) .ne. 1 ) then
!                scannum = scannum + 1
!            end if
!        end do
!        allocate( iscan(scannum) )
!        
!        j=1
!        do i=1,pn,  1
!            if( stepnum(i) .ne. 1) iscan(j) = i
!            j = j+1
!        end do
!        
!        do i=1, scannum, 1
!            do j=1, stepnum(iscan(i)), 1
!                paramwrite(iscan(i),(j-1)*num_ini/stepnum(iscan(i))+1:j*num_ini/stepnum(iscan(i))) &
!                &  = parameters(iscan(i)) + (j-1)*steplength(iscan(i))
!            end do
!        end do        
        
        
        allocate(paramwrite(num_ini,pn))
        paramwrite = 0.0
!        do i=1, pn, 1
!            if( stepnum(i) .eq. 1 ) paramwrite(i,:) = parameters(i)
!        end do
        j=1
        
        do i1=1, stepnum(1), 1
        do i2=1, stepnum(2), 1
        do i3=1, stepnum(3), 1
        do i4=1, stepnum(4), 1
        do i5=1, stepnum(5), 1
        do i6=1, stepnum(6), 1
        do i7=1, stepnum(7), 1
        do i8=1, stepnum(8), 1
        do i9=1, stepnum(9), 1
        do i10=1, stepnum(10), 1
        do i11=1, stepnum(11), 1
            paramwrite(j,1) = parameters(1) + (i1-1)*steplength(1)
            paramwrite(j,2) = parameters(2) + (i2-1)*steplength(2)
            paramwrite(j,3) = parameters(3) + (i3-1)*steplength(3)
            paramwrite(j,4) = parameters(4) + (i4-1)*steplength(4)
            paramwrite(j,5) = parameters(5) + (i5-1)*steplength(5)
            paramwrite(j,6) = parameters(6) + (i6-1)*steplength(6)
            paramwrite(j,7) = parameters(7) + (i7-1)*steplength(7)
            paramwrite(j,8) = parameters(8) + (i8-1)*steplength(8)
            paramwrite(j,9) = parameters(9) + (i9-1)*steplength(9)
            paramwrite(j,10) = parameters(10) + (i10-1)*steplength(10)
            paramwrite(j,11) = parameters(11) + (i11-1)*steplength(11)
            j=j+1
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        
        do i=1, num_ini, 1
            write(*,'(11F7.1)') paramwrite(i,:)
        end do
        
       call writesamples(paramwrite,num_ini,pn)

    elseif ( atype .eq. 'Param' ) then

        write(*,*) 'Sample mode: Param'
   	    call sampledebuglength(trim(adjustl(filename)),pn)
   	    pn = pn-1
        if ( pn .ne. 11 ) write(*,*) 'Warning! there should be 11 parameters'
   	    allocate(parameters(pn), coef(pn), options(pn))
   	    call readsampleparam(trim(adjustl(filename)), pn, parameters,startp,stepl,stepn,coef,options)
   	    allocate(paramwrite(stepn,pn))
   	    paramwrite = 0.0
        do i=1, pn, 1
            if( options(i) .eq. 'F' ) paramwrite(:,i) = parameters(i)
            if( options(i) .eq. 'S' ) then
                do j=1,stepn , 1
                    par = startp + stepl*(j-1)
                    paramwrite(j,i) = parameters(i) + coef(i)*par
                end do
            end if
        end do   	    

        do i=1, stepn, 1
            write(*,'(11F7.1)') paramwrite(i,:)
        end do
        
       call writesamples(paramwrite,stepn,pn)   	    
   	       	    
	elseif ( atype .eq. 'Latin' ) then
	
	    write(*,*) 'Sample mode: Latin'
	    call samplelatinlength(trim(adjustl(filename)),pn,num_ini)
   	    allocate(paramwrite(num_ini,pn))
	    call readsamplelatin(paramwrite,pn,num_ini)
	    call EulAngCheck
	    call writesamples(paramwrite,num_ini,pn)
	end if
	close(444)


	
	end program
	
!=========================================================
!=========================================================
! 
    subroutine sampledebuglength(filename,pn)
    implicit none
    
    integer :: pn
    character(len=*) :: filename

    call system('wc -l '//filename//' > sampletemp.dat')
    open(unit=445,file='sampletemp.dat')
    read(445,*) pn
    close(445) 
    call system('rm sampletemp.dat')

    pn = pn - 1  
!    write(*,*) pn
    end subroutine sampledebuglength

!=========================================================
!=========================================================
!  
    subroutine readsampledebug(filename,pn,parameters,steplength,stepnum)  
    implicit none
    
    character(len=*) :: filename
    real(kind=8) :: parameters(pn), steplength(pn)
    integer :: stepnum(pn), pn,i
    character :: option, junk
    character(256) :: aline
    
    parameters = 0.0
    steplength = 0.0
    stepnum = 1
    
    do i=1, pn, 1
        read(444,'(A)') aline
        read(aline, *) option
        if ( option .eq. 'F') read(aline,*)  junk, parameters(i)
        if ( option .eq. 'S') read(aline,*) junk,parameters(i),steplength(i), stepnum(i)
    end do

!    do i=1, pn, 1
!        write(*,*) parameters(i),steplength(i), stepnum(i)
!    end do
 
    end subroutine readsampledebug

!=========================================================================================
!=========================================================================================
!
    subroutine writesamples(paramwrite,num_ini,pn)
    implicit none
    
    integer :: num_ini, pn, i,j
    real(kind=8) :: paramwrite(num_ini,pn)
    character(3) :: i_c

    do i=1, num_ini, 1
        write(i_c,'(I3)') i
        open(unit=i+499,file='ColDynGen_'//trim(adjustl(i_c))//'.ini')
        write(i+499,*) 'Debug'
        do j=1, pn, 1
            write(i+499,*) paramwrite(i,j)
        end do
        close(i+499)
    end do

    end subroutine writesamples

!=========================================================================================
!=========================================================================================
!
    subroutine readsampleparam(filename,pn,parameters,startp,stepl,stepn,coef,options)
    implicit none
    
    character(len=*) :: filename
    integer :: pn,i, stepn
    real(kind=8) :: parameters(pn), coef(pn)
    character :: option, junk, options(pn)
    character(256) :: aline
    real(kind=8) :: stepl,startp
    
    parameters = 0.0
    options = ''
    coef = 0.0
    
    do i=1, pn, 1
        read(444,'(A)') aline
        read(aline, *) option
        if ( option .eq. 'F') read(aline,*)  options(i), parameters(i)
        if ( option .eq. 'S') read(aline,*) options(i),parameters(i),coef(i)
    end do
    
    read(444,*) junk, startp,stepl, stepn

    end subroutine readsampleparam    

!====================================================================================
!====================================================================================
!
    subroutine samplelatinlength(filename,pn,num_ini)
    implicit none
    
    integer :: pn, num_ini, i
    character(len=*) :: filename

    call system('wc -l '//filename//' > sampletemp.dat')
    open(unit=445,file='sampletemp.dat')
    read(445,*) num_ini
    close(445) 
    call system('rm sampletemp.dat')
    num_ini = num_ini - 2
    read(444,*) pn

    end subroutine samplelatinlength

!=====================================================================================
!=====================================================================================
!
    subroutine readsamplelatin(paramwrite,pn,num_ini)
    implicit none
    
    integer :: pn, num_ini,i
    real(kind=8) :: paramwrite(num_ini,pn)
    
    paramwrite(:,1) = 0.5
    paramwrite(:,pn)= 0.5
    do i=1, num_ini, 1
        read(444,*) paramwrite(i,2:pn-1)
    end do
    paramwrite = 2*paramwrite - 1.0
    
    paramwrite(:,1) = paramwrite(:,1)+1600
    paramwrite(:,2) = paramwrite(:,2)*90
    paramwrite(:,3) = paramwrite(:,3)*90
    paramwrite(:,4) = paramwrite(:,4)*90 + 90
    paramwrite(:,5) = paramwrite(:,5)*90
    paramwrite(:,6) = paramwrite(:,6)*90
    paramwrite(:,7) = paramwrite(:,7)*90 + 90
    paramwrite(:,8) = paramwrite(:,8) !? find the reference point!!!
    paramwrite(:,9) = paramwrite(:,8) !?
    paramwrite(:,10) = paramwrite(:,8) !?
    paramwrite(:,11) = paramwrite(:,8) !?
    
    end subroutine
    
!!==================================================================================
!!==================================================================================
!
    subroutine EulAngCheck
    implicit none
    
    end subroutine    






