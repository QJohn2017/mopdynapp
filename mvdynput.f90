! move ./relax/collision_xx.dynput ColDynGen_xx.ini to ./xx/
! generate job in which the last line should be mopdyn13 collision_xx.dynput
! working directory is the parent directory of relax/ and xx/
! ./mvdynput

    program mvdynput
    implicit none
    
    integer :: lines,i,n1,n2
    integer,allocatable :: num(:)
!    character(2) :: temp
    character(2),allocatable :: num_c(:)
    character(256),allocatable :: filenames(:)
!    character(256) :: tempfilename
    
    call system('find ./relax/ -name "ColDynGen_*.ini" > inifilename.dat')
    call system('wc -l inifilename.dat > inifilenamelength.dat')
    
    open(unit=444,file='inifilenamelength.dat')
    read(444,*) lines
    close(444)
    call system('rm inifilenamelength.dat')
    
    write(*,*) 'there are',lines,'ColDynGen_xx.ini files in ./relax/'
    
    allocate( num(lines), num_c(lines), filenames(lines) )
    num=0
    num_c=''
    
    open(unit=445,file='inifilename.dat')
    
    do i=1, lines, 1
        read(445,'(A)') filenames(i)
        write(*,*) trim(adjustl(filenames(i)))
    end do
    
    close(445)
    
    do i=1, lines, 1
        n1=scan(filenames(i),'_')
        n2=scan(filenames(i),'.',.True.)
!        write(*,*) n1,n2
        if ( n2-n1 .gt. 3) stop 'ColDynGen_xx.ini, where xx should be smaller than 100'
        num_c(i) = filenames(i)(n1+1:n2-1)
        read(num_c(i),*) num(i)
!        write(*,*) num(i), num_c(i)
        call system('mkdir ./'//trim(adjustl(num_c(i)))//'/')
        call system('mkdir ./'//trim(adjustl(num_c(i)))//'/Debug/')
        write(*,*) 'make directory: ./'//trim(adjustl(num_c(i)))//'/'
        call system('mv ./relax/ColDynGen_'//trim(adjustl(num_c(i)))//'.ini  ./'//trim(adjustl(num_c(i)))//'/')
        write(*,*) 'move ini and dynput files'
        call system('mv ./relax/collision_'//trim(adjustl(num_c(i)))//'.dynput  ./'//trim(adjustl(num_c(i)))//'/')
        call generate_job(num_c(i))
        call system('mv job'//trim(adjustl(num_c(i)))//' ./'//trim(adjustl(num_c(i)))//'/')
        write(*,*) 'generate job file and move to ./'//trim(adjustl(num_c(i)))//'/'
        call system('mv ./relax/Debug/collision_initial_'//trim(adjustl(num_c(i)))//'.xyz ./'//trim(adjustl(num_c(i)))//'/Debug/')
        write(*,*) 'move collision_initial.xyz files'
    end do
    
    end program
    
!=================================================================================
!=================================================================================
!  generate job, a bsub script

    subroutine generate_job(num_c)
    implicit none
    
    character(2) :: num_c
    
    open(unit=446,file='job'//trim(adjustl(num_c)))    
    
    write(446,*) '#BSUB -q hpc_linux'
    write(446,*) '#BSUB -a openmp'
    write(446,*) '#BSUB -n 12'
    write(446,*) '#BSUB -R "span[ptile=12]"'
    write(446,*) '#BSUB -o output.%J'
    write(446,*) '#BSUB -e error.%J'
    write(446,*) 'mopdyn13 collision_'//trim(adjustl(num_c))//'.dynput'
    
    close(446)
    
    end subroutine generate_job
    
    
    
