! a script to run:
!    mopdynapp 1 xx
! if there is a ColDynGen_xx.ini in folder ./relax/
! working directory is the parent directory of relax/ and xx/
! ./dynputscript

    program dynputscript
    implicit none
    
    integer :: lines,i,n1,n2
    integer,allocatable :: num(:)
    character(2),allocatable :: num_c(:)
    character(256),allocatable :: filenames(:)
    character(256) :: pathmopdynapp

    pathmopdynapp = '/home/youxq/mopdynapp/'
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
        call system('')
!        call system('which mopdynapp > pathmopdynapp.dat')
        call system('cd ./relax/ ; '//trim(adjustl(pathmopdynapp))//'mopdynapp 1 '//num_c(i))
    end do
    
    end program
