! This is a program to assist a MD package MOPDYN13
! Hong-Bo Zhang, Oct 2013
! usage: mopdynapp 1 12  
!        mopdynapp 2 12 


	program mopdynapp
	implicit none
	
	integer :: argue_count, argue_int,argue2_int
	character(len=256) :: argue, argue2

! read the command
	argue_count = command_argument_count()
	if ( (argue_count .ne. 1) .and. (argue_count .ne. 2) ) then
		write(*,*) 'Number of command argument should be one or two at present'
		stop
	end if
	argue = ''
	call get_command_argument(1,argue)
	read(argue,'(i2)') argue_int
    call get_command_argument(2,argue2)
    read(argue2,'(i2)') argue2_int

! if ./mopdynapp 0
	if ( argue_int .eq. 0) then
		write(*,*) 'Do 1, 2 and 3 in readme.txt'
		write(*,*) 'Out of date!!! use "plot" instead of "mopdynapp 0" '
		
! if ./mopdynapp 1, do 6 in readme.txt
	elseif ( argue_int .eq. 1 ) then
	    write(*,*) 'ColDynGen_',trim(adjustl(argue2)),'.ini'
	    write(*,*) 'Do 6 in readm.txt'
		write(*,*) 'generating *.dynput for collision'
		write(*,*) ''
		call collision_dynput_generate(argue2_int)
		write(*,*) ''
		write(*,*) '****************'
		write(*,*) '6 in readme.txt has been don'
		write(*,*) ''

		
! if ./mopdynapp 2, do 4 and 5 in readme.txt
	elseif ( argue_int .eq. 2) then
	    write(*,*) 'collision_',trim(adjustl(argue2)),'.coord/veloc/accel/bond'
		write(*,*) 'in developing. Do 4 and 5 in readme.txt'
		write(*,*) 'dynamics of monomers in dimerization'
		write(*,*) 'bond orders'
		call MassCenterDyn(argue2_int)
		write(*,*) ''
		write(*,*) '*****************'

		
! if ./mopdynapp 3, what about 7 8 and 9 in readme.txt
	elseif (argue_int .eq. 3 ) then
		write(*,*) 'The 7, 8 and 9 in readme.txt is still in developing.'

		
	else
		write(*,*) 'only 0, 1, 2 and 3 are allowed.'
	endif

	write(*,*) "Done"

	end program
