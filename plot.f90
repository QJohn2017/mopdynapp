!==================================================================
!==================================================================
! working directory: ./xx/, which is the same path as collision_xx.xxxx


    program mdplot
    implicit none

	integer :: argue_count, argue2_int, n, n1,n2, num
	character :: argue
	character(2) :: argue2, num_c
	character(256) :: temp, prefix, aline, paramline
    integer :: ind,i,NTIME,NSNAP, timeinps
    real(kind=8) :: DT, simultime,param(11)

! read the command
	argue_count = command_argument_count()
	if ( argue_count .ne. 1 ) then
		write(*,*) 'usage: mdplot [m|d]'
		stop
	end if
	argue = ''
	call get_command_argument(1,argue)
	if ( (argue .ne. 'm') .and. (argue .ne. 'd') ) stop 'usage: mdplot [m|d]'
!    call get_command_argument(2,argue2)
!    read(argue2,'(i2)') argue2_int  
    
    call system('mkdir ./result/')
    temp=''
    prefix=''
    call system('ls *.dynput > prefix.dat')
    open(unit=444,file='prefix.dat')
    read(444,'(A)') temp
    close(444)
    call system('rm prefix.dat')
    n = scan(temp,'.')
    prefix = adjustl(temp(1:n-1))
    
    open(unit=444,file=trim(adjustl(prefix))//'.dynput')
    read(444,'(A)') aline
    do i=1,9, 1
        read(444,*)
    end do
    read(444,*) DT,NTIME
    close(444)
    
    call substituteString(aline,'_','-')
    
    call system( 'wc -l '//trim(adjustl(prefix))//'.itinf > Nsnap.temp' )
    open(unit=446,file='Nsnap.temp')
    read(446,*) NSNAP
    close(446)
    call system( 'rm Nsnap.temp')
    
    simultime = DT*NTIME*NSNAP
    timeinps = floor(simultime/1000.0)+1
    
    if( argue .eq. 'm' ) then
        call plotTemp(trim(adjustl(prefix)),trim(adjustl(aline)),timeinps)
        call plotKEPETE(trim(adjustl(prefix)),trim(adjustl(aline)),timeinps)
    elseif ( argue .eq. 'd' ) then    
        n1=scan(prefix,'_')
        if ( len(trim(prefix))-n1 .gt. 3) stop 'ColDynGen_xx.ini, where xx should be smaller than 100'
        num_c = prefix(n1+1:len(trim(prefix)))
!        read(num_c,*) num
        open(unit=446,file='ColDynGen_'//trim(adjustl(num_c))//'.ini')
            read(446,*)
            do i=1, 11, 1
                read(446,*) param(i)
            end do
        close(446)
        write(paramline,'(I4,"|",3I3,"|",3I3,"|",I2,I1,I3,"|",I2)') int(param(:))
        call plotTemp_D(trim(adjustl(prefix)),trim(adjustl(aline)),timeinps,trim(adjustl(paramline)))
        call plotKEPETE_D(trim(adjustl(prefix)),trim(adjustl(aline)),timeinps,trim(adjustl(paramline)))
        call plotBondOrder(trim(adjustl(prefix)),trim(adjustl(aline)),timeinps,trim(adjustl(paramline)))
        call plotMassCenter_FTF(trim(adjustl(prefix)),trim(adjustl(aline)),timeinps,trim(adjustl(paramline))) ! face to face, b = 0
!        call plotMassCenter3D(trim(adjustl(prefix)),trim(adjustl(aline)),timeinps)
    end if    
    
    
    end program
!===================================================================
!===================================================================
! find b in a, and substitute b by c, where b and c are character(len=1)
! easy to change b and c with arbitrary length
!         if ( a(i:len(b) .eq. b ) a(i:len(b))=c
    subroutine substituteString(a,b,c)
    implicit none
    
    character(len=*) :: a
    character :: b, c
    integer :: i
    
    do i=1, len(a), 1
        if ( a(i:i) .eq. b ) a(i:i)=c
    end do
    
    end subroutine substituteString    
!===================================================================
!===================================================================
! plot temperature of the whole system
    subroutine plotTemp(prefix,title,timeinps)
    implicit none
    
    character(256) :: filename
    character(len=*) :: prefix, title
    integer :: n, ind,i,NTIME,NSNAP, timeinps
    real(kind=8) :: DT, simultime
    
!    n = scan(prefix,'_')
!    if ( n .ne. 0) then
        
    
    filename = './result/'//trim(adjustl(prefix))//'_Tem.plot'
    open( unit=445,file=trim(adjustl(filename)) )
    
    write(445,*) "set term post enh color"
    write(445,*) "set out '"//trim(adjustl(prefix))//"_Tem.eps'"
    write(445,*) "set key bottom right"
    write(445,*) "set key box"
    write(445,*) "set title '"//trim(adjustl(title))//" Temperature'"
    write(445,*) "set autoscale"
    write(445,*) "set encoding iso_8859_1"
!    write(445,*) "set xtics nomirror"
!    write(445,*) "set ytics nomirror"
!    write(445,*) "set y2tics"
!    write(445,*) "set ytics 0.05"
!    write(445,*) "set y2tics 0.05"
!    write(445,*) "set y2range [0:0.55]"
    write(445,*) "unset xlabel"
    write(445,*) "unset ylabel"
    write(445,*) "unset y2label"    
    write(445,*) "unset xtics"
    write(445,*) "unset y2tics"
    write(445,*) "set xlabel '{/Symbol t} ({\327}} 10^3 fs)'"
    write(445,*) "set ylabel 'Temperature (Kelvin)'"
    write(445,*) "set xtics 0, 1,", timeinps
    write(445,*) "plot '../"//trim(adjustl(prefix))//".itinf' u ($5/1000):4 '%lf,%lf,%lf,%lf,%lf' t ' T' w d linewidth 5"
    
!    plot 'r51.itinf' u 5:2 "%lf,%lf,%lf,%lf,%lf" t ' P.E.' w d linewidth 5, \
!                    'r51.itinf' u 5:1 "%lf,%lf,%lf,%lf,%lf" t ' K.E.' w d linewidth 5 axes x1y2, \
!                	'r51.itinf' u 5:3 "%lf,%lf,%lf,%lf,%lf" t ' T.E.' w d linewidth 5    
    
    close(445)
    end subroutine plotTemp    
    
!==============================================================================
!==============================================================================
! plot Kinetic energy; potential energy and total energy of the whole system
    subroutine plotKEPETE(prefix,title,timeinps)
    implicit none
    
    character(256) :: filename
    character(len=*) :: prefix, title
    integer :: timeinps
    
    filename = './result/'//trim(adjustl(prefix))//'_KPT.plot'
    open(unit=445, file=trim(adjustl(filename)))

    write(445,*) "set term post enh color"
    write(445,*) "set out '"//trim(adjustl(prefix))//"_KPT.eps'"
    write(445,*) "set key bottom right"
    write(445,*) "set key box"
    write(445,*) "set title '"//trim(adjustl(title))//" KE PE TE'"
    write(445,*) "set autoscale"
    write(445,*) "set encoding iso_8859_1"
    write(445,*) "unset xlabel"
    write(445,*) "unset ylabel"
    write(445,*) "unset y2label"    
    write(445,*) "unset xtics"
!    write(445,*) "unset ytics"
    write(445,*) "set xtics nomirror"
    write(445,*) "set ytics nomirror"
    write(445,*) "set y2tics"
!    write(445,*) "set ytics 0.05"
!    write(445,*) "set y2tics 0.05"
!    write(445,*) set y2range [0:0.55]
    write(445,*) "set xlabel '{/Symbol t} ({\327}} 10^3 fs)'"
    write(445,*) "set ylabel 'Potential Energy (Hartrees)'"
    write(445,*) "set y2label 'Kinetic Energy (Hartrees)'"
    write(445,*) "set xtics 0, 1,", timeinps
    write(445,*) "plot '../"//trim(adjustl(prefix))//".itinf' u ($5/1000):2 '%lf,%lf,%lf,%lf,%lf' t ' P.E.' w d lw 5, \"
    write(445,*) "    '../"//trim(adjustl(prefix))//".itinf' u ($5/1000):1 '%lf,%lf,%lf,%lf,%lf' t ' K.E.' w d lw 5 axes x1y2, \"
    write(445,*) "    '../"//trim(adjustl(prefix))//".itinf' u ($5/1000):3 '%lf,%lf,%lf,%lf,%lf' t ' T.E.' w d lw 5"
    	
    close(445)
    end subroutine plotKEPETE
!===================================================================
!===================================================================
! plot temperature of the whole system
    subroutine plotTemp_D(prefix,title,timeinps,paramline)
    implicit none
    
    character(256) :: filename
    character(len=*) :: prefix, title,paramline
    integer :: n, ind,i,NTIME,NSNAP, timeinps
    real(kind=8) :: DT, simultime
    
!    n = scan(prefix,'_')
!    if ( n .ne. 0) then
        
    
    filename = './result/'//trim(adjustl(prefix))//'_Tem.plot'
    open( unit=445,file=trim(adjustl(filename)) )
    
    write(445,*) "set term post enh color"
    write(445,*) "set out '"//trim(adjustl(prefix))//"_Tem.eps'"
    write(445,*) "set key bottom right"
    write(445,*) "set key box"
    write(445,*) "set title '"//trim(adjustl(title))//" Temperature "//trim(adjustl(paramline))//"'"
    write(445,*) "set autoscale"
    write(445,*) "set encoding iso_8859_1"
!    write(445,*) "set xtics nomirror"
!    write(445,*) "set ytics nomirror"
!    write(445,*) "set y2tics"
!    write(445,*) "set ytics 0.05"
!    write(445,*) "set y2tics 0.05"
!    write(445,*) "set y2range [0:0.55]"
    write(445,*) "unset xlabel"
    write(445,*) "unset ylabel"
    write(445,*) "unset y2label"    
    write(445,*) "unset xtics"
    write(445,*) "unset y2tics"
    write(445,*) "set xlabel '{/Symbol t} ({\327}} 10^3 fs)'"
    write(445,*) "set ylabel 'Temperature (Kelvin)'"
    write(445,*) "set xtics 0, 1,", timeinps
    write(445,*) "plot '../"//trim(adjustl(prefix))//".itinf' u ($5/1000):4 '%lf,%lf,%lf,%lf,%lf' t ' T' w d linewidth 5"
    
!    plot 'r51.itinf' u 5:2 "%lf,%lf,%lf,%lf,%lf" t ' P.E.' w d linewidth 5, \
!                    'r51.itinf' u 5:1 "%lf,%lf,%lf,%lf,%lf" t ' K.E.' w d linewidth 5 axes x1y2, \
!                	'r51.itinf' u 5:3 "%lf,%lf,%lf,%lf,%lf" t ' T.E.' w d linewidth 5    
    
    close(445)
    end subroutine plotTemp_D    
    
!==============================================================================
!==============================================================================
! plot Kinetic energy; potential energy and total energy of the whole system
    subroutine plotKEPETE_D(prefix,title,timeinps,paramline)
    implicit none
    
    character(256) :: filename
    character(len=*) :: prefix, title,paramline
    integer :: timeinps
    
    filename = './result/'//trim(adjustl(prefix))//'_KPT.plot'
    open(unit=445, file=trim(adjustl(filename)))

    write(445,*) "set term post enh color"
    write(445,*) "set out '"//trim(adjustl(prefix))//"_KPT.eps'"
    write(445,*) "set key bottom right"
    write(445,*) "set key box"
    write(445,*) "set title '"//trim(adjustl(title))//" KE PE TE "//trim(adjustl(paramline))//"'"
    write(445,*) "set autoscale"
    write(445,*) "set encoding iso_8859_1"
    write(445,*) "unset xlabel"
    write(445,*) "unset ylabel"
    write(445,*) "unset y2label"    
    write(445,*) "unset xtics"
!    write(445,*) "unset ytics"
    write(445,*) "set xtics nomirror"
    write(445,*) "set ytics nomirror"
    write(445,*) "set y2tics"
!    write(445,*) "set ytics 0.05"
!    write(445,*) "set y2tics 0.05"
!    write(445,*) set y2range [0:0.55]
    write(445,*) "set xlabel '{/Symbol t} ({\327}} 10^3 fs)'"
    write(445,*) "set ylabel 'Potential Energy (Hartrees)'"
    write(445,*) "set y2label 'Kinetic Energy (Hartrees)'"
    write(445,*) "set xtics 0, 1,", timeinps
    write(445,*) "plot '../"//trim(adjustl(prefix))//".itinf' u ($5/1000):2 '%lf,%lf,%lf,%lf,%lf' t ' P.E.' w d lw 5, \"
    write(445,*) "    '../"//trim(adjustl(prefix))//".itinf' u ($5/1000):1 '%lf,%lf,%lf,%lf,%lf' t ' K.E.' w d lw 5 axes x1y2, \"
    write(445,*) "    '../"//trim(adjustl(prefix))//".itinf' u ($5/1000):3 '%lf,%lf,%lf,%lf,%lf' t ' T.E.' w d lw 5"
    	
    close(445)
    end subroutine plotKEPETE_D

!========================================================================
!========================================================================
!  plot bond order. BOB: Inter-monomers Bond Order
    subroutine plotBondOrder(prefix,title,timeinps,paramline) 
    implicit none
    
    character(256) :: filename
    character(len=*) :: prefix, title,paramline
    integer :: timeinps
    
    filename = './result/'//trim(adjustl(prefix))//'_IBO.plot'   
    open(unit=445,file=trim(adjustl(filename)))
    
    write(445,*) "set term post enh color"
    write(445,*) "set out '"//trim(adjustl(prefix))//"_IBO.eps'"
    write(445,*) "set key bottom right"
    write(445,*) "set key box"
    write(445,*) "set title '"//trim(adjustl(title))//" Inter-monomers Bond Order "//trim(adjustl(paramline))//"'"
    write(445,*) "set autoscale"
    write(445,*) "set encoding iso_8859_1"
    write(445,*) "unset xlabel"
    write(445,*) "unset ylabel"
    write(445,*) "unset y2label"    
    write(445,*) "unset xtics"
!    write(445,*) "unset ytics"
    write(445,*) "set xlabel '{/Symbol t} ({\327}} 10^3 fs)'"
    write(445,*) "set ylabel 'Bond Order'"
    write(445,*) "set xtics 0, 1,", timeinps
    write(445,*) "plot '"//trim(adjustl(prefix))//".ibo' u ($2/1000):3 t ' IBO' w d lw 5"
    
    close(445)
    end subroutine plotBondOrder
    
!==========================================================================
!==========================================================================
! plot Center of Mass Coordinate(CMC), Center of Mass Velocity(CMV) and Center of Mass Acceleration(CMA)

    subroutine plotMassCenter_FTF(prefix,title,timeinps,paramline)
    implicit none
    
    character(256) :: filename
    character(len=*) :: prefix, title,paramline
    integer :: timeinps

!----------------------------    
    filename = './result/'//trim(adjustl(prefix))//'_CMC.plot'   
    open(unit=445,file=trim(adjustl(filename)))    
    
    write(445,*) "set term post enh color"
    write(445,*) "set out '"//trim(adjustl(prefix))//"_CMC.eps'"
    write(445,*) "set key bottom right"
    write(445,*) "set key box"
    write(445,*) "set title '"//trim(adjustl(title))//" Center of Mass Coordinates "//trim(adjustl(paramline))//"'"
    write(445,*) "set autoscale"
    write(445,*) "set encoding iso_8859_1"
    write(445,*) "unset xlabel"
    write(445,*) "unset ylabel"
    write(445,*) "unset y2label"    
    write(445,*) "unset xtics"
!    write(445,*) "unset ytics"
    write(445,*) "set xlabel '{/Symbol t} ({\327}} 10^3 fs)'"
    write(445,*) "set ylabel 'z-coordinates of CM {\305}'"
    write(445,*) "set xtics 0, 1,", timeinps
    write(445,*) "plot '"//trim(adjustl(prefix))//".cm' u ($2/1000):5 t ' CMC1' w l lw 5, \"
    write(445,*) "    '"//trim(adjustl(prefix))//".cm' u ($2/1000):8  t ' CMC2' w l lw 5"    
    
    close(445)
! ---------------------------    
    filename = './result/'//trim(adjustl(prefix))//'_CMV.plot'   
    open(unit=445,file=trim(adjustl(filename)))    
    
    write(445,*) "set term post enh color"
    write(445,*) "set out '"//trim(adjustl(prefix))//"_CMV.eps'"
    write(445,*) "set key bottom right"
    write(445,*) "set key box"
    write(445,*) "set title '"//trim(adjustl(title))//" Center of Mass Velocities "//trim(adjustl(paramline))//"'"
    write(445,*) "set autoscale"
    write(445,*) "set encoding iso_8859_1"
    write(445,*) "unset xlabel"
    write(445,*) "unset ylabel"
    write(445,*) "unset y2label"    
    write(445,*) "unset xtics"
!    write(445,*) "unset ytics"
    write(445,*) "set xlabel '{/Symbol t} ({\327}} 10^3 fs)'"
    write(445,*) "set ylabel 'v_z of CM {\305}/fs^2'"
    write(445,*) "set xtics 0, 1,", timeinps
    write(445,*) "plot '"//trim(adjustl(prefix))//".cm' u ($2/1000):11 t ' CMV1' w l lw 5, \"
    write(445,*) "    '"//trim(adjustl(prefix))//".cm' u ($2/1000):15  t ' CMV2' w l lw 5"       
    
    close(445)    
!---------------------------    
    filename = './result/'//trim(adjustl(prefix))//'_CMA.plot'   
    open(unit=445,file=trim(adjustl(filename)))    
    
    write(445,*) "set term post enh color"
    write(445,*) "set out '"//trim(adjustl(prefix))//"_CMA.eps'"
    write(445,*) "set key bottom right"
    write(445,*) "set key box"
    write(445,*) "set title '"//trim(adjustl(title))//" Center of Mass Accelerations "//trim(adjustl(paramline))//"'"
    write(445,*) "set autoscale"
    write(445,*) "set encoding iso_8859_1"
    write(445,*) "unset xlabel"
    write(445,*) "unset ylabel"
    write(445,*) "unset y2label"    
    write(445,*) "unset xtics"
!    write(445,*) "unset ytics"
    write(445,*) "set xlabel '{/Symbol t} ({\327}} 10^3 fs)'"
    write(445,*) "set ylabel 'a_z of CM {\305}/fs^2'"
    write(445,*) "set xtics 0, 1,", timeinps
    write(445,*) "plot '"//trim(adjustl(prefix))//".cm' u ($2/1000):19 t ' CMA1' w l lw 5, \"
    write(445,*) "    '"//trim(adjustl(prefix))//".cm' u ($2/1000):23  t ' CMA2' w l lw 5"           
    
    close(445)    
!--------------------------
    filename = './result/'//trim(adjustl(prefix))//'_CM.plot'   
    open(unit=445,file=trim(adjustl(filename)))    
    
    write(445,*) "set term post enh color"
    write(445,*) "set out '"//trim(adjustl(prefix))//"_CM2D.eps'"
    write(445,*) "set key bottom right"
    write(445,*) "set key box"
    write(445,*) "set title '"//trim(adjustl(title))//" CM 2D XVA "//trim(adjustl(paramline))//"'"
    write(445,*) "set autoscale"
    write(445,*) "set encoding iso_8859_1"
    write(445,*) "unset xlabel"
    write(445,*) "unset ylabel"
    write(445,*) "unset y2label"    
    write(445,*) "unset xtics"
!    write(445,*) "unset ytics"
    write(445,*) "set xlabel '{/Symbol t} ({\327}} 10^3 fs)'"
    write(445,*) "set ylabel 'z-coordinates {\305}; v_z*100{\305}/fs; a_z*1000{\305}/fs^2'"
    write(445,*) "set xtics 0, 1,", timeinps
    write(445,*) "plot '"//trim(adjustl(prefix))//".cm' u ($2/1000):5 t ' CMC1' w l lw 5 axes x1y1, \"
!    write(445,*) "     '"//trim(adjustl(prefix))//".cm' u ($2/1000):8  t ' CMC2' w l lw 5 axes x1y1, \"    
    write(445,*) "     '"//trim(adjustl(prefix))//".cm' u ($2/1000):($11*100) t ' CMV1' w l lw 5, \"
!    write(445,*) "     '"//trim(adjustl(prefix))//".cm' u ($2/1000):($15*100)  t ' CMV2' w l lw 5, \"       
    write(445,*) "     '"//trim(adjustl(prefix))//".cm' u ($2/1000):($19*1000) t ' CMA1' w l lw 5" !, \"
!    write(445,*) "     '"//trim(adjustl(prefix))//".cm' u ($2/1000):($23*1000)  t ' CMA2' w l lw 5"   

    close(445)        
    
    end subroutine plotMassCenter_FTF    
    
    
