! bond orders
! the atom number should be an integer multiply of 3

    subroutine bondorder(Natom,Nsnap,Delta,atomsymbol,atomsymbol_c,arg,LDebug)
    implicit none   
    
    integer :: Natom, Nsnap, atomsymbol(Natom),arg,i,j,k,ii,iii,junki,cou
    real(kind=8) :: Delta
    character :: atomsymbol_c(Natom), junkc
    logical :: LDebug
    real(kind=8) :: order_all(Nsnap,Natom*2,Natom*2) !order_all(isnap, hang|row, lie|column)
    real(kind=8) :: totalorder(Nsnap)
    character(256) :: arg_c, filename, filename1
    
    write(arg_c,'(I2)') arg
    filename = 'collision_'//trim(adjustl(arg_c))//'.bond'    
    
! read bond order matrix    
    order_all = 0.0
    totalorder = 0.0
    call system("grep -E '([CH]    [0-9]|[CH]   [0-9][0-9])' "//trim(adjustl(filename))//" > temp.bond")

    open(unit=449,file="temp.bond")
    
    do i=1, Nsnap, 1 !isnap
    
        do j=1, Natom*2/6,1  ! Natom is an integer multiply of 3 in R51 case !column
            cou=0
            do k=1+(j-1)*6, Natom*2, 1 !row
                cou = cou + 1
                if( cou < 6 ) then
                    read(449,*) junkc,junki, order_all(i,k,((j-1)*6+1):((j-1)*6+cou))
                else
                    read(449,*)  junkc,junki, order_all(i,k,((j-1)*6+1):j*6)
                endif
            end do
        end do
        
        do ii=Natom+1, 2*Natom, 1
            do iii=1, Natom,1
                totalorder(i) = totalorder(i)+order_all(i,ii,iii)
            end do
        end do    

    end do
    
    close(449)
  
! e.g. order_all(47,14) is the bond order of atom15 and atom47       
! order_all(Natom+1~2*Natom,1~Natom): bond order between atom in M1 and atom in M2
! lower triangle of order_all(1~Natom, 1~Natom): bond order between two atoms in M1
! lower triangle of order_all(Natom+1~2*Natom,Natom+1~2*Natom): bond order between two atoms in M2

! for brief. we define b_ij : bond order between the ith atom in M1 and jth atom in M2
! that is to say i \in {1, Natom}, and j \in {Natom+1, 2*Natom}.
! totalorder = \sum_i \sum_j b_ij

! if one want to calculate the summation of bond order between ith atom in M1 and all atoms in M2
!   order(isnap,ith atom) = 0.0
!   do j=Natom+1, 2*Natom, 1
!       order(isnap, ith atom) = order(isnap, ith atom) + order_all(isnap, j , ith atom)
!   end do
! 
!   and actually totalorder is
!   totalorder(isnap) = 0.0
!   do i = 1 ,Natom, 1
!       totalorder(isnap) = totalorder(isnap) + order(isnap, i)
!   end do 

    filename1 = 'collision_'//trim(adjustl(arg_c))//'.ibo'    

    open(unit=450,file='./result/'//trim(adjustl(filename1)))
    
    do i=1, Nsnap, 1
        write(450,*) i, i*Delta, totalorder(i)
    end do
    
    close(450)
    
    open(unit=451,file='./result/collision_'//trim(adjustl(arg_c))//'.dob',position='APPEND')   
        write(451,*) 'time of forming bonds / fs ;      total simulation time / fs'
        write(451,*) count(totalorder>0.5), Nsnap
        write(451,*) Delta*count(totalorder>0.5), Delta*Nsnap
     close(451)
   
    end subroutine bondorder
    
!===============================================================================
!===============================================================================
! find the line beginning with pattern: '[CH] {4}\d' or '[CH] {3}\d{2}'
! '[CH]    [0-9]' or '[CH]   [0-9][0-9]'
! ([CH]    [0-9]|[CH]   [0-9][0-9])
! grep -E '([CH]    [0-9]|[CH]   [0-9][0-9])' collision.bond > temp.bond


!===============================================================================
!===============================================================================
! 
