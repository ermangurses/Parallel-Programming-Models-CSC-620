#include "setup.h"
subroutine tecOutput(data,coords,nx,ny,nz,nv,filename)
  !
  implicit none
  !
  Real*8, intent(in)        ::  data  (1:nv,1:nx,1:ny,1:nz)
  Real*8, intent(in)        ::  coords(1:nx,1:ny,1:nz,1:2+iins3d)
  Integer, intent(in)       ::  nx,ny,nz,nv
  Character*120,intent(in)  ::  filename
  ! local
  Integer                   ::  funit
  Integer                   ::  i,j,k,ivar
  Character*200             ::  VariableString
  !
  funit=100
  ! open file and write header                                                                                                                                                                  
  open(unit=funit,file=trim(filename),status='unknown')
  !                                                                                                                                                                                             
  ! write file title                                                                                                                                                                            
  write(funit,'(A)')'TITLE = "Test Problem"'
  !                                                                                                                                                                                             
  ! assemble variable string          
#if (iins3d==0)                                                                                                                                                          
  VariableString = 'VARIABLES = "x","y"'
#else
  VariableString = 'VARIABLES = "x","y","z"'
#endif
  !                                                                                                                                                                                             
  do ivar = 1,nv
     !                                                                                                                                                                                          
!     write(VariableString,'(A,A,A,A)')trim(VariableString),',"',trim(VarNames(ivar)),'"'
     if (ivar.lt.10) then
        write(VariableString,'(A,A,I1,A)')trim(VariableString),',"var ',ivar,'"'
     else if (ivar.lt.100) then
        write(VariableString,'(A,A,I2,A)')trim(VariableString),',"var ',ivar,'"'
     else
        write(*,'(A)') "[E] linSolve: Too many variable names for output file"
     endif        
     !                                                                                                                                                                                          
  end do
  !                                                                                                                                                                                             
  ! write variable names to file                                                                                                                                                                
  write(funit,'(A)') trim(VariableString)
  !  
  write(funit,'(A,I6,A,I6,A,I6)')'ZONE I = ',nx,', J = ',ny,', K = ',nz
  !                                                                                                                                                                                             
  ! loop over coordinates of block                                                                                                                                                              
  do k=1,nz
     do j=1,ny
        do i=1,nx
           !
#if (iins3d==1)
           write(funit,'(3E27.16)',Advance='No') coords(i,j,k,1:3)
#else
           write(funit,'(2E27.16)',Advance='No') coords(i,j,k,1:2)
#endif
           !
           do ivar=1,nv-1
              write(funit,'(E27.16)',Advance='No') data(ivar,i,j,k)
           enddo
           write(funit,'(E27.16)') data(nv,i,j,k)
           !
        end do
     end do
  end do
  close(funit)
end subroutine tecOutput
