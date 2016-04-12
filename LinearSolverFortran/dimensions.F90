module dimensions
  !
  implicit none
  ! Parameter
  Real*8                                 ::    pi
  parameter (pi = 3.141592653589793238462643383279)
  ! Parameters defining computational domain
  Real*8,dimension(:,:,:,:),allocatable  ::    coords
  Real*8                                 ::    Lx,Ly,Lz
  Real*8                                 ::    dx,dy,dz
  Integer                                ::    BCregions(6,3,2)
  Integer                                ::    BCType(6)
  Integer                                ::    nxb,nyb,nzb
  Integer                                ::    il_bnd,iu_bnd
  Integer                                ::    jl_bnd,ju_bnd
  Integer                                ::    kl_bnd,ku_bnd
  Integer                                ::    nguard,numVars
  
  ! Parameters defining test problem
  contains
    !
    subroutine computeCoordinates
      !
      implicit none
      !
      Integer  :: i,j,k
      !
      iu_bnd=nxb+2*nguard
      ju_bnd=nyb+2*nguard
#if (iins3d==1)
      ku_bnd=nzb+2*nguard
#else
      ku_bnd=1
#endif
      !
      allocate(coords(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:3))
      !
      dx=Lx/dble(nxb)
      dy=Ly/dble(nyb)
#if (iins3d==1)  
      dz=Lz/dble(nzb)
#endif
      !
      do k=kl_bnd,ku_bnd
         do j=jl_bnd,ju_bnd
            do i=il_bnd,iu_bnd
               !
               coords(i,j,k,1)=dx*(real(i-nguard)-0.5)
               coords(i,j,k,2)=dy*(real(j-nguard)-0.5)
#if (iins3d==1)
               coords(i,j,k,3)=dz*(real(k-nguard)-0.5)
#else
               coords(i,j,k,3)=0.d0
#endif         
               !
            end do
         end do
      end do
      !
      ! define domain ghost cell regions
      ! regular x-region
      BCregions(1:6,1,1)=nguard+1
      BCregions(1:6,1,2)=nguard+nxb
      ! regular y-region
      BCregions(1:6,2,1)=nguard+1
      BCregions(1:6,2,2)=nguard+nyb
      ! regular z-region
      BCregions(1:6,3,1)=nguard+1
      BCregions(1:6,3,2)=nguard+nzb
      ! left
      BCregions(1,1,1)=nguard
      BCregions(1,1,2)=1
      ! right
      BCregions(2,1,1)=nxb+nguard+1
      BCregions(2,1,2)=nxb+2*nguard
      ! bottom
      BCregions(3,2,1)=nguard
      BCregions(3,2,2)=1
      ! right
      BCregions(4,2,1)=nyb+nguard+1
      BCregions(4,2,2)=nyb+2*nguard
#if (iins3d==1)
      ! front
      BCregions(5,3,1)=nguard
      BCregions(5,3,2)=1
      ! back
      BCregions(6,3,1)=nzb+nguard+1
      BCregions(6,3,2)=nzb+2*nguard
#else
      ! in 2D
      BCregions(:,3,:)=1
#endif
      !
    end subroutine computeCoordinates
    !
    subroutine mkdir(dir)
!      use typedef
      !                                                                                                                                                                                                      
      implicit none
      !                                                                                                                                                                                                      
!      include 'mpif.h'
      !                                                                                                                                                                                                      
      character*120,intent(in)    ::   dir
      character*120               ::   cmd
      logical                     ::   dir_exists
      integer                     ::   ierr,errnum
      !                                                                                                                                                                                                      
 !     if (mypeno==0) then
         !                                                                                                                                                                                                   
      INQUIRE(FILE=trim(dir), EXIST=dir_exists)
      !                                                                                                                                                                                                   
      if (.not.dir_exists) then
         !                                                                                                                                                                                                
         write(*,'(3A)') "[I] Creating non-existing",trim(dir)," directory"
         write(cmd,'(2A)') "mkdir -p ",trim(dir)
         call SYSTEM(trim(cmd))
      else
         write(*,'(3A)') "[I] ",trim(dir)," directory exists"
      endif
      !      endif
      !                                                                                                                                                                                                      
      !  call MPI_Barrier(MPI_COMM_WORLD,ierr)
      !                                                                                                                                                                                                      
    end subroutine mkdir
    !
end module dimensions
