#include "setup.h"
module LinearSolver
  !
  use dimensions, only: coords,Lx,Ly,Lz,&
       dx,dy,dz,nguard,&
       nxb,nyb,nzb,&
       il_bnd,iu_bnd,&
       jl_bnd,ju_bnd,&
       kl_bnd,ku_bnd
  !
  implicit none
  !
  Logical                                      ::      errorOutput
  Integer                                      ::      maxIt,maxN,nVars,funitLin
  Real*8                                       ::      relaxFactor,errTol=1e-6
  Real*8,dimension(:,:,:,:,:,:),allocatable    ::      matLHS!(1:5+2*iins3d,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars,1:nVars)  
  Real*8,dimension(:,:,:,:),allocatable        ::      rhs!(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars)
  Real*8,dimension(:,:,:,:),allocatable        ::      sol!(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars)
  Real*8,dimension(:,:,:,:),allocatable        ::      error!(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars)
  !
  private 
  !
  public :: sol,assignRHS,assignLHS,solve,solveLinearSystem_ini,assignSol
  public :: error,errorOutput,evaluateRHS,rhs
  !
contains
  !
  subroutine solveLinearSystem_ini(maxItIn,errTolIn,relaxFactorIn,nVarsIn,nx,ny,nz)
    !
    use dimensions, only: il_bnd,iu_bnd,jl_bnd,ju_bnd,kl_bnd,ku_bnd
    !
    implicit none
    !
    Integer, intent(in)  ::  maxItIn
    Real*8,intent(in)    ::  errTolIn
    Real*8,intent(in)    ::  relaxFactorIn
    Integer, intent(in)  ::  nVarsIn
    Integer, intent(in)  ::  nx,ny,nz
    !
    !local
    Character*120        ::  filename
    !    
    maxIt       =maxItIn
    errTol      =errTolIn
    relaxFactor =relaxFactorIn
    nVars       =nVarsIn
    errorOutput =.true.
    maxN=max(nx,ny)    
#if (iins3d==1)
    maxN=max(MaxN,nz)
#endif
    !
    allocate(matLHS(1:5+2*iins3d,il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars,1:nVars))
    allocate(rhs(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars))
    allocate(sol(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars))
    if (errorOutput) allocate(error(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars))
    !
#if (OUTPUTRES==1)
    funitLin=101
    filename="post/stats.dat"
    open(unit=funitLin,file=trim(filename),status='unknown')
#endif    
    !
  end subroutine solveLinearSystem_ini
  !
  subroutine solve(scheme)
    !
    implicit none
    !
    Integer,intent(in)    ::    scheme
    !
    if (scheme.eq.1) then
       !
       call jacobiPointSolve
       !
    else if (scheme.eq.2) then
       !
       call jacobiLineSolve
       !
    else if (scheme.eq.3) then	
       !
       call jacobiPointSolveOMP       
    else
       !
       write(*,'(A,I2)') "[E] LinearSolve::solve: Linear solver not supported! scheme = ",scheme
       stop
       !
    endif
    !
  end subroutine solve
  !
  subroutine jacobiPointSolve
    !
    implicit none
    real*8             ::  solm1  (1:nxb+2*nguard,nyb+2*nguard,nzb+2*nguard*iins3d,1:nVars)
    real*8             ::  l2norm (1:nVars),linfnorm (1:nVars),dsol(1:nVars)
    real*8             ::  rhstmp2(1:nVars)
    real*8             ::  vec_tmp(1:nVars),sol_tmp(1:nVars)
    real*8             ::  Dinv(1:nVars,1:nVars)
    integer            ::  i,j,k,b,iter,v
    !
    do iter=1,maxIt
       !
       solm1(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)=sol(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)
       !
       l2norm(1:nVars)=0.d0
       !
       do k=kl_bnd+nguard*iins3d,ku_bnd-nguard*iins3d
          !
          do j=jl_bnd+nguard,ju_bnd-nguard
             !
             do i=il_bnd+nguard,iu_bnd-nguard
                ! invert diagonal block
                call inv(matLHS(1,i,j,k,1:nVars,1:nVars),Dinv(1:nVars,1:nVars),nVars)
                ! compute modified rhs
                rhstmp2(1:nVars)=rhs(i,j,k,1:nVars)
                !
                ! i-1,j,k
                vec_tmp(1:nVars)=matmul(matLHS(2,i,j,k,1:nVars,1:nVars),solm1(i-1,j,k,1:nVars))
                rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)
                !
                ! i+1,j,k
                vec_tmp(1:nVars)=matmul(matLHS(3,i,j,k,1:nVars,1:nVars),solm1(i+1,j,k,1:nVars))
                rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)
                !
                ! i,j-1,k
                vec_tmp(1:nVars)=matmul(matLHS(4,i,j,k,1:nVars,1:nVars),solm1(i,j-1,k,1:nVars))
                rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)
                !
                ! i,j+1,k
                vec_tmp(1:nVars)=matmul(matLHS(5,i,j,k,1:nVars,1:nVars),solm1(i,j+1,k,1:nVars))
                rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)
                !
#if (iins3d==1)
                ! i,j,k-1
                vec_tmp(1:nVars)=matmul(matLHS(6,i,j,k,1:nVars,1:nVars),solm1(i,j,k-1,1:nVars))
                rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)
                ! i,j,k+1
                vec_tmp(1:nVars)=matmul(matLHS(7,i,j,k,1:nVars,1:nVars),solm1(i,j,k+1,1:nVars))
                rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)                      
#endif
                !
                sol_tmp(1:nVars)=matmul(Dinv(1:nVars,1:nVars),rhstmp2(1:nVars))
                !
                dsol(1:nVars)=sol_tmp(1:nVars)-solm1(i,j,k,1:nVars)
                !
#if (OUTPUTRES==1)
                do v=1,nVars
                   l2norm  (v)=l2norm  (v)+dsol(v)*dsol(v)
!                   linfnorm(v)=max(linfnorm(v),abs(dsol(v)))
                enddo
#endif
                !
                sol (i,j,k,1:nVars)=solm1(i,j,k,1:nVars)+relaxFactor*dsol(1:nVars)
                !
             end do
          end do
       end do
       !
       call exchange_ghosts
       !
#if (OUTPUTRES==1)
       do v=1,nVars
          l2norm  (v)=sqrt(l2norm(v)/nVars)
       enddo
       !
       write(funitLin,*) iter,l2norm(1:nVars)!,linfnorm(1:nVars)
#endif
       !
    enddo
    !
  end subroutine jacobiPointSolve
   
subroutine jacobiPointSolveOMP
    !
    implicit none
    real*8             ::  solm1  (1:nxb+2*nguard,nyb+2*nguard,nzb+2*nguard*iins3d,1:nVars)
    real*8             ::  l2norm (1:nVars),linfnorm (1:nVars),dsol(1:nVars)
    real*8             ::  rhstmp2(1:nVars)
    real*8             ::  vec_tmp1(1:nVars)
    real*8             ::  vec_tmp2(1:nVars) 
    real*8             ::  vec_tmp3(1:nVars)
    real*8             ::  vec_tmp4(1:nVars) 
    real*8             ::  vec_tmp5(1:nVars)
    real*8             ::  vec_tmp6(1:nVars)       
    real*8             ::  sol_tmp(1:nVars)
    real*8             ::  Dinv(1:nVars,1:nVars)
    integer            ::  i,j,k,b,iter,v
    !
    do iter=1,maxIt
       !
       solm1(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)=sol(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)
       !
       l2norm(1:nVars)=0.d0
       !
       !$omp  parallel do        
       !$omp& shared (nVars, rhstmp2, solm1, matLHS)
       !$omp& private (i,j,k,vec_tmp1, vec_tmp2, vec_tmp3, vec_tmp4, vec_tmp5, vec_tmp6) 
       do k=kl_bnd+nguard*iins3d,ku_bnd-nguard*iins3d
          !
          do j=jl_bnd+nguard,ju_bnd-nguard
             !
             do i=il_bnd+nguard,iu_bnd-nguard
                ! invert diagonal block
                call inv(matLHS(1,i,j,k,1:nVars,1:nVars),Dinv(1:nVars,1:nVars),nVars)
                ! compute modified rhs
                rhstmp2(1:nVars)=rhs(i,j,k,1:nVars)
                !
                ! i-1,j,k
                vec_tmp1(1:nVars)=matmul(matLHS(2,i,j,k,1:nVars,1:nVars),solm1(i-1,j,k,1:nVars))
                !rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)
                !
                ! i+1,j,k
                vec_tmp2(1:nVars)=matmul(matLHS(3,i,j,k,1:nVars,1:nVars),solm1(i+1,j,k,1:nVars))
                !rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)
                !
                ! i,j-1,k
                vec_tmp3(1:nVars)=matmul(matLHS(4,i,j,k,1:nVars,1:nVars),solm1(i,j-1,k,1:nVars))
                !rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)
                !
                ! i,j+1,k
                vec_tmp4(1:nVars)=matmul(matLHS(5,i,j,k,1:nVars,1:nVars),solm1(i,j+1,k,1:nVars))
                !rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)
                !
#if (iins3d==1)
                ! i,j,k-1
                vec_tmp5(1:nVars)=matmul(matLHS(6,i,j,k,1:nVars,1:nVars),solm1(i,j,k-1,1:nVars))
                !rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)
                ! i,j,k+1
                vec_tmp6(1:nVars)=matmul(matLHS(7,i,j,k,1:nVars,1:nVars),solm1(i,j,k+1,1:nVars))
                !rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp(1:nVars)                      
#endif
                rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp2(1:nVars)&
                                                 -vec_tmp3(1:nVars)&
                                                 -vec_tmp4(1:nVars)&
                                                 -vec_tmp5(1:nVars)
#if (iins3d==1)                                                 
                rhstmp2(1:nVars)=rhstmp2(1:nVars)-vec_tmp6(1:nVars)&
                                                 -vec_tmp7(1:nVars)
#endif                
                !
                sol_tmp(1:nVars)=matmul(Dinv(1:nVars,1:nVars),rhstmp2(1:nVars))
                !
                dsol(1:nVars)=sol_tmp(1:nVars)-solm1(i,j,k,1:nVars)
                !
#if (OUTPUTRES==1)
                do v=1,nVars
                   l2norm  (v)=l2norm  (v)+dsol(v)*dsol(v)
!                   linfnorm(v)=max(linfnorm(v),abs(dsol(v)))
                enddo
#endif
                !
                sol (i,j,k,1:nVars)=solm1(i,j,k,1:nVars)+relaxFactor*dsol(1:nVars)
                !
             end do
          end do
       end do
       !$omp end parallel do  
       !
       call exchange_ghosts
       !
#if (OUTPUTRES==1)
       do v=1,nVars
          l2norm  (v)=sqrt(l2norm(v)/nVars)
       enddo
       !
       write(funitLin,*) iter,l2norm(1:nVars)!,linfnorm(1:nVars)
#endif
       !
    enddo
    !
  end subroutine jacobiPointSolveOMP
  !
  subroutine jacobiLineSolve
    !
    implicit none
    Real*8             ::  solm1  (1:nxb+2*nguard,nyb+2*nguard,nzb+2*nguard*iins3d,1:nVars)
    Real*8             ::  l2norm (1:nVars),linfnorm (1:nVars),dsol(1:nVars)
    Real*8             ::  rhstmp2(1:nVars,1:maxN)
    Real*8             ::  vec_tmp(1:nVars)
    Real*8             ::  sol_tmp(1:nVars,1:maxN)
    Real*8             ::  CoefMat_loc(1:3,1:nVars,1:nVars,1:maxN)
    integer            ::  i,j,k,b,iter,v,count,numPts
    !
    do iter=1,maxIt
       !
       solm1(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)=sol(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)
       !
       l2norm(1:nVars)=0.d0
       !
       ! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       !
       !       do lb=1,lnblocks ! loop over all local blocks 
       !
       ! x-solve
       !
       do k=kl_bnd+nguard*iins3d,ku_bnd-nguard*iins3d
          !
          do j=jl_bnd+nguard,ju_bnd-nguard
             !
             !add first and last coefficients
             !
             count=1 
             i=il_bnd+nguard
             !
             CoefMat_loc(1,1:nVars,1:nVars,count)=matLHS(2,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(2,1:nVars,1:nVars,count)=matLHS(1,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(3,1:nVars,1:nVars,count)=matLHS(3,i,j,k,1:nVars,1:nVars)
             !
             rhstmp2(1:nVars,count)=rhs(i,j,k,1:nVars)
             ! i-1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(2,i,j,k,1:nVars,1:nVars),solm1(i-1,j,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i,j-1,k
             vec_tmp(1:nVars)=matmul(matLHS(4,i,j,k,1:nVars,1:nVars),solm1(i,j-1,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i,j+1,k
             vec_tmp(1:nVars)=matmul(matLHS(5,i,j,k,1:nVars,1:nVars),solm1(i,j+1,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
             !
#if (iins3d==1)
             ! i,j,k-1
             vec_tmp(1:nVars)=matmul(matLHS(6,i,j,k,1:nVars,1:nVars),solm1(i,j,k-1,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i,j,k+1
             vec_tmp(1:nVars)=matmul(matLHS(7,i,j,k,1:nVars,1:nVars),solm1(i,j,k+1,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
#endif
             !
             count=2
             !
             do i=il_bnd+nguard+1,iu_bnd-nguard-1
                ! 
                CoefMat_loc(1,1:nVars,1:nVars,count)=matLHS(2,i,j,k,1:nVars,1:nVars)
                CoefMat_loc(2,1:nVars,1:nVars,count)=matLHS(1,i,j,k,1:nVars,1:nVars)
                CoefMat_loc(3,1:nVars,1:nVars,count)=matLHS(3,i,j,k,1:nVars,1:nVars)
                !
                rhstmp2(1:nVars,count)=rhs(i,j,k,1:nVars)
                vec_tmp(1:nVars)=matmul(matLHS(4,i,j,k,1:nVars,1:nVars),solm1(i,j-1,k,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
                !
                vec_tmp(1:nVars)=matmul(matLHS(5,i,j,k,1:nVars,1:nVars),solm1(i,j+1,k,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                      
                !
#if (iins3d==1)
                ! i,j,k-1
                vec_tmp(1:nVars)=matmul(matLHS(6,i,j,k,1:nVars,1:nVars),solm1(i,j,k-1,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
                ! i,j,k+1
                vec_tmp(1:nVars)=matmul(matLHS(7,i,j,k,1:nVars,1:nVars),solm1(i,j,k+1,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
#endif
                !
                count=count+1
                !
             enddo
             !
             i=iu_bnd-nguard
             !
             CoefMat_loc(1,1:nVars,1:nVars,count)=matLHS(2,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(2,1:nVars,1:nVars,count)=matLHS(1,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(3,1:nVars,1:nVars,count)=matLHS(3,i,j,k,1:nVars,1:nVars)
             !
             rhstmp2(1:nVars,count)=rhs(i,j,k,1:nVars)
             ! i+1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(3,i,j,k,1:nVars,1:nVars),solm1(i+1,j,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i,j-1,k
             vec_tmp(1:nVars)=matmul(matLHS(4,i,j,k,1:nVars,1:nVars),solm1(i,j-1,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i,j+1,k
             vec_tmp(1:nVars)=matmul(matLHS(5,i,j,k,1:nVars,1:nVars),solm1(i,j+1,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                                                         
#if (iins3d==1)
             ! i,j,k-1
             vec_tmp(1:nVars)=matmul(matLHS(6,i,j,k,1:nVars,1:nVars),solm1(i,j,k-1,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i,j,k+1
             vec_tmp(1:nVars)=matmul(matLHS(7,i,j,k,1:nVars,1:nVars),solm1(i,j,k+1,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
#endif
             !
             numPts=count
             !
             call blk_tridiag_solve(CoefMat_loc(1,1:nVars,1:nVars,1:numPts),&
                  CoefMat_loc(2,1:nVars,1:nVars,1:numPts),&
                  CoefMat_loc(3,1:nVars,1:nVars,1:numPts),&
                  sol_tmp(1:nVars,1:numPts),&
                  rhstmp2(1:nVars,1:numPts),&
                  numPts,nVars)
             !
             count=1
             do i=il_bnd+nguard,iu_bnd-nguard
                !
                dsol(1:nVars)=sol_tmp(1:nVars,count)-solm1(i,j,k,1:nVars)
                !
                sol (i,j,k,1:nVars)=solm1(i,j,k,1:nVars)+relaxFactor*dsol(1:nVars)
                !
#if (OUTPUTRES==1)
                do v=1,nVars
                   l2norm  (v)=l2norm  (v)+dsol(v)*dsol(v)
!                   linfnorm(v)=max(linfnorm(v),abs(dsol(v)))
                enddo
#endif
                !
                count=count+1
                !
             enddo
             !                   
          end do
       end do
       !       end do
       !
#if (OUTPUTRES==1)
       do v=1,nVars
          l2norm  (v)=sqrt(l2norm(v)/nVars)
       enddo
       !
       write(funitLin,*) iter,l2norm(1:nVars)!,linfnorm(1:nVars)
#endif
       !
       call exchange_ghosts
       !
       solm1(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)=sol(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)
       !
#if (OUTPUTRES==1)
       l2norm(1:nVars)=0.d0
#endif
       !
       ! YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
       !
       !       do lb=1,lnblocks ! loop over all local blocks 
       !
       ! y-solve
       !
       do k=kl_bnd+nguard*iins3d,ku_bnd-nguard*iins3d
          !
          do i=il_bnd+nguard,iu_bnd-nguard
             !
             !add first and last coefficients
             !
             count=1 
             j=jl_bnd+nguard
             !
             CoefMat_loc(1,1:nVars,1:nVars,count)=matLHS(4,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(2,1:nVars,1:nVars,count)=matLHS(1,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(3,1:nVars,1:nVars,count)=matLHS(5,i,j,k,1:nVars,1:nVars)
             !
             rhstmp2(1:nVars,count)=rhs(i,j,k,1:nVars)
             ! i,j-1,k
             vec_tmp(1:nVars)=matmul(matLHS(4,i,j,k,1:nVars,1:nVars),solm1(i,j-1,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i-1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(2,i,j,k,1:nVars,1:nVars),solm1(i-1,j,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i+1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(3,i,j,k,1:nVars,1:nVars),solm1(i+1,j,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
             !
#if (iins3d==1)
             ! i,j,k-1
             vec_tmp(1:nVars)=matmul(matLHS(6,i,j,k,1:nVars,1:nVars),solm1(i,j,k-1,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i,j,k+1
             vec_tmp(1:nVars)=matmul(matLHS(7,i,j,k,1:nVars,1:nVars),solm1(i,j,k+1,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
#endif
             !
             count=2
             !
             do j=jl_bnd+nguard+1,ju_bnd-nguard-1
                ! 
                CoefMat_loc(1,1:nVars,1:nVars,count)=matLHS(4,i,j,k,1:nVars,1:nVars)
                CoefMat_loc(2,1:nVars,1:nVars,count)=matLHS(1,i,j,k,1:nVars,1:nVars)
                CoefMat_loc(3,1:nVars,1:nVars,count)=matLHS(5,i,j,k,1:nVars,1:nVars)
                !
                rhstmp2(1:nVars,count)=rhs(i,j,k,1:nVars)
                ! i-1,j,k
                vec_tmp(1:nVars)=matmul(matLHS(2,i,j,k,1:nVars,1:nVars),solm1(i-1,j,k,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
                ! i+1,j,k
                vec_tmp(1:nVars)=matmul(matLHS(3,i,j,k,1:nVars,1:nVars),solm1(i+1,j,k,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
#if (iins3d==1)
                ! i,j,k-1
                vec_tmp(1:nVars)=matmul(matLHS(6,i,j,k,1:nVars,1:nVars),solm1(i,j,k-1,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
                ! i,j,k+1
                vec_tmp(1:nVars)=matmul(matLHS(7,i,j,k,1:nVars,1:nVars),solm1(i,j,k+1,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
#endif
                !
                count=count+1
                !
             enddo
             !
             j=ju_bnd-nguard
             !
             CoefMat_loc(1,1:nVars,1:nVars,count)=matLHS(4,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(2,1:nVars,1:nVars,count)=matLHS(1,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(3,1:nVars,1:nVars,count)=matLHS(5,i,j,k,1:nVars,1:nVars)
             !
             rhstmp2(1:nVars,count)=rhs(i,j,k,1:nVars)
             ! i,j+1,k
             vec_tmp(1:nVars)=matmul(matLHS(5,i,j,k,1:nVars,1:nVars),solm1(i,j+1,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i-1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(2,i,j,k,1:nVars,1:nVars),solm1(i-1,j,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i+1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(3,i,j,k,1:nVars,1:nVars),solm1(i+1,j,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                                                         
             !
#if (iins3d==1)
             ! i,j,k-1
             vec_tmp(1:nVars)=matmul(matLHS(6,i,j,k,1:nVars,1:nVars),solm1(i,j,k-1,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i,j,k+1
             vec_tmp(1:nVars)=matmul(matLHS(7,i,j,k,1:nVars,1:nVars),solm1(i,j,k+1,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
#endif
             !
             numPts=count
             !
             call blk_tridiag_solve(CoefMat_loc(1,1:nVars,1:nVars,1:numPts),&
                  CoefMat_loc(2,1:nVars,1:nVars,1:numPts),&
                  CoefMat_loc(3,1:nVars,1:nVars,1:numPts),&
                  sol_tmp(1:nVars,1:numPts),&
                  rhstmp2(1:nVars,1:numPts),&
                  numPts,nVars)
             !
             count=1
             do j=jl_bnd+nguard,ju_bnd-nguard
                !
                dsol(1:nVars)=sol_tmp(1:nVars,count)-solm1(i,j,k,1:nVars)
                !
#if (OUTPUTRES==1)
                do v=1,nVars
                   l2norm  (v)=l2norm  (v)+dsol(v)*dsol(v)
!                   linfnorm(v)=max(linfnorm(v),abs(dsol(v)))
                enddo
#endif
                !
                sol (i,j,k,1:nVars)=solm1(i,j,k,1:nVars)+relaxFactor*dsol(1:nVars)
                !
                count=count+1
                !
             enddo
             !                   
          end do
       end do
       !
       !       end do
       !     
      call exchange_ghosts
      !
       solm1(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)=sol(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)
       !
#if (iins3d==1)
       !
       ! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
       !
       !       do lb=1,lnblocks ! loop over all local blocks 
       !

       !
       ! z-solve
       !
       do i=il_bnd+nguard,iu_bnd-nguard
          do j=jl_bnd+nguard,ju_bnd-nguard
             !
             !add first and last coefficients
             !
             count=1 
             k=kl_bnd+nguard
             !
             CoefMat_loc(1,1:nVars,1:nVars,count)=matLHS(6,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(2,1:nVars,1:nVars,count)=matLHS(1,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(3,1:nVars,1:nVars,count)=matLHS(7,i,j,k,1:nVars,1:nVars)
             !
             rhstmp2(1:nVars,count)=rhs(i,j,k,1:nVars)
             ! i,j,k-1
             vec_tmp(1:nVars)=matmul(matLHS(6,i,j,k,1:nVars,1:nVars),solm1(i,j,k-1,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i-1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(2,i,j,k,1:nVars,1:nVars),solm1(i-1,j,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i+1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(3,i,j,k,1:nVars,1:nVars),solm1(i+1,j,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
             ! i,j-1,k
             vec_tmp(1:nVars)=matmul(matLHS(4,i,j,k,1:nVars,1:nVars),solm1(i,j-1,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i,j+1,k
             vec_tmp(1:nVars)=matmul(matLHS(5,i,j,k,1:nVars,1:nVars),solm1(i,j+1,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
             !
             count=2
             !
             do k=kl_bnd+nguard*iins3d+1,ku_bnd-nguard*iins3d-1
                ! 
                CoefMat_loc(1,1:nVars,1:nVars,count)=matLHS(6,i,j,k,1:nVars,1:nVars)
                CoefMat_loc(2,1:nVars,1:nVars,count)=matLHS(1,i,j,k,1:nVars,1:nVars)
                CoefMat_loc(3,1:nVars,1:nVars,count)=matLHS(7,i,j,k,1:nVars,1:nVars)
                !
                rhstmp2(1:nVars,count)=rhs(i,j,k,1:nVars)
                ! i-1,j,k
                vec_tmp(1:nVars)=matmul(matLHS(2,i,j,k,1:nVars,1:nVars),solm1(i-1,j,k,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
                ! i+1,j,k
                vec_tmp(1:nVars)=matmul(matLHS(3,i,j,k,1:nVars,1:nVars),solm1(i+1,j,k,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
                ! i,j-1,k
                vec_tmp(1:nVars)=matmul(matLHS(4,i,j,k,1:nVars,1:nVars),solm1(i,j-1,k,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
                ! i,j+1,k
                vec_tmp(1:nVars)=matmul(matLHS(5,i,j,k,1:nVars,1:nVars),solm1(i,j+1,k,1:nVars))
                rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
                !
                count=count+1
                !
             enddo
             !
             k=ku_bnd-nguard
             !
             CoefMat_loc(1,1:nVars,1:nVars,count)=matLHS(6,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(2,1:nVars,1:nVars,count)=matLHS(1,i,j,k,1:nVars,1:nVars)
             CoefMat_loc(3,1:nVars,1:nVars,count)=matLHS(7,i,j,k,1:nVars,1:nVars)
             !
             rhstmp2(1:nVars,count)=rhs(i,j,k,1:nVars)
             ! i,j,k+1
             vec_tmp(1:nVars)=matmul(matLHS(7,i,j,k,1:nVars,1:nVars),solm1(i,j,k+1,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i-1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(2,i,j,k,1:nVars,1:nVars),solm1(i-1,j,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i+1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(3,i,j,k,1:nVars,1:nVars),solm1(i+1,j,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                                                         
             ! i,j-1,k
             vec_tmp(1:nVars)=matmul(matLHS(4,i,j,k,1:nVars,1:nVars),solm1(i,j-1,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)
             ! i,j+1,k
             vec_tmp(1:nVars)=matmul(matLHS(5,i,j,k,1:nVars,1:nVars),solm1(i,j+1,k,1:nVars))
             rhstmp2(1:nVars,count)=rhstmp2(1:nVars,count)-vec_tmp(1:nVars)                   
             !
             numPts=count
             !
             call blk_tridiag_solve(CoefMat_loc(1,1:nVars,1:nVars,1:numPts),&
                  CoefMat_loc(2,1:nVars,1:nVars,1:numPts),&
                  CoefMat_loc(3,1:nVars,1:nVars,1:numPts),&
                  sol_tmp(1:nVars,1:numPts),&
                  rhstmp2(1:nVars,1:numPts),&
                  numPts,nVars)
             !
             count=1
             do k=kl_bnd+nguard,ku_bnd-nguard
                !
                dsol(1:nVars)=sol_tmp(1:nVars,count)-solm1(i,j,k,1:nVars)
                !
#if (OUTPUTRES==1)
                do v=1,nVars
                   l2norm  (v)=l2norm  (v)+dsol(v)*dsol(v)
!                   linfnorm(v)=max(linfnorm(v),abs(dsol(v)))
                enddo
#endif
                !
                sol (i,j,k,1:nVars)=solm1(i,j,k,1:nVars)+relaxFactor*dsol(1:nVars)
                !
                count=count+1
                !
             enddo
             !                   
          end do
       end do
       !
       !       end do
       !
       call exchange_ghosts
       solm1(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)=sol(1:nxb+2*nguard,1:nyb+2*nguard,1:nzb+2*nguard*iins3d,1:nVars)
       !
#endif       
       !
#if (OUTPUTRES==1)
       do v=1,nVars
          l2norm  (v)=sqrt(l2norm(v)/nVars)
       enddo
       !
       write(funitLin,*) iter,l2norm(1:nVars)!,linfnorm(1:nVars)
       flush(funitLin)
#endif
       !
    enddo
#if (OUTPUTRES==1)
    close(funitLin)
#endif    
    !
  end subroutine jacobiLineSolve
  !
  subroutine blk_tridiag_solve(bm,bd,bp,x,b,jmax,neqs)
    !*********************************************************************
    implicit none
    integer  :: jmax,j,neqs
    Real(8)  :: bm(neqs,neqs,jmax)
    Real(8)  :: bd(neqs,neqs,jmax)
    Real(8)  :: bp(neqs,neqs,jmax)
    Real(8)  :: x(neqs,jmax)
    Real(8)  :: b(neqs,jmax)
    Real(8)  :: T(neqs,neqs),Dinv(neqs,neqs),y(neqs)
    Real(8)  :: SD(neqs,neqs),sb(neqs),sx(neqs)
    !-----
    !     Solve blk_tridiagonal system
    !-----
    do j = 2,jmax
       call inv(bd(:,:,j-1),Dinv,neqs)
       call matmat(bm(:,:,j),Dinv,T,neqs)
       call matmat(T,bp(:,:,j-1),SD,neqs)
       call matvec(T,b(:,j-1),sb,neqs)
       bd(:,:,j) = bd(:,:,j) - SD
       b(:,j)    = b(:,j)    - sb
    enddo
    call lu_solve(bd(:,:,jmax),x(:,jmax),b(:,jmax),neqs)
    do j = jmax-1,1,-1
       call matvec(bp(:,:,j),x(:,j+1),sx,neqs)
       y = b(:,j) - sx
       call lu_solve(bd(:,:,j),x(:,j),y,neqs)
    enddo
    !-----
    !     End blk_tridiag_solve
    !-----
    return
  end subroutine blk_tridiag_solve
  !
  subroutine inv(A,B,neqs)
    !*********************************************************************
    implicit none
    integer :: i,k,neqs
    Real(8)  :: A(neqs,neqs),B(neqs,neqs),C(neqs,neqs)
    !-----
    !     Initialize B = I and C = A
    !-----
    do k = 1,neqs
       do i = 1,neqs
          B(i,k) = 0.0
          C(i,k) = A(i,k)
       enddo
       B(k,k) = 1.0
    enddo
    !-----
    !     Row reduce
    !-----
    do k = 1,neqs-1
       B(k,:) = B(k,:)/C(k,k)
       C(k,:) = C(k,:)/C(k,k)
       do i = k+1,neqs
          B(i,:) = B(i,:) - C(i,k)*B(k,:)
          C(i,:) = C(i,:) - C(i,k)*C(k,:)
       enddo
    enddo
    B(neqs,:) = B(neqs,:)/C(neqs,neqs)
    C(neqs,:) = C(neqs,:)/C(neqs,neqs)
    do k = neqs,2,-1
       do i = 1,k-1
          B(i,:) = B(i,:) - C(i,k)*B(k,:)
          C(i,:) = C(i,:) - C(i,k)*C(k,:)
       enddo
    enddo
    !-----
    !     End inv
    !-----
    return
  end subroutine inv
  !
  !*********************************************************************
  subroutine lu_solve(A,x,b,neqs)
    !*********************************************************************

    implicit none
    integer :: neqs
    Real(8) :: A(neqs,neqs),L(neqs,neqs),U(neqs,neqs)
    Real(8) :: x(neqs),y(neqs),b(neqs)
    !-----
    !     Factor L*U = A
    !-----
    call lu(A,L,U,neqs)
    !-----
    !     Solve lower triangular system
    !-----
    call lower_tri_solve(L,y,b,neqs)
    !-----
    !     Solve upper triangular system
    !-----
    call upper_tri_solve(U,x,y,neqs)
    !-----
    !     End lu_solve
    !-----
    return
  end subroutine lu_solve
  !
  !*********************************************************************
  subroutine lu(A,L,U,neqs)
    !*********************************************************************

    implicit none
    integer :: i,j,k,neqs
    Real(8)  :: A(neqs,neqs),T(neqs,neqs)
    Real(8)  :: L(neqs,neqs),U(neqs,neqs)
    Real(8)  :: mult
    !-----
    !     Set T = A
    !-----
    do j = 1,neqs
       do i = 1,neqs
          T(i,j) = A(i,j)
       enddo
    enddo
    !-----
    !     Factor A
    !-----
    do k = 1,neqs-1
       mult = 1.0/T(k,k)
       do i = k+1,neqs
          T(i,k) = mult*T(i,k)
          do j = k+1,neqs
             T(i,j) = T(i,j) - T(i,k)*T(k,j)
          enddo
       enddo
    enddo
    !-----
    !     Store factors in L and U
    !-----
    do i = 1,neqs
       do j = 1,neqs
          if(i.eq.j)then
             L(i,j) = 1.0
             U(i,j) = T(i,j)
          elseif(i.lt.j)then
             L(i,j) = 0.0
             U(i,j) = T(i,j)
          else
             L(i,j) = T(i,j)
             U(i,j) = 0.0
          endif
       enddo
    enddo
    !-----
    !     End lu
    !-----
    return
  end subroutine lu
  !
  !*********************************************************************
  subroutine upper_tri_solve(U,x,b,neqs)
    !*********************************************************************

    implicit none
    integer :: i,j,neqs
    Real(8)  :: U(neqs,neqs),x(neqs),b(neqs)
    Real(8)  :: sum
    !-----
    !     Back solve
    !-----
    x(neqs) = b(neqs)/U(neqs,neqs)
    do i = neqs-1,1,-1
       sum = 0.0
       do j = i+1,neqs
          sum = sum + U(i,j)*x(j)
       enddo
       x(i) = (b(i) - sum)/U(i,i)
    enddo
    !-----
    !     End upper_tri_solve
    !-----
    return
  end subroutine upper_tri_solve
  !
  !*********************************************************************
  subroutine lower_tri_solve(L,x,b,neqs)
    !*********************************************************************

    implicit none
    integer :: i,j,neqs
    Real(8)  :: L(neqs,neqs),x(neqs),b(neqs)
    Real(8)  :: sum
    !-----
    !     Foward solve
    !-----
    x(1) = b(1)/L(1,1)
    do i = 2,neqs
       sum = 0.0
       do j = 1,i-1
          sum = sum + L(i,j)*x(j)
       enddo
       x(i) = (b(i) - sum)/L(i,i)
    enddo
    !-----
    !     End lower_tri_solve
    !-----
    return
  end subroutine lower_tri_solve
  !
  !*********************************************************************
  subroutine matvec(A,b,c,n)
    !*********************************************************************
    implicit none
    integer :: n
    Real(8)  :: A(n,n),B(n),C(n)
    Real(8)  :: s
    integer :: i,k
    !-----
    !     Compute c = A*b
    !-----
    do i = 1,n
       s = 0.0
       do k = 1,n
          s = s + A(i,k)*b(k)
       enddo
       c(i) = s
    enddo
    !-----
    !     End matvec
    !-----
    return
  end subroutine matvec
  !
  !*********************************************************************
  subroutine matmat(A,B,C,n)
    !*********************************************************************
    implicit none
    integer :: n
    Real(8)  :: A(n,n),B(n,n),C(n,n)
    Real(8)  :: s
    integer :: i,j,k
    !-----
    !     Compute C = A*B
    !-----
    do j = 1,n
       do i = 1,n
          s = 0.0
          do k = 1,n
             s = s + A(i,k)*B(k,j)
          enddo
          C(i,j) = s
       enddo
    enddo
    !-----
    !     End matmat
    !-----
    return
  end subroutine matmat
  !*********************************************************************
  !
  subroutine assignLHS
    !
    use dimensions, only: il_bnd,iu_bnd,jl_bnd,ju_bnd,kl_bnd,ku_bnd
    !
    implicit none
    !
    Integer      ::  i,j,k,d
    Real*8       ::  rdx2,rdy2,rdz2
    !
    !    do lb=1,lnblocks ! loop over all local blocks 
    !
    do k=kl_bnd+nguard*iins3d,ku_bnd-nguard*iins3d
       !
       do j=jl_bnd+nguard,ju_bnd-nguard
          !
          do i=il_bnd+nguard,iu_bnd-nguard
             !
             rdx2=1.d0/(dx*dx)
             rdy2=1.d0/(dy*dy)
             rdz2=1.d0/(dz*dz)
             !
             matLHS(1:5+iins3d*2,i,j,k,1:nVars,1:nVars)=0.
             do d=1,nVars
                !                      
                matLHS(1  ,i,j,k,d,d)=-2.d0*rdx2-2.d0*rdy2-dble(iins3d)*2.d0*rdz2
                matLHS(2:3,i,j,k,d,d)=1.d0*rdx2
                matLHS(4:5,i,j,k,d,d)=1.d0*rdy2
#if (iins3d==1)
                matLHS(6:7,i,j,k,d,d)=1.d0*rdz2
#endif
                !
             end do
             !
          end do
       end do
    end do
    ! end do
    !
  end subroutine assignLHS
  !
  subroutine assignSol
    !
    use dimensions, only: il_bnd,iu_bnd,jl_bnd,ju_bnd,kl_bnd,ku_bnd,coords
    !
    implicit none
    !
    Integer      ::  i,j,k,d
    Real*8       ::  xyz(3)
    !
    !  do lb=1,lnblocks ! loop over all local blocks 
    !
    !
    do k=kl_bnd,ku_bnd
       !
       do j=jl_bnd,ju_bnd
          !
          do i=il_bnd,iu_bnd
             !
             xyz(1:3)=coords(i,j,k,1:3)
             !
             !call userSpecifiedFunction(sol(i,j,k,1:nVars),nVars,xyz(1:3),Lx,Ly,Lz)
             sol(i,j,k,1:nVars)=0.
          end do
       end do
    end do
    call exchange_ghosts
    !  end do
    !
  end subroutine assignSol
  !
  subroutine userSpecifiedFunction(data,nv,xyz,Lx,Ly,Lz)
    !
    use dimensions, only: pi
    !
    implicit none
    !
    Real*8,intent(out)        ::    data(nv)
    Integer,intent(in)        ::    nv
    Real*8,intent(in)         ::    xyz(3)
    Real*8,intent(in)         ::    Lx,Ly,Lz
    ! local
    Real*8                    ::    lambdax,lambday,lambdaz
    !
    lambdax=Lx
    lambday=Ly
    lambdaz=Lz
    !
#if (iins3d==1)
!    data(1:nv)=sin(2.*pi*xyz(1)/lambdax)*sin(2.*pi*xyz(2)/lambday)*sin(2.*pi*xyz(3)/lambdaz)
    data(1:nv)=cos(2.*pi*xyz(1)/lambdax)*cos(2.*pi*xyz(2)/lambday)*cos(2.*pi*xyz(3)/lambdaz)
#else
!    data(1:nv)=sin(2.*pi*xyz(1)/lambdax)*sin(2.*pi*xyz(2)/lambday)
    data(1:nv)=cos(2.*pi*xyz(1)/lambdax)*cos(2.*pi*xyz(2)/lambday)
#endif  
    !
  end subroutine userSpecifiedFunction
  !
  subroutine assignRHS
    !
    use dimensions, only: pi
    !
    implicit none
    !
    Integer      ::  i,j,k,idir,m
    Real*8       ::  xyz(3),rdx2,rdy2,rdz2
    Real*8       ::  smallNumber,phixx(nVars),data(10,nVars),xyz_eps(3)
    !
    smallNumber=1e-4
    !
    !  do lb=1,lnblocks ! loop over all local blocks 
    !
    do k=kl_bnd+nguard*iins3d,ku_bnd-nguard*iins3d
       !
       do j=jl_bnd+nguard,ju_bnd-nguard
          !
          do i=il_bnd+nguard,iu_bnd-nguard
             !
             xyz(1:3)=coords(i,j,k,1:3)
             !                   
             xyz_eps=xyz
! 
! rhs(i,j,k,1:nVars)=-8.d0*pi*pi*sin(2.d0*pi*xyz(1))*sin(2.d0*pi*xyz(2))
! using numerical differentation instead of analytical solution to it more general
!             
             rhs(i,j,k,1:nVars)=0.
             !
             do idir=1,2+iins3d
                !
                do m=1,7
                   xyz_eps(idir)=xyz(idir)+dble(m-4)*smallNumber                   
                   call userSpecifiedFunction(data(m,1:nVars),nVars,xyz_eps(1:3),Lx,Ly,Lz)
                end do
                !
                phixx(1:nvars)=( 1.d0 /90.d0*(data(1,1:nVars)+data(7,1:nVars))&
                     -3.d0 /20.d0*(data(2,1:nVars)+data(6,1:nVars))&
                     +3.d0 / 2.d0*(data(3,1:nVars)+data(5,1:nVars))&
                     -49.d0/18.d0*(data(4,1:nVars)))/(smallNumber*smallNumber)
                !
                rhs(i,j,k,1:nVars)=rhs(i,j,k,1:nVars)+phixx(1:nVars)
                !
             end do
             !
          end do
       end do
    end do
    !       end do
    !
  end subroutine assignRHS
  !
  subroutine evaluateRHS(data,rhstmp)
    !
    implicit none
    !
    Real*8, intent(in)     ::   data  (il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars)
    Real*8, intent(out)    ::   rhstmp(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars)
    ! local
    Integer                ::   i,j,k
    Real*8                 ::   vec_tmp(1:nVars)
    !
    do k=kl_bnd+nguard*iins3d,ku_bnd-nguard*iins3d
       !
       do j=jl_bnd+nguard,ju_bnd-nguard
          !
          do i=il_bnd+nguard,iu_bnd-nguard
             !
             rhstmp(i,j,k,1:nVars)=0.d0
             !
             ! i,j,k
             vec_tmp(1:nVars)=matmul(matLHS(1,i,j,k,1:nVars,1:nVars),data(i,j,k,1:nVars))
             rhstmp(i,j,k,1:nVars)=rhstmp(i,j,k,1:nVars)-vec_tmp(1:nVars)
             ! i-1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(2,i,j,k,1:nVars,1:nVars),data(i-1,j,k,1:nVars))
             rhstmp(i,j,k,1:nVars)=rhstmp(i,j,k,1:nVars)-vec_tmp(1:nVars)
             ! i+1,j,k
             vec_tmp(1:nVars)=matmul(matLHS(3,i,j,k,1:nVars,1:nVars),data(i+1,j,k,1:nVars))
             rhstmp(i,j,k,1:nVars)=rhstmp(i,j,k,1:nVars)-vec_tmp(1:nVars)
             ! i,j-1,k
             vec_tmp(1:nVars)=matmul(matLHS(4,i,j,k,1:nVars,1:nVars),data(i,j-1,k,1:nVars))
             rhstmp(i,j,k,1:nVars)=rhstmp(i,j,k,1:nVars)-vec_tmp(1:nVars)
             ! i,j+1,k
             vec_tmp(1:nVars)=matmul(matLHS(5,i,j,k,1:nVars,1:nVars),data(i,j+1,k,1:nVars))
             rhstmp(i,j,k,1:nVars)=rhstmp(i,j,k,1:nVars)-vec_tmp(1:nVars)                   
             !
#if (iins3d==1)
             ! i,j,k-1
             vec_tmp(1:nVars)=matmul(matLHS(6,i,j,k,1:nVars,1:nVars),data(i,j,k-1,1:nVars))
             rhstmp(i,j,k,1:nVars)=rhstmp(i,j,k,1:nVars)-vec_tmp(1:nVars)
             ! i,j,k+1
             vec_tmp(1:nVars)=matmul(matLHS(7,i,j,k,1:nVars,1:nVars),data(i,j,k+1,1:nVars))
             rhstmp(i,j,k,1:nVars)=rhstmp(i,j,k,1:nVars)-vec_tmp(1:nVars)                   
#endif
             !
          end do
       end do
    end do
    !
  end subroutine evaluateRHS
  !
  subroutine exchange_ghosts
    !
    use dimensions, only: BCregions,BCType
    !
    implicit none
    ! local
    Integer    :: idir,fb,bid,dd
    Integer    :: i1,i2,j1,j2,k1,k2

    ! loop over domain ghosts
    !
    do idir=1,2+iins3d
       do fb=1,2
          ! boundary ID
          bid=fb+(idir-1)*2
          dd=-3+2*fb
          !
          i1=abs(min(dd*BCregions(bid,1,1),dd*BCregions(bid,1,2)))
          i2=abs(max(dd*BCregions(bid,1,1),dd*BCregions(bid,1,2)))
          j1=abs(min(dd*BCregions(bid,2,1),dd*BCregions(bid,2,2)))
          j2=abs(max(dd*BCregions(bid,2,1),dd*BCregions(bid,2,2)))
          k1=abs(min(dd*BCregions(bid,3,1),dd*BCregions(bid,3,2)))
          k2=abs(max(dd*BCregions(bid,3,1),dd*BCregions(bid,3,2)))
          !
          if (BCType(bid).eq.1) then
             !
             !
             call dirichletBC(sol(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars),nVars,i1,i2,j1,j2,k1,k2,dd)
             !
          else if (BCType(bid).eq.2) then
             !
             call neumannBC(sol(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nVars),nVars,i1,i2,j1,j2,k1,k2,dd,idir)
             !
          else
             !
             write(*,'(A,I3)') "[E] BC-type not implemented! BCtype = ",BCType(bid)
             stop
             !
          endif
          !
       enddo
       !
    enddo
    !
  end subroutine exchange_ghosts
  !
  subroutine dirichletBC(data,nv,i1,i2,j1,j2,k1,k2,dd)
    !
    implicit none
    !
    Real*8, intent(inout)   :: data(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nv)
    Integer, intent(in)     :: nv,i1,i2,j1,j2,k1,k2,dd
    ! local
    Integer                 :: i,j,k
    Real*8                  :: xyz(3)
    !
    do k=k1,k2,dd
       do j=j1,j2,dd
          do i=i1,i2,dd
             !
             xyz(1:3)=coords(i,j,k,1:3)
             !
             call userSpecifiedFunction(data(i,j,k,1:nv),nv,xyz(1:3),Lx,Ly,Lz)
             !
          end do
       end do
    end do
    !
  end subroutine dirichletBC
  !
    subroutine neumannBC(data,nv,i1,i2,j1,j2,k1,k2,dd,dir)
    !
    implicit none
    !
    Real*8, intent(inout)   :: data(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:nv)
    Integer, intent(in)     :: nv,i1,i2,j1,j2,k1,k2,dd,dir
    ! local
    Integer                 :: i,j,k
    Integer                 :: iloc,jloc,kloc
    Integer                 :: npoints,dirVec(3)
    Real*8                  :: xyz(3)
    !
    dirVec=0
    dirVec(dir)=1
    !
    do k=k1,k2,dd
       do j=j1,j2,dd
          do i=i1,i2,dd
             !
             npoints=nxb*dirVec(1)+nyb*dirVec(2)+nzb*dirVec(3)
             !                                                                                                                                                                                                
             iloc=(2*nguard+1+(1+dd)*(npoints)-i)*dirVec(1)+(1-dirVec(1))*i
             jloc=(2*nguard+1+(1+dd)*(npoints)-j)*dirVec(2)+(1-dirVec(2))*j
#if (iins3d==1)
             kloc=(2*nguard+1+(1+dd)*(npoints)-k)*dirVec(3)+(1-dirVec(3))*k
#else
             kloc=1
#endif
             !
             data(i,j,k,1:nv)=data(iloc,jloc,kloc,1:nv)
             !
          end do
       end do
    end do
    !
  end subroutine neumannBC
  !
end module LinearSolver
