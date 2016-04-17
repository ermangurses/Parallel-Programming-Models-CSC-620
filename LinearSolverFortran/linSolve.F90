#include "setup.h"
program linSolve
  !
  use omp_lib
  !use FTTimerClass
  use dimensions
  use LinearSolver, only: sol,assignSol,assignRHS,assignLHS,solve,solveLinearSystem_ini,&
                          error,errorOutput,evaluateRHS,rhs
  !
  implicit none
  !
  !TYPE(FTTimer)  ::     Timer
  Real*8         ::     fstart, fend
  Real*8         ::     ostart,oend  
  Integer        ::     k,ioerror,scheme
  Integer        ::     maxIterations
  Real*8         ::     errorTolerance
  Real*8         ::     relaxationFactor
  Character*120  ::     filename,arg
  !
  ! define domain
  !
  nxb=1
  nyb=1
  nzb=1
  numVars=4
  !
  il_bnd=1
  jl_bnd=1
  kl_bnd=1
  iu_bnd=1
  ju_bnd=1
  ku_bnd=1
  !
  Lx=1.d0
  Ly=1.d0
  Lz=1.d0
  !
  dx=1.d0
  dy=1.d0
  dz=1.d0
  !
  if (iargc().gt.0) then
    DO k = 1, iargc()
       CALL getarg(k, arg)
    END DO
    filename=trim(arg)
 else
    write(filename,'(A)') 'input.dat'
 endif
 !
 write(*,'(A,A)') "[I] Read input file = ",trim(filename)
 open(unit=99, file=trim(filename), status='old', iostat=ioerror)
 !
 read(99,*) scheme
   write(*,'(A,I2)')    "[I] Scheme          = ",scheme
  read(99,*) BCType(1:6)
  write(*,'(A,6I2)')    "[I] BCType          = ",BCType(1:6)
  read(99,*) numVars                         
  write(*,'(A,I2)')     "[I] numVars         = ",numVars
#if (iins3d==0)                              
  !                                          
  read(99,*) Lx,Ly                           
  write(*,'(A,2E13.7)') "[I] Lx,Ly           = ",Lx,Ly
  read(99,*) nxb,nyb                         
  write(*,'(A,2I2)')    "[I] nxb,nyb         = ",nxb,nyb
#else                                        
  read(99,*) Lx,Ly,Lz                        
  write(*,'(A,3E13.7)') "[I] Lx,Ly,Lz        = ",Lx,Ly
  read(99,*) nxb,nyb,nzb                     
  write(*,'(A,3I2)')    "[I] nxb,nyb,nzb     = ",nxb,nyb
#endif                                       
  read(99,*) nguard                          
  write(*,'(A,I2)')     "[I] numVars         = ",nguard
  read(99,*) maxIterations                   
  write(*,'(A,I6)')     "[I] maxIterations   = ",maxIterations
  read(99,*) errorTolerance                  
  write(*,'(A,E13.7)')  "[I] errorTolerance  = ",errorTolerance
  read(99,*) relaxationFactor                
  write(*,'(A,E13.7)')  "[I] relaxationFact  = ",relaxationFactor
  close(99)
  !
  filename="post"
  call mkdir(filename)
  !
  ! define coordinates
  call computeCoordinates
  ! initialize linear system
  write(*,'(A)') "[I] Initialize linear system"
  call solveLinearSystem_ini(maxIterations,errorTolerance,relaxationFactor,&
                             numVars,nxb,nyb,nzb)
  !
  ! assign discretization matrix (left-hand-side)
  write(*,'(A)') "[I] Assign left-hand-side"
  call assignLHS
  !
  ! assign right-hand-side
  write(*,'(A)') "[I] Assign right-hand-side"
  call assignRHS
  !
  ! define initial condition
  call assignSol
  !
  filename="post/initialCondition.dat"
  !
  if (errorOutput) error=sol
  !
  write(*,'(A,A)') "[I] Writing output file = ",trim(filename)
  !
  call tecOutput(sol(il_bnd+nguard:iu_bnd-nguard,jl_bnd+nguard:ju_bnd-nguard,kl_bnd+nguard*iins3d:ku_bnd,1:numVars),&
                 coords(il_bnd+nguard:iu_bnd-nguard,jl_bnd+nguard:ju_bnd-nguard,kl_bnd+nguard*iins3d:ku_bnd-nguard*iins3d,1:3),&       
                 nxb,nyb,nzb,numVars,filename)
  !
  ! solve linear system
  write(*,'(A)') "[I] Solve linear system"
#if (RHSTEST==0)
  call omp_set_num_threads( omp_get_max_threads ( ) )
  write ( *, '(a,i8)' ) 'The number of processors available = ', omp_get_num_procs ( )
  write ( *, '(a,i8)' ) 'The number of threads available    = ', omp_get_max_threads ( )
!  call Timer % start()
   call cpu_time (fstart)
   ostart = omp_get_wtime() 
   !
   call solve(scheme)
   !
   call cpu_time (fend)
   oend = omp_get_wtime() 
   write(*,*) 'Fortran CPU time elapsed', fend-fstart
   write(*,*) 'OpenMP Walltime elapsed', oend-ostart 
  !
!  call Timer % stop()  
!  write(*,'(A,F17.8,A)') "[T] LinSolve took ",Timer%elapsedTime(TC_SECONDS) 
#endif
  !
  ! output results to tecplot file
  !
  if (errorOutput) then
     !
     error=error-sol
     filename="post/error.dat"
     !
     write(*,'(A,A)') "[I] Writing output file = ",trim(filename)
     !
     call tecOutput(error(il_bnd+nguard:iu_bnd-nguard,jl_bnd+nguard:ju_bnd-nguard,kl_bnd+nguard*iins3d:ku_bnd,1:numVars),&
                 coords(il_bnd+nguard:iu_bnd-nguard,jl_bnd+nguard:ju_bnd-nguard,kl_bnd+nguard*iins3d:ku_bnd-nguard*iins3d,1:3),&       
                 nxb,nyb,nzb,numVars,filename)
     !
#if (RHSTEST==1)
     call evaluateRHS(sol  (il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:numVars),&
                      error(il_bnd:iu_bnd,jl_bnd:ju_bnd,kl_bnd:ku_bnd,1:numVars))
     !
     filename="post/errorRHS.dat"
     !
     write(*,'(A,A)') "[I] Writing output file = ",trim(filename)
     !
     error=error+rhs
     !
     call tecOutput(error(il_bnd+nguard:iu_bnd-nguard,jl_bnd+nguard:ju_bnd-nguard,kl_bnd+nguard*iins3d:ku_bnd,1:numVars),&
                 coords(il_bnd+nguard:iu_bnd-nguard,jl_bnd+nguard:ju_bnd-nguard,kl_bnd+nguard*iins3d:ku_bnd-nguard*iins3d,1:3),&       
                 nxb,nyb,nzb,numVars,filename)
#endif
     !
  endif
  !
  filename="post/solution.dat"
  !
  write(*,'(A,A)') "[I] Writing output file = ",trim(filename)
  !
  call tecOutput(sol(il_bnd+nguard:iu_bnd-nguard,jl_bnd+nguard:ju_bnd-nguard,kl_bnd+nguard*iins3d:ku_bnd,1:numVars),&
                 coords(il_bnd+nguard:iu_bnd-nguard,jl_bnd+nguard:ju_bnd-nguard,kl_bnd+nguard*iins3d:ku_bnd-nguard*iins3d,1:3),&       
                 nxb,nyb,nzb,numVars,filename)
  !
end program linSolve
  
