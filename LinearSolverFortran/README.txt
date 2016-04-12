Routines that can be ignored for point solve algorithm:

subroutine blk_tridiag_solve(bm,bd,bp,x,b,jmax,neqs)
subroutine lu_solve(A,x,b,neqs)
subroutine lu(A,L,U,neqs)
subroutine upper_tri_solve(U,x,b,neqs)
subroutine lower_tri_solve(L,x,b,neqs)
subroutine matvec(A,b,c,n)
subroutine matmat(A,B,C,n)
