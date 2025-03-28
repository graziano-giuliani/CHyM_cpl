!-----------------------------------------------------------------------
!     CHyM main program
!-----------------------------------------------------------------------
program chymmain
!
!-----------------------------------------------------------------------
!     Used module declarations
!-----------------------------------------------------------------------
!
  use mod_chym_iface
  use mod_chym_param, only : nstep, pstep
!
  implicit none

!-----------------------------------------------------------------------
!     Local variable declarations
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!     Initialize
!-----------------------------------------------------------------------
!
  call chym_init

  call chym_start(1,1)
!
!-----------------------------------------------------------------------
! Run
!-----------------------------------------------------------------------
!
  if (pstep == 0) pstep = pstep+1
  call chym_run(pstep, nstep, .false., 1, 1)

end program chymmain
