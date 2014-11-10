module parameters
  implicit none
  ! Adimensional parameters
  double precision:: TAU,Rt,Pt,C,Le,Rc,eta
  double precision:: ri, ro
  ! Dimensions
  integer::NT,M0,ND,LMIN,LD,NE, NSMAX, LCALC
  double precision:: LowerLimit,UpperLimit,StepSize,DRt,ABSE,RELE, DRc
contains
   subroutine setDefaults()
      implicit none
      NE    = 2
      LCALC = 2
      Rt   = 4.D3
      DRt  = Rt*0.1d0
      Rc   = 4.D3
      DRc  = Rc*0.1D0
      TAU  = 100.0D0
      LowerLimit  = TAU
      UpperLimit  = 10*LowerLimit
      Pt   = 1.D-1
      Le   = 1.D00
      ETA  = 0.4D0
      NT   = 3
      M0   = 6
      ABSE = 1.0d-6
      RELE = 1.0d-6
      NSMAX = 100
   end subroutine
end module parameters
