module parameters
  implicit none
  ! Adimensional parameters
  double precision:: TAU,Rt,Pt,RI,C,pL,Rc,eta
  ! Dimensions
  integer::NT,M0,ND,LMIN,LD,NE, NSMAX, LCALC
  double precision:: TTA,TTF,TTSTEP,DRt,ABSE,RELE, DRc
contains
   subroutine setDefaults()
      implicit none
      NE    = 2
      LCALC = 2
      Rt       = 4.D3
      DRt      = Rt*0.1d0
      Rc    = 4.D3
      DRc   = Rc*0.1D0
      TAU  = 100.0D0
      TTA  = TAU
      TTF  = 10*TTA
      Pt   = 1.D-1
      pL   = 1.D00
      ETA  = 0.4D0
      NT   = 3
      M0   = 6
      ABSE = 0.D0
      RELE = 1.D-6
      NSMAX = 100
   end subroutine
end module parameters
