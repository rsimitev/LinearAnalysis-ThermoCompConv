module parameters
  implicit none
  ! Adimensional parameters
  double precision:: TAU,RA,PR,RI,C,pL,Rconc,eta
  ! Dimensions
  integer::NT,M0,ND,LMIN,LD,NE, NSMAX, LCALC
  double precision:: TTA,TTF,TTSTEP,DRA,ABSE,RELE
contains
   subroutine setDefaults()
      implicit none
      NE    = 2
      LCALC = 2
      RA       = 4.D3
      RAOLD    = 4.D3
      DRA      = RA*0.1d0
      Rconc    = 4.D3
      RconcOLD = 4.D3
      DRconc   = Rconc*0.1D0
      TAU  = 100.0D0
      TAU1 = 100.0D0
      TAU2 = 100.0D0
      TTA  = TAU
      TTF  = 10*TTA
      PR   = 1.D-1
      pL   = 1.D00
      ETA  = 0.4D0
      NT   = 3
      M0   = 6
      ABSE = 0.D0
      RELE = 1.D-6
      NSMAX = 100
   end subroutine
end module parameters
