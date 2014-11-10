!****** Dimension of Matrix: *******************************************
!       NT=10:  NMAX=165  
!       NT=12:  NMAX=234  
!       NT=14:  NMAX=316  
!       NT=16:  NMAX=408  
!       NT=18:  NMAX=513  
!       NT=20:  NMAX=630
!       NT=22:  NMAX=759
#ifndef NMAX
#define NMAX 900
#endif

#ifdef __hpux
#define fdate fdate_
#endif

!******** limit for CPU time in sec, for C routine clock(): *************
#ifdef __hpux
!******  medium queue:   4 h = 14000 sec
!******  long   queue:  16 h = 57600 sec
#define CPUTIMELIM         52000.0D0
#else
#define CPUTIMELIM         40000.0D0
#endif
#define CLOCKS_PER_SECOND  1000000.0D0
#define NCPUMAXVALUE       2147483647

!****** DQS 3.0: Signal, when process exceeds CPU time limit: ********** 
!****** (doesn't work actually, because signal is send only to 
!******  shell, not to the working process.
!******  Only the Bourneshell can catch the signal, but reacts
!******  not before one shell command is finished (especially the
!******  working process).
#define SIGRESTART  80

!******  defines from <sys/signal.h> *********************
#  define SIGABRT       6       /* Process abort signal */
#  define SIGUSR1       16      /* user defined signal 1 */
#  define SIGUSR2       17      /* user defined signal 2 */

!****** Process return values (2...255) ********************************
#define NO_INFILE      100
#define ERR_IN_INFILE  101
#define DIM_TO_SMALL   102
#define ERR_WRT_OUTFILE 103
#define NO_RA_FOUND    120
#define START_NEXT_RUN 190
#define FINISHED       199

#define DOUBLEMAXVALUE  1.0D306
#define VERBOSE 1
#define QUIET   0
