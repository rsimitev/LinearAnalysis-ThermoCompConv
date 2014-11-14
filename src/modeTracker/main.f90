!> Reads a file with the number of modes and a list of files 
!! containing said modes for different parameters, 
!! one file per line.
program modeTracker
   implicit none
   integer, parameter:: unitOut=16, unitIn=15
   character(len=128):: listFileName, fileToRead
   character(len=256):: line
   double precision, target:: par(3)
   double precision, allocatable, target:: modes(:,:)
   double precision:: newPar, projectedModeVal
   double precision, pointer:: newModeSet(:)
   double precision, pointer:: modesPtr(:), parPtr
   integer:: nModes !< the number of modes to be read from each file
   integer:: nFile !< The index of the modes file read
   integer:: info, i, j


   call init(info)

   ! We need to read one file at a time and recompute the fit.
   do
      call readLineOfFilesToRead(unitIn, line, info)
      if(info.ne.0) exit
      Write(*,*) line
      nFile = nFile + 1
      ! Shift the parameter and mode pointers to the next value
      parPtr => par(mod(nFile,3)+1)
      modesPtr => modes(:,mod(nFile,3)+1)
      read(line,*) newPar, fileToRead
      ! Read one file at a time
      call readModesFile(trim(fileToRead), newModeSet)
      ! Now that we have the next parameter value and the
      ! next set of modes, we reorder it.
      do i=1, nModes-1
         ! For each mode, we compute the projected mode value
         projectedModeVal = nextPoint(par,modes(i,:),newPar)
         ! Reorder the newly read mode set according to the distance to the
         ! projected value
         do j=i+1, nModes
            if(abs(projectedModeVal-newModeSet(j)).lt.abs(projectedModeVal-newModeSet(i))) then
               call swap(newModeSet(i),newModeSet(j))
            endif
         enddo
      enddo
      ! Now that things are in proper order, we need to save the new set of modes to file
      call saveReorderedModeSet(newModeSet)
      ! and then add it to the modes that are going to be used for interpolation 
      ! in the next iteration.
      parPtr   = newPar
      modesPtr = newModeSet
   enddo

   call cleanup()
contains

   !**********************************************************************
   !> Initialize stuff
   !**********************************************************************
   subroutine init(info)
      implicit none
      integer:: info

      call getarg(1,listFileName)
      if (listFileName.eq.' ' .or. listFileName.eq.'-h') then
         print*, 'Usage : '
         print*, 'modeTracker <list file>'
         stop
      endif

      open(unit=unitIn, file=trim(listFileName), status='OLD', iostat=info)
      if (info.ne.0) call MTabort(info)

      ! Read the number of modes
      read(unitIn,*) nModes

      allocate(modes(nModes,3),newModeSet(nModes))
      modes = 0.0d0
      par = 3.0d0
      par(1)=5.0

      ! Read three files so we can make the first quadratic fit.
      nFile = 0
      do
         call readLineOfFilesToRead(unitIn, line, info)
         if(info.ne.0) exit
         Write(*,*) line
         ! Read one file at a time
         nFile = nFile + 1
         parPtr => par(nFile)
         modesPtr => modes(:,nFile)
         read(line,*) parPtr, fileToRead
         Write(*,*) nFile, parPtr, fileToRead
         call readModesFile(fileToRead, modesPtr)
         call saveReorderedModeSet(modesPtr)
         if (nFile==3) exit
      enddo
   end subroutine

   !**********************************************************************
   !> Reads a file one line at a time until it finds one that
   !! does not start with a '#'.
   !**********************************************************************
   subroutine readLineOfFilesToRead(unitIn,line, info)
      implicit none
      character(len=*), intent(out):: line
      integer, intent(in):: unitIn
      integer, intent(out):: info
      info=0
      do
         read(unitIn,'(A)', iostat=info) line
         if (info.ne.0) exit
         line = adjustl(line)
         if (line(1:1).ne.'#') exit
      enddo
   end subroutine

   !**********************************************************************
   !> Read all the growth rates from one file.
   !**********************************************************************
   subroutine readModesFile(fileToRead, modes)
      implicit none
      character(len=*), intent(in):: fileToRead
      double precision, pointer, intent(out):: modes(:)
      character(len=128):: line
      double precision:: modeVal, modeFreq
      integer:: n, readStatus

      readStatus = 0
      open(unit=20, file=fileToRead, status='OLD')

      do
         ! Read one line at a time
         read(20,'(A)', iostat=readStatus) line
         ! Exit on error
         if(readStatus.ne.0) exit
         ! Check if this is a comment and cycle if it is.
         line = adjustl(line)
         if (line(1:1).eq.'#') cycle
         ! Read the mode index and the value
         read(line,*) n, modeFreq, modeVal
         modes(n) = modeVal
      enddo

      close(unit=20)
   end subroutine

   !**********************************************************************
   !> Given three points of coordinates xOld, yOld, extimate the next point at x.
   !**********************************************************************
   function nextPoint(xOld,yOld,x) result(y)
      implicit none
      double precision:: y
      double precision, intent(in):: xOld(3), yOld(3), x
      double precision:: a,b,c, x1, x2, x3, y1, y2, y3

      x1 = xOld(1)
      x2 = xOld(2)
      x3 = xOld(3)
      y1 = yOld(1)
      y2 = yOld(2)
      y3 = yOld(3)

      A = ((Y2-Y1)*(X1-X3) + (Y3-Y1)*(X2-X1))/((X1-X3)*(X2**2-X1**2) + (X2-X1)*(X3**2-X1**2))
      B = ((Y2 - Y1) - A*(X2**2 - X1**2)) / (X2-X1)
      C = Y1 - A*X1**2 - B*X1

      y = a*x**2 + b*x + c
   end function

   !**********************************************************************
   !> Swap a with b
   !**********************************************************************
   subroutine swap(a,b)
      implicit none
      double precision, intent(inout):: a,b
      double precision:: c
      c=a
      a=b
      b=c
   end subroutine
   
   !**********************************************************************
   !> 
   !**********************************************************************
   subroutine saveReorderedModeSet(modeSet)
      implicit none
      double precision, intent(in):: modeSet(:)
      integer:: i
      open(unit=33, file=trim(fileToRead)//'-sorted', status='NEW')
      do i=1, nModes
         Write(33,*) i, modeSet(i)
      enddo
      close(unit=33)
   end subroutine

   !**********************************************************************
   !> Save, cleanup and shut down.
   !**********************************************************************
   subroutine cleanup()
      implicit none
      close(unit=unitIn)
   end subroutine

   !**********************************************************************
   !> Cleanly aborts operations Writing an error message.
   !**********************************************************************
   subroutine MTabort(info)
      implicit none
      integer, intent(in):: info
      Write(*,*) "Aborting: error = ", info
      stop 1
   end subroutine

end program modeTracker
! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
