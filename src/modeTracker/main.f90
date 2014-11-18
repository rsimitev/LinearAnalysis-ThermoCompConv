!> Reads a file with the number of modes and a list of files
!! containing said modes for different parameters,
!! one file per line. For each file, outputs a reordered set
!! of values based on a quadratic extrapolation through the 
!! previous 3 files.
program modeTracker
   implicit none
   integer, parameter:: unitFileList=15, nPoints=3
   character(len=128):: listFileName, fileToRead
   character(len=256):: line
   double precision:: par(nPoints)
   double complex, allocatable:: modes(:,:), projectedModeVals(:)
   double precision:: newPar
   double complex, allocatable:: newModeSet(:)
   integer:: nModes !< the number of modes to be read from each file
   integer:: nFile !< The index of the modes file read
   integer:: info, i, j


   call init(info)

   ! We need to read one file at a time and recompute the fit.
   do
      call readLineOfFilesToRead(unitFileList, line, info)
      if(info.ne.0) exit
      nFile = nFile + 1
      ! Shift the parameter and mode pointers to the next value
      read(line,*) newPar, fileToRead
      Write(*,*) "Processing ", trim(fileToRead)
      ! Read one file at a time
      call readModesFile(trim(fileToRead), newModeSet)
      ! For each mode, we compute the projected mode value
      do i=1, nModes
         projectedModeVals(i) = nextPoint3(par, modes(i,:), newPar)
      enddo
      ! Now that we have the next parameter value and the
      ! next set of modes, we reorder it.
      do i=1, nModes
         ! Reorder the newly read mode set according to the distance to the
         ! projected value.
         do j=1, nModes
            if(abs(projectedModeVals(i)-newModeSet(j)).lt. &
               abs(projectedModeVals(i)-newModeSet(i))) then
               ! if j>i always swap
               if (j>i) then
                  call swap(newModeSet(i),newModeSet(j))
               else !if j<i, swap only when it makes a difference to both.
                  if(abs(projectedModeVals(j)-newModeSet(i)).lt. &
                     abs(projectedModeVals(j)-newModeSet(j))) then
                     call swap(newModeSet(i),newModeSet(j))
                  endif
               endif
            endif
         enddo
      enddo
      par(mod(nFile-1,nPoints)+1)     = newPar
      modes(:,mod(nFile-1,nPoints)+1) = newModeSet(:)
      ! Now that things are in proper order, we need to save the new set of modes to file
      call saveReorderedModeSet(newModeSet, nModes)
      ! and then add it to the modes that are going to be used for interpolation
      ! in the next iteration.
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

      open(unit=unitFileList, file=trim(listFileName), status='OLD', iostat=info)
      if (info.ne.0) call MTabort(info)

      ! Read the number of modes
      read(unitFileList,*) nModes

      allocate(modes(nModes,nPoints),newModeSet(nModes))
      allocate(projectedModeVals(nModes))
      modes = 0.0d0

      ! Read three files so we can make the first quadratic fit.
      nFile = 0
      do
         call readLineOfFilesToRead(unitFileList, line, info)
         if(info.ne.0) exit
         ! Read one file at a time
         nFile = nFile + 1
         read(line,*) par(nFile), fileToRead
         Write(*,*) "Processing ", trim(fileToRead)
         call readModesFile(fileToRead, modes(:,nFile))
         call saveReorderedModeSet(modes(:,nFile), nModes)
         if (nFile==nPoints) exit
      enddo
   end subroutine

   !**********************************************************************
   !> Reads a file one line at a time until it finds one that
   !! does not start with a '#'.
   !**********************************************************************
   subroutine readLineOfFilesToRead(unitFileList,line, info)
      implicit none
      character(len=*), intent(out):: line
      integer, intent(in):: unitFileList
      integer, intent(out):: info
      info=0
      do
         read(unitFileList,'(A)', iostat=info) line
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
      double complex, intent(out):: modes(:)
      character(len=128):: line
      double precision:: modeVal, modeFreq
      integer:: n, readStatus

      readStatus = 0
      modes = 0.0d0
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
         if (n.gt.nModes) exit
         modes(n) = cmplx(modeFreq, modeVal)
      enddo

      close(unit=20)
   end subroutine

   !**********************************************************************
   !> Given three points of coordinates xOld, yOld, extimate the next point at x.
   !**********************************************************************
   function nextPoint3(xOld,yOld,x) result(y)
      implicit none
      double complex:: y
      double precision, intent(in):: xOld(nPoints), x
      double complex, intent(in):: yOld(nPoints)
      double precision:: d, x1, x2, x3
      double complex:: a,b,c,y1, y2, y3

      x1 = xOld(1)
      x2 = xOld(2)
      x3 = xOld(3)
      y1 = yOld(1)
      y2 = yOld(2)
      y3 = yOld(3)

      d = (x2-x1)*(x3-x1)*(x3-x2)
      a = ( (x2-x1)*y3 - (x3-x1)*y2 + (x3-x2)*y1 )/d
      b =-( (x2**2-x1**2)*y3 - (x3**2-x1**2)*y2 + (x3**2-x2**2)*y1 )/d
      c = ( (x2 - x1)*x1*x2*y3 - (x3 - x1)*x1*x3*y2 + (x3 - x2)*x2*x3*y1 )/d

      y = a*x**2 + b*x + c
   end function


   !**********************************************************************
   !> Swap a with b
   !**********************************************************************
   subroutine swap(a,b)
      implicit none
      double complex, intent(inout):: a,b
      double complex:: c
      c=a
      a=b
      b=c
   end subroutine

   !**********************************************************************
   !>
   !**********************************************************************
   subroutine saveReorderedModeSet(modeSet, nModes)
      implicit none
      integer, intent(in):: nModes
      double complex, intent(in):: modeSet(nModes)
      integer:: i
      open(unit=33, file=trim(fileToRead)//'-sorted', status='NEW')
      do i=1, nModes
         Write(33,*) i, dble(modeSet(i)), dimag(modeSet(i))
      enddo
      close(unit=33)
   end subroutine

   !**********************************************************************
   !> Save, cleanup and shut down.
   !**********************************************************************
   subroutine cleanup()
      implicit none
      close(unit=unitFileList)
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
