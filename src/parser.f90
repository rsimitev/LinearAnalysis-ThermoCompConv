!     Copyright (C) 2005-2015 Luis Silva <lacsilva@gmail.com>
!
!     This program is free software: you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation, either version 3 of the License, or
!     (at your option) any later version.
!
!     This program is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!>   This module provides a simple parser for input files of the type
!!   'key = val'
!!   eventually separated by '[Sections]'
!! 
!!   Typical use of this module goes like:
!! ~~~~~~~ 
!! open(unit=444, file='myconfig.conf', status='old', iostat=error)
!! if (error.ne.0) return
!! do
!!    call parse(444, varname, line, error)
!!    select case(varname)
!!       case('var1')
!!          call read_val(line, var1)
!!       case('var2')
!!          call read_val(line, var2)
!!       case default
!!          cycle
!!    end select
!! enddo
!! close(444)
!! ~~~~~~
!
module parser
   implicit none
   private
   type vars
      character(len=15):: variable
      logical:: set
      logical:: required
   end type
   interface read_val
      module procedure read_val_int
      module procedure read_val_real
      module procedure read_val_double
      module procedure read_val_logical
      module procedure read_val_string
   end interface

   public:: read_val, parse, vars

contains

!> Reads one line from the specified unit and outputs the variable name and possible values as a string.
!! @par variable is the variable name.
!! @par line is the rhs of the attribution.
!! @par unit_in is the unit to read from
!! @par error is an error code.
subroutine parse(unit_in, variable, line, error)
   implicit none
   integer, intent(in):: unit_in
   character(len=*), intent(out):: variable
   character(len=256), intent(out):: line
   character(len=256):: line2
   character:: first
   integer, intent(out):: error
   integer:: eq
   variable = ""
   error    = 0
   do
      read(unit_in,'(A)',iostat=error) line
      if (error.ne.0) then
         return
      endif
      line2 = adjustl(line)
      first = line2(1:1)
      if((first=='#').or.(first=='*')) then
         cycle
      elseif(first=='[') then
         eq = scan(line2,"]")
         variable = 'SECTION'
         line = trim(adjustl(line2(1:eq-1)))
      else
         eq = scan(line2,"=")
         variable = trim(adjustl(line2(1:eq-1)))
         line = adjustl(line2(eq+1:))
         exit
      endif
   enddo
end subroutine parse
!***************************************************
subroutine read_val_int(line, value)
   implicit none
   character(len=256), intent(in):: line
   integer, intent(out):: value
   read(line,*) value
end subroutine read_val_int
!***************************************************
subroutine read_val_real(line, value)
   implicit none
   character(len=256), intent(in):: line
   real, intent(out):: value
   read(line,*) value
end subroutine read_val_real
!***************************************************
subroutine read_val_double(line, value)
   implicit none
   character(len=256), intent(in):: line
   double precision, intent(out):: value
   read(line,*) value
end subroutine read_val_double
!***************************************************
subroutine read_val_logical(line, value)
   implicit none
   character(len=256), intent(in):: line
   logical, intent(out):: value
   read(line,*) value
end subroutine read_val_logical
!***************************************************
subroutine read_val_string(line, value)
   implicit none
   character(len=256), intent(in):: line
   character(len=*), intent(out):: value
   integer:: spc
   spc = scan(line," ")
   value = line(:spc-1)
end subroutine read_val_string
end module

! vim: tabstop=3:softtabstop=3:shiftwidth=3:expandtab
