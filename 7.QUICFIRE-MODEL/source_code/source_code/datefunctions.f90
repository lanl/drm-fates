      subroutine dateToUnix(idate, utime)
         implicit none
         integer(8),intent(in)  :: idate(6)
         integer(8),intent(out) :: utime
         integer(8) :: i,days(12)
			integer(8),parameter ::			&
				four_hundred = int(400,8), &
				one_hundred = int(100,8),	&
				four = int(4,8)
      	
         days = (/31,28,31,30,31,30,31,31,30,31,30,31/)
         !initially, set the timestamp to be zero
         utime = 0
         
         !determine if the current year is a leap year, which will potentially modify the 
         !number of days in February
         if(mod(idate(1), four_hundred) == 0) then
            days(2) = 29
         elseif(mod(idate(1), one_hundred) == 0) then
            days(2) = 28
         elseif(mod(idate(1), four) == 0) then 
            days(2) = 29
         endif
         
         !then check if the input values are valid
         if(idate(1) .lt. 1901 .or. idate(1) .gt. 2038)then
            print*, "dateToUnix: invalid year specified", idate(1)
            return
         elseif(idate(2) .lt. 1 .or. idate(2) .gt. 12)then
            print*, "dateToUnix: invalid month specified", idate(2)
            return
         elseif(idate(3) .lt. 1 .or. idate(3) .gt. days(idate(2)))then
            print*, "dateToUnix: invalid day specified", idate(3)
            return
         elseif(idate(4) .lt. 0 .or. idate(4) .gt. 23)then
            print*, "dateToUnix: invalid hour specified", idate(4)
            return
         elseif(idate(5) .lt. 0 .or. idate(5) .gt. 59)then
            print*, "dateToUnix: invalid minute specified", idate(5)
            return
         elseif(idate(6) .lt. 0 .or. idate(6) .gt. 59)then
            print*, "dateToUnix: invalid second specified", idate(6)
            return
         endif
         
         !determine if the specified date is before or after the beginning of the epoch
         if(idate(1) .lt. 1970)then
            !then calculate the number of seconds between jan 1, 1970, 12:00am and the beginning
            !of the specified year, moving backwards (before beginning of epoch)
            do i=1970, idate(1) + 1, -1
               !leap years are years which are multiples of 400 or years that are multiples of 
               ! 4 and not 100
               if(mod(i,four_hundred) == 0) then
                  utime = utime - 31622400 !366 * 24 * 60 * 60 = 31622400s
               elseif(mod(i,one_hundred) == 0) then
                  utime = utime - 31536000 !365 * 24 * 60 * 60 = 31536000s
               elseif(mod(i,four) == 0) then 
                  utime = utime - 31622400 !366 * 24 * 60 * 60 = 31622400s
               else
                  utime = utime - 31536000 !365 * 24 * 60 * 60 = 31536000s
               endif
            enddo
            if(days(2) == 29)then
               utime = utime - 86400 ! on a leap year, the beginning of the year is a day further away from 0
            endif
         else
            !then calculate the number of seconds between jan 1, 1970, 12:00am and the beginning
            !of the specified year (after beginning of epoch)
            do i=1970, idate(1) - 1
               !leap years are years which are multiples of 400 or years that are multiples of 
               ! 4 and not 100
               if(mod(i, four_hundred) == 0) then
                  utime = utime + 31622400 !366 * 24 * 60 * 60 = 31622400s
               elseif(mod(i, one_hundred) == 0) then
                  utime = utime + 31536000 !365 * 24 * 60 * 60 = 31536000s
               elseif(mod(i, four) == 0) then 
                  utime = utime + 31622400 !366 * 24 * 60 * 60 = 31622400s
               else
                  utime = utime + 31536000 !365 * 24 * 60 * 60 = 31536000s
               endif
            enddo
         endif
         
         !then calculate the number of seconds between the beginning of the month specified
         !and the beginning of the year specified
         do i=1,idate(2) - 1
            utime = utime + days(i) * 86400 ! 24 * 60 * 60 = 86400s
         enddo
         
         !From here on out, everything should be standard.  no leap year stuff to deal with, 
         !and leap seconds are not accounted for in UNIX time
         utime = utime + (idate(3) - 1) * 86400 ! 24 * 60 * 60 = 86400s
         utime = utime + idate(4) * 3600 ! 60 * 60 = 3600s
         utime = utime + idate(5) * 60
         utime = utime + idate(6) 
      end subroutine
      
      ! unix2c Converts Unix system time to date/time integer array.
      subroutine unix2c(utime, idate)
         implicit none
         integer(8) utime, idate(6), days(12)
         ! utime  input  Unix system time, seconds since 1970.0
         ! idate  output Array: 1=year, 2=month, 3=date, 4=hour, 5=minute, 6=secs
         ! -Author  Clive Page, Leicester University, UK.   1995-MAY-2
         integer mjday, nsecs
         real day
         ! Note the MJD algorithm only works from years 1901 to 2099.
         mjday    = int(utime/86400 + 40587)
         idate(1) = 1858 + int( (mjday + 321.51) / 365.25)
         day      = aint( mod(mjday + 262.25, 365.25) ) + 0.5
         idate(2) = 1 + int(mod(day / 30.6 + 2.0, 12.0) )
         idate(3) = 1 + int(mod(day,30.6))
         nsecs    = int(mod(utime, int(86400,8)),4)
         idate(6) = mod(nsecs, 60)
         nsecs    = nsecs / 60
         idate(5) = mod(nsecs, 60)
         idate(4) = nsecs / 60
         
         days = (/31,28,31,30,31,30,31,31,30,31,30,31/)
         !initially, set the timestamp to be zero
         utime = 0
         
         !determine if the current year is a leap year, which will potentially modify the 
         !number of days in February
         if(mod(idate(1),int(400,8)) == 0) then
            days(2) = 29
         elseif(mod(idate(1),int(100,8)) == 0) then
            days(2) = 28
         elseif(mod(idate(1),int(4,8)) == 0) then 
            days(2) = 29
         endif
         
         !normalize seconds
         do while(idate(6) .lt. 0)
            idate(5) = idate(5) - 1
            idate(6) = idate(6) + 60
         enddo
         do while(idate(6) .gt. 59)
            idate(5) = idate(5) + 1
            idate(6) = idate(6) - 60
         enddo
         !normalize minutes
         do while(idate(5) .lt. 0)
            idate(4) = idate(4) - 1
            idate(5) = idate(5) + 60
         enddo
         do while(idate(5) .gt. 59)
            idate(4) = idate(4) + 1
            idate(5) = idate(5) - 60
         enddo
         !normalize hours
         do while(idate(4) .lt. 0)
            idate(3) = idate(3) - 1
            idate(4) = idate(4) + 24
         enddo
         do while(idate(4) .gt. 23)
            idate(3) = idate(3) + 1
            idate(4) = idate(4) - 24
         enddo
         !normalize days
         !note that day is a special case.  It can increment or decrement month, and in order
         !to get the number of days to add/subract, we need to know the normalized month.  As
         !such, the month has to be normalized in the middle of normalizing the day
         do while(idate(3) .lt. 1)
            idate(2) = idate(2) - 1
            !normalize month
            do while(idate(2) .lt. 1)
               idate(1) = idate(1) - 1
               idate(2) = idate(2) + 12
            enddo
            do while(idate(2) .gt. 12)
               idate(1) = idate(1) + 1
               idate(2) = idate(2) - 12
            enddo
            idate(3) = idate(3) + days(idate(2))
         enddo
         do while(idate(3) .gt. days(idate(2)))
            idate(2) = idate(2) + 1
            !normalize month
            do while(idate(2) .lt. 1)
               idate(1) = idate(1) - 1
               idate(2) = idate(2) + 12
            enddo
            do while(idate(2) .gt. 12)
               idate(1) = idate(1) + 1
               idate(2) = idate(2) - 12
            enddo
            idate(3) = idate(3) - days(idate(2))
         enddo
         !normalize month 
         do while(idate(2) .lt. 1)
            idate(1) = idate(1) - 1
            idate(2) = idate(2) + 12
         enddo
         do while(idate(2) .gt. 12)
            idate(1) = idate(1) + 1
            idate(2) = idate(2) - 12
         enddo
      end
