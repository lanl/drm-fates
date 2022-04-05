      program rhofinfo

      parameter(n=200,m=200,l=41)
      real, dimension (1640000) :: livemoist
      real, dimension (1640000) :: moist
      real, dimension(100) :: a
      integer :: count1 
      integer :: count2
      integer :: count3
      integer :: x 
      integer :: y
      integer :: z
      integer :: ii
      integer :: nn
      character(len=40) ::file1
      character(len=40) ::file2
      character(len=16) :: str1
      character(len=4) :: str2
      character(len=4) :: str3
      character(len=16) :: str4
      character(len=4) :: str5

c      GETTING LIVE FULE MOISTURE FROM ORIGONAL treesmoist.dat
      open(1,file='treesmoist.dat',form='unformatted',
     +      status='old')
      read (1) livemoist 
      close(1)

       nn = 1
       str1 = 'MoistureContent.'
       str2 = '.txt'
       str4 = 'treesmoist.'
       str5 = '.dat'
        
      do nn=0,288 
        write (str3, '(I4)') nn
        str3 = adjustl(str3)

        file1 = adjustl(trim(str1))//adjustl(trim(str3))//str2
        file2 = adjustl(trim(str4))//adjustl(trim(str3))//str5
        file1 = adjustl(trim(file1))
        file2 = adjustl(trim(file2))
        file1 = trim(file1)
        file2 = trim(file2)
c        print *, trim(str1),str3, str2
        print *, file1, file2

      open(1,file=file1,form='formatted' 
     & ,status="old")
      do i=1,1640000
        read(1,*) moist(i)
c       Adding intercepted and live fuel moisture together 
        moist(i) = moist(i) + livemoist(i)    
        count1 = 1 + MOD(i,n)
        count2 = 1 + MOD((i/n),m)
        count3 = 1 + MOD((i/(n*m)),l)     
c        print *, count1, count2, count3, moist(i)
      enddo
      close(1)

      open (1,file=file2,form='unformatted' 
     & ,status='unknown')
      write (1) moist
      close (1)



      enddo


c      open(1,file="MoistureContent.19.txt",form='formatted' 
c     & ,status="old")
c      do i=1,1640000
c        read(1,*) moist(i)
c        count1 = 1 + MOD(i,n)
c        count2 = 1 + MOD((i/n),m)
c        count3 = 1 + MOD((i/(n*m)),l)     
cc        print *, count1, count2, count3, moist(i)
c      enddo
c      close(1)
c
c      open (1,file='treesmoist.19.dat',form='unformatted' 
c     & ,status='unknown')
c      write (1) moist
c      close (1)



      end
