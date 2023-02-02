      PROGRAM test
      implicit none
      integer i, n
      real, allocatable, dimension (:) :: a
      namelist /input_deck/a
      n = 5
      allocate (a(n))
      open(unit=1,file='test.dat',status='old')
      read(unit=1,nml=input_deck)
      close(unit=1)

      do i = 1, n
        print*, a(i)
      end do
      END PROGRAM
