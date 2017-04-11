!!$program test_rand
!!$  use ifport
!!$  integer,parameter :: seed = 86456
!!$  
!!$  call srand(seed)
!!$  print *, rand(), rand(), rand(), rand()
!!$  print *, rand(seed), rand(), rand(), rand()
!!$end program test_rand

  program test_rand
    use ifport
    
  print *, rand(), rand(), rand(), rand()
  print *, rand(), rand(), rand()
end program test_rand
