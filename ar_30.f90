program ar
  use omp_lib 
  implicit none
  integer(8) :: N, i, idx, ic, is
  real(8) :: dt, A, k, B, F, e, kb
  real(8) :: t, C, t_max, x1,y1,t1
  integer :: temp, last,down,ns, seed_size,ir,ik,sx,sy,tx,ty,k_coef
  real(8) :: start_time, end_time
  real(8) :: x, y, ax, ay, xo, xn, yo, yn
  real(8) :: n_x, n_y,x_1,y_1,x2,y2
  character (len=20) :: filename,fn,fn1,fn2
  real(8) :: u1, u2, v1, v2,pi, dx, dy 
  integer :: count,V0, ios, ti, tc, sum, tid, c_x, c_y, cpu=36
  real(8), dimension(0:1000) :: bx, by
  integer, dimension(8) :: values
  integer, dimension(:), allocatable :: seed
  character(10) :: D1,D2,D3
  
  call date_and_time (D1,D2,D3,values)
  print '(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)', values(1), '-',  values(2),&
 '-', values(3), '     ', values(5), ':', values(6), ':', values(7)
  start_time=omp_get_wtime()

  ! Initial conditions
  pi=acos(-1.0)
  N=1000000000
  ns=50
  down=1000
!  last=N/down-2*ns
!  last=N/down
  dt=1.0e-14
!  k_coef=10               !meV/rad^2
  t_max=N*dt
  V0=30
!  temp=40
  e=1.60217663e-19
  kb=1.380649e-23
  A=1.0e-42
  B=5.0e-31                            !5.0e-31
  F=V0*e/1000
!  C=sqrt(2*B*kb*temp/dt)
  k_coef=0
!  ik=0
  do ik=1,4  
  

  call random_seed(size=seed_size)
  allocate(seed(seed_size))
  call omp_set_num_threads (cpu)
  !$omp parallel private(i,t,x2,y2,tid,u1,u2,v1,v2,n_x,n_y,ax,ay,x,y,xo,yo,xn,yn,x1,y1,x_1,y_1,fn,seed,count,dx,dy,c_x,c_y,k,temp,C)
  
  tid=omp_get_thread_num() + 1
  write(fn, '(A, I3.3, A)') 'ar_', k_coef+tid-1, '.txt'
  open(unit=600+tid, file=fn, status='replace')
  write(600+tid,'(A,E8.1,A,I5,A,I5,A,E8.1,A,E8.1)')'dt= ',dt, ' s   V0= ',V0,&
' mV   k= ',k_coef+tid-1, ' meV/rad^2   t_max= ',t_max, ' s   win= ',ns*dt    
  
  temp=70
  do while (temp<=140)
  C=sqrt(2*B*kb*temp/dt)
  call system_clock(count=seed(1))
  seed=seed+tid*101
  call random_seed(put=seed)
!  write (fn, '(A,I2.2,A)') 'rot_',tid,'.txt'
!  open (unit=100+tid, file=fn, status='replace')
!  write(fn1, '(A,I2.2,A)') 'smooth_',tid,'.txt' 
!  open (unit=300+tid, file=fn1, status='replace')
  
  k=(k_coef+tid-1)*e/1000
  x = 0.0
  y = 0.0
  xo= 0.0
  yo= 0.0
  x_1=0.0
  y_1=0.0
  x2=0.0
  y2=0.0
  count=0
  c_x=0
  c_y=0
  do i = 1, N
    call random_number(u1)
    call random_number(u2)
    call random_number(v1)
    call random_number(v2)
    n_x = sqrt(-2.0 * log(u1+1e-15)) * cos(2.0 * pi * u2)
    n_y = sqrt(-2.0 * log(v1+1e-15)) * cos(2.0 * pi * v2)
    
    ax = ( -(F*6*0.5*sin(6*x) + k*6*sin(6*x+6*y)) - B*(x-xo)/dt + C*n_x ) / A
    ay = ( -(F*6*0.5*sin(6*y) + k*6*sin(6*x+6*y)) - B*(y-yo)/dt + C*n_y ) / A    
    xn = 2.0d0*x - xo + ax*dt**2 
    yn = 2.0d0*y - yo + ay*dt**2
    
    x1=(pi/3.0)*nint(x*3.0/pi) 
    y1=(pi/3.0)*nint(y*3.0/pi)
    
    if (abs(x1-x)> 0.2) then
 !          d=sign(x1-x, 1.0d0) 
           x1=x2
    else
           x1=(pi/3.0)*nint(x*3.0/pi) 
           x2=x1 
    end if 
    if (abs (y1-y)> 0.2) then
            y1=y2
    else
            y1=(pi/3.0)*nint(y*3.0/pi)
            y2=y1
    end if

!    if (mod(i,down) == 0) then
!        write(100+tid,'(F16.13, F20.10, F20.10)') i*dt, x,  y
!        write(300+tid,'(F16.13, F20.10, F20.10)') i*dt, x1, y1
!    end if
    
!    if (mod(i,100)==0) then
!            dx=x1 - x_1
!            dy=y1 - y_1
!            if (abs(dx/1.047)>0.6 .and. abs(dy/1.047)>0.6 .and. dx*dy<0.0) then
!                count=count+1
!            end if
!            x_1=x1
!            y_1=y1
!    end if

    if (mod(i,ns)==0) then
            dx=x1 - x_1
            dy=y1 - y_1
            if (abs(dy/1.047)>0.6) c_y=c_y+1
            if (abs(dx/1.047)>0.6) then
                c_x=c_x+1
                if (abs(dy/1.047)>0.6 .and. dx*dy<0.0) count=count+1
            end if
            x_1=x1
            y_1=y1
    end if
    
    xo = x
    yo = y
    x = xn
    y = yn

  end do
  write(600+tid,'(I5,I12,I12,I12)') temp,count,c_x,c_y 
  temp=temp+10
  end do
  flush(600+tid)
  close(600+tid) 

!  flush(100+tid)
!  close(100+tid)
!  flush(300+tid)
!  close(300+tid)
 

  !$omp end parallel
  k_coef=k_coef+cpu
  deallocate(seed)
  end do

  end_time=omp_get_wtime()
  call date_and_time(D1,D2,D3,values)
  print '(I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)',&
values(1),'-',values(2),'-',values(3),&
         '     ',values(5),':',values(6),':',values(7) 

  print'(A,F15.1,A)', 'Total execution time:', end_time - start_time, ' seconds'
end program ar

