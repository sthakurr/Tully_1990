program md
  implicit none
  real*8 pi,hbar
  complex*16 iota
  real*8 mass,AA,BB,CC,DD
  real*8 dt,curr_time,total_time
  integer nsteps
  integer,parameter :: nclass=1,nquant=2
  real*8 x(nclass),v(nclass),acc(nclass)
  complex*16 ci(nquant)
  integer state
  real*8 potential(nquant),energy
  real*8 dij(nquant,nquant,nclass)
  real*8 vdotd(nquant,nquant)
  real*8 psi(nquant,nquant),grad_H(nquant,nquant,nclass)
  real*8 transmission_lower,transmission_upper,reflection_lower,reflection_upper
  real*8 momentum

  call setup_parameters
  call main

  !-----------------------------------------------------------------  
  contains
  !-----------------------------------------------------------------  

  subroutine main
    implicit none
    integer i,j

    do i=1,30
      momentum= 1.d-1+30.d0*(i-1)/29.d0!10.d0 !2.d0+30.d0*(i-1)/19.d0
      transmission_lower=0.d0
      transmission_upper=0.d0
      reflection_lower=0.d0
      reflection_upper=0.d0

      do j=1, 500!2000
        call initial_conditions
        call evolve
      enddo

     call segregating_trajectories
 
    enddo

  end subroutine main
  !-----------------------------------------------------------------  
   subroutine segregating_trajectories
   implicit none
   
      transmission_lower=transmission_lower/500.d0
      transmission_upper=transmission_upper/500.d0
      reflection_lower=reflection_lower/500.d0
      reflection_upper=reflection_upper/500.d0

      write(15,*) momentum,transmission_lower
      write(16,*) momentum,transmission_upper
      write(17,*) momentum,reflection_lower
      write(18,*) momentum,reflection_upper

   end subroutine segregating_trajectories

 !-------------------------------------------------------------

  subroutine setup_parameters
    implicit none

    pi=dacos(-1.d0)
    iota=(0.d0,1.d0)
    hbar=1.d0

    mass=2000.d0

    AA=0.01
    BB=1.6
    CC=0.005
    DD=1.0

    dt=5.d0
    total_time=1000*dt

    nsteps=int(total_time/dt)

  end subroutine setup_parameters
  !-----------------------------------------------------------------  

  subroutine initial_conditions
    implicit none

    x(1)=-10.d0
    !momentum=20.d0
    v(1)=momentum/mass

    ci(1)=1.d0
    ci(2)=0.d0

    state=1

    curr_time=0.d0

  end subroutine initial_conditions
  !-----------------------------------------------------------------  

  subroutine evolve
    implicit none
    integer i 
!    write(6,*) nsteps,dt,total_time, ci
    call compute_acceleration
    call compute_energy

!call check_acceleration
!call draw_pes

    do i=1,nsteps
      write(10,*) curr_time,x,v,energy
      write(11,*) curr_time,abs(ci)**2,sum(abs(ci)**2)
      write(12,*) curr_time,state
      write(13,*) curr_time,vdotd(1,2)
      write(14,*) x,potential,dij(1,2,1)

      call evolve_1step_classical  !! evolves x,v
      !! Termintion condition - if x(1)<-10 or x(1)>10.d0  then terminate
      call evolve_1step_quantum    !! evolves ci
      call hop                     !! evolves state
      if(abs(x(1))>10.01d0)exit
    enddo

    !! Check where trajectory ended up
    !! if it ends with x<-10 and state=1, then 
    call checking_for_trajectory    

  end subroutine evolve
  !----------------------------------------------------------------  
   subroutine checking_for_trajectory
   implicit none

    if(x(1)>5.d0.and.state==1) transmission_lower=transmission_lower+1.d0
    if(x(1)>5.d0.and.state==2) transmission_upper=transmission_upper+1.d0
    if(x(1)<-5.d0.and.state==1) reflection_lower=reflection_lower+1.d0
    if(x(1)<-5.d0.and.state==2) reflection_upper=reflection_upper+1.d0

   end subroutine checking_for_trajectory
!--------------------------------------------------------------------------------

  subroutine evolve_1step_classical
    !! Evolving using velocity verlet
    implicit none

    x=x+v*dt+0.5*acc*dt*dt
    v=v+0.5*acc*dt
    call compute_acceleration
    v=v+0.5*acc*dt
    call compute_energy
    curr_time=curr_time+dt

  end subroutine evolve_1step_classical
  !-----------------------------------------------------------------  

  subroutine evolve_1step_quantum
    !! i hbar cj dot = Vj cj - i hbar/sum_k (v*d)_jk ck
    !! https://en.wikipedia.org/wiki/Runge–Kutta_methods
    implicit none
    complex*16 k1(nquant),k2(nquant),k3(nquant),k4(nquant)

    call compute_vdotd

    k1=derivative(ci)
    k2=derivative(ci + dt*k1/2.d0)
    k3=derivative(ci + dt*k2/2.d0)
    k4=derivative(ci + dt*k3)

    ci = ci + 1/6.d0 * dt *(k1+2*k2+2*k3+k4)
!print *, 'ci', ci
  end subroutine evolve_1step_quantum
  !-----------------------------------------------------------------  

  subroutine hop
    implicit none
    complex*16 a_kj(nquant)
    real*8 g_kj(nquant),b_kl(nquant)
    real*8 rnd,prob
    integer i

    a_kj = ci*conjg(ci(state))
    b_kl = -2*real(conjg(a_kj)*vdotd(:,state))
    g_kj = dt * b_kl/real(a_kj(state))

    call random_number(rnd)
    prob=0.d0
    do i=1,nquant
      if(i.ne.state.and.g_kj(i)>0.d0) then
        prob=prob+g_kj(i)
        if(rnd<prob) then
          call perform_hop(i)
          exit
        endif
      endif
    enddo

  end subroutine hop
  !-----------------------------------------------------------------  

  subroutine perform_hop(state_tentative)
    !! FOR 1 classical DIMENSION ONLY
    implicit none
    integer,intent(in) :: state_tentative
    real*8 apple

    apple=2/mass*(1/2.d0*mass*sum(v*v)+potential(state)-potential(state_tentative))
    if(apple>=0.d0) then
      state=state_tentative
      v=dsqrt(apple)*v/abs(v)
      call compute_acceleration
      call compute_energy
    endif

    if(apple<0.d0) then
      !! Frustrated hop

    endif

  end subroutine perform_hop
  !-----------------------------------------------------------------  

  function derivative(ci)
    implicit none
    complex*16,intent(in):: ci(nquant)
    complex*16 derivative(nquant)

    derivative=-iota*potential*ci/hbar - matmul(vdotd,ci)
    !print *, 'derivative=', derivative
    
 !print *, 'vdotd', vdotd
 !print *, 'ci', ci, nquant
 !print *, 'matmul(vdotd,ci)=', matmul(vdotd,ci)
   !print *, 'sum(matmul)=', sum(matmul(vdotd,ci))
  end function derivative
  !-----------------------------------------------------------------  

  subroutine compute_vdotd
    implicit none
    integer i,j

    call compute_dij
    do i=1,nquant
      do j=1,nquant
        vdotd(i,j)=sum(v*dij(i,j,:))
!	print *, 'vdotd(i,j)', i, j, vdotd(i,j)
      enddo
    enddo

  end subroutine compute_vdotd
  !-----------------------------------------------------------------  

  subroutine compute_dij
    implicit none
    integer i,j,k

    do i=1,nquant
      do j=1,nquant
        do k=1,nclass
          !! Look out for sign
          if(i.ne.j)dij(i,j,k)=sum(psi(:,i)*matmul(grad_H(:,:,k),psi(:,j)))/(potential(j)-potential(i))
!	print *, 'psi(:,i)=', psi(:,i)
!	print *, 'psi(:,j)=', psi(:,j)
!	print *, 'grad_H', grad_H(:,:,k), i,j,k
!	print *, 'matmul(grad_H,psi)=', matmul(grad_H(:,:,k),psi(:,j))
 !       print *, 'sum=', sum(psi(:,i)*matmul(grad_H(:,:,k),psi(:,j)))
          if(i==j)dij(i,j,k)=0.d0
        enddo
      enddo
    enddo

  end subroutine compute_dij
  !-----------------------------------------------------------------  

  subroutine compute_acceleration
    implicit none
    real*8 Hamil(nquant,nquant)
    integer i

    if(x(1)>0.d0) then
      Hamil(1,1)=AA*(1-exp(-BB*x(1)))
      grad_H(1,1,1)=AA*BB*exp(-BB*x(1))
    else
      Hamil(1,1)=-AA*(1-exp(BB*x(1)))
      grad_H(1,1,1)=AA*BB*exp(BB*x(1))
    endif

    Hamil(2,2)=-Hamil(1,1)
    grad_H(2,2,1)=-grad_H(1,1,1)

    Hamil(1,2)=CC*exp(-DD*x(1)**2)
    grad_H(1,2,1)=-2*CC*DD*x(1)*exp(-DD*x(1)**2)
    Hamil(2,1)=Hamil(1,2)
    grad_H(2,1,1)=grad_H(1,2,1)

    call diag(Hamil,psi,potential)

!    delv=0.d0
!    do k=1,nquant
!      do kp=1,nquant
!        do l=1,nclass
!          delV(l)=delv(l)+psi(k,state)*grad_H(k,kp,l)*psi(kp,state)
!        enddo
!      enddo
!    enddo

    do i=1,nclass
      acc(i)=-1.d0/mass * sum(psi(:,state)*matmul(grad_H(:,:,i),psi(:,state)))
    enddo

  end subroutine compute_acceleration
  !-----------------------------------------------------------------  

  subroutine diag(M,psi,potential)
    implicit none
    real*8,intent(in):: M(nquant,nquant)
    real*8,intent(out) :: psi(2,2),potential(2)

    potential(1)=0.5*(M(1,1)+M(2,2)-dsqrt((M(1,1)-M(2,2))**2+4*M(1,2)**2))
    potential(2)=0.5*(M(1,1)+M(2,2)+dsqrt((M(1,1)-M(2,2))**2+4*M(1,2)**2))

    psi(1,1)=-M(1,2)/dsqrt(M(1,2)**2+(M(1,1)-potential(1))**2)
    psi(2,1)=(M(1,1)-potential(1))/dsqrt(M(1,2)**2+(M(1,1)-potential(1))**2)
    psi(1,2)=-M(1,2)/dsqrt(M(1,2)**2+(M(1,1)-potential(2))**2)
    psi(2,2)=(M(1,1)-potential(2))/dsqrt(M(1,2)**2+(M(1,1)-potential(2))**2)
 
  end subroutine diag
  !-----------------------------------------------------------------  

  subroutine compute_energy
    implicit none
    energy=0.5*mass*sum(v*v) + potential(state)	!print *, 'v', v, 'v*v', v*v, 'sum(v*v)', sum(v*v)
  end subroutine compute_energy
  !-----------------------------------------------------------------  

  subroutine draw_pes
    implicit none
    integer i

   ! state=1
   ! x(1)=20.d0
   ! call compute_acceleration
   ! call compute_dij
   ! stop

    state=1
    do i=1,1000
      x(1)=-10.d0+i*30.d0/1000.d0
      call compute_acceleration
      call compute_dij
      write(20,*)x(1),potential
      write(21,*)x(1),dij(1,2,1)
      write(22,*) x(1), potential,dij(1,2,1)/50.d0
    enddo
    stop

  end subroutine draw_pes
  !-----------------------------------------------------------------  

  subroutine check_acceleration
    implicit none
    integer i
    real*8 delx
    real*8 pot_sav,acc_sav(nclass),acc(nclass)

    state=1
    x(1)=1.d0
    call compute_acceleration
    pot_sav=potential(1)
    acc_sav=acc

    delx=1.d-6
    x(1)=x(1)+delx
    call compute_acceleration
    acc=-1.d0/mass*((potential(1)-pot_sav)/delx)
print *, acc, acc_sav
  end subroutine check_acceleration
  !-----------------------------------------------------------------  

end program md
