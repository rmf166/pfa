    module global

      implicit none

      integer(4),    parameter     ::                    &
          sp = kind(1.0),                                &
          dp = selected_real_kind(2*precision(1.0_sp)),  &
          qp = selected_real_kind(2*precision(1.0_dp))
      integer(4),    parameter     :: kr=qp
      integer(4)                   :: sweeps
      logical(4)                   :: linpro

    end module global

    program main

      implicit none

      !$ call omp_set_num_threads(4)

      call drive_fa

      contains

      subroutine drive_fa

        use global
        !$ use omp_lib, only: omp_set_num_threads

        implicit none

        integer(4)                   :: i
        integer(4)                   :: k
        integer(4)                   :: sol
        integer(4)                   :: p
        integer(4)                   :: s
        integer(4)                   :: xn
        integer(4)                   :: jmax
        integer(4),    parameter     :: n=16
        integer(4),    parameter     :: kmax=1000
        real(kind=kr)                :: x
        real(kind=kr)                :: tau
        real(kind=kr)                :: c(1:4)=(/0.8000_kr,0.9000_kr,0.9900_kr,0.9999_kr/)
        real(kind=kr)                :: h(19)
        real(kind=kr)                :: rho(19,4,6,4,2)
        character(1)                 :: cm
        character(1)                 :: sw
        character(2)                 :: cols
        character(2)                 :: solopt(0:3)=(/'DD','SC','LD','LC'/)
        character(10)                :: caserun
        character(20)                :: datafile
        character(132)               :: datfmt

        write(cols,'(i2)') 1+4
        write(datfmt,'(a)') '(' // trim(adjustl(cols)) // '(es12.5))'

        !$ call omp_set_num_threads(n/2)

        h  =0.0_kr
        rho=0.0_kr
        do k=1,2
          linpro=.false.
          if (k == 2) linpro=.true.
          do sol=0,3
            do i=1,4
              sweeps=i
              do p=1,6
                do s=1,4
                  x=50.0_kr
                  tau=100.0_kr*(0.6_kr**18)
                  do xn=1,19
                    if (0.1_kr <= tau .and. tau < 1.0_kr) then
                      x=50.0_kr
                    elseif (1.0_kr <= tau .and. tau < 10.0_kr) then
                      x=500.0_kr
                    elseif (10.0_kr <= tau .and. tau < 500.0_kr) then
                      x=5000.0_kr
                    endif
                    h(xn)=tau
                    jmax=x/h(xn)
                    jmax=(jmax/p)*p
                    h(xn)=x/jmax
                    if (mod(jmax,p) == 1) stop 'Incorrect JMAX value'
                    call solve_slab_fa(sol,p,c(s),n,kmax,jmax,h(xn),rho(xn,s,p,i,k))
                    tau=tau*(1.0_kr/0.6_kr)
                  enddo
                enddo
              enddo
            enddo
            do i=1,4
              do p=1,6
                write(cm,'(i1)') p
                write(sw,'(i1)') i
                if (linpro) then
                  write(caserun,'(a)') '-lp' // trim(adjustl(cm)) // '-s' // trim(adjustl(sw)) // '-' // solopt(sol)
                else
                  write(caserun,'(a)') '-p' // trim(adjustl(cm)) // '-s' // trim(adjustl(sw)) // '-' // solopt(sol)
                endif
                write(datafile,'(a)') 'numres' // trim(adjustl(caserun)) // '.dat'
                open(unit=1,file=datafile,action='write',status='unknown')
                do xn=1,19
                  write(1,datfmt) h(xn),(rho(xn,s,p,i,k),s=1,4)
                enddo
                close(1)
              enddo
            enddo
          enddo
        enddo

      end subroutine drive_fa

      subroutine solve_slab_fa(sol,p,c,n,kmax,jmax,h,rho)

        use global

        implicit none

        integer(4),    intent(in)    :: sol
        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: p
        integer(4),    intent(in)    :: jmax
        real(kind=kr), intent(in)    :: c
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(out)   :: rho

        integer(4)                   :: bc(2)
        real(kind=kr)                :: eps
        real(kind=kr)                :: mu(n/2)
        real(kind=kr)                :: w (n/2)
        real(kind=kr), allocatable   :: phi (:)
        real(kind=kr), allocatable   :: jnet(:)
        real(kind=kr), allocatable   :: jp(:)
        real(kind=kr), allocatable   :: jm(:)
        real(kind=kr)                :: q
        real(kind=kr)                :: sigt
        real(kind=kr)                :: sigs

      ! dynamic allocation of arrays

        allocate(phi(jmax))
        allocate(jnet(jmax+1))
        allocate(jp(jmax+1))
        allocate(jm(jmax+1))
        phi=0.0_kr
        jnet=0.0_kr
        jp=0.0_kr
        jm=0.0_kr

      ! build source based on options

        q=1.0_kr
        bc(1)=0
        bc(2)=1

      ! set cross sections

        sigt=1.0_kr
        sigs=c*sigt

      ! set quadrature

        call quad(n,mu,w)

      ! solve fixed-source problem

        eps=1.0e-06
        if (sol == 0) then
          call solve_dd(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet,jp,jm,p,rho)
        elseif (sol == 1) then
          call solve_sc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet,jp,jm,p,rho)
        elseif (sol == 2) then
          call solve_ld(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet,jp,jm,p,rho)
        elseif (sol == 3) then
          call solve_lc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet,jp,jm,p,rho)
        else
          write(0,'(a)') ' Incorrect solution scheme selected.'
          stop
        endif

      ! clean up arrays

        deallocate(phi)
        deallocate(jnet)
        deallocate(jp)
        deallocate(jm)

      end subroutine solve_slab_fa

      subroutine quad(n,mu,w)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        real(kind=kr), intent(out)   :: mu(n/2)
        real(kind=kr), intent(out)   :: w (n/2)

        integer(4)                   :: j
        integer(4),    parameter     :: nmaxp=300
        real(kind=kr)                :: xnew(nmaxp)
        real(kind=kr)                :: wnew(nmaxp)

        xnew=0.0_kr
        wnew=0.0_kr
        call gauleg(-1.0_kr,1.0_kr,xnew,wnew,n)

        do j=1,n/2
          mu(j)=xnew(j)
          w(j) =wnew(j)
        enddo

      end subroutine quad

      subroutine gauleg(x1,x2,x,w,n)
 
        use global

        implicit none

        integer(4),    intent(in)    :: n
        real(kind=kr), intent(in)    :: x1
        real(kind=kr), intent(in)    :: x2
        real(kind=kr), intent(inout) :: x(n)
        real(kind=kr), intent(inout) :: w(n)

        integer(4)                   :: i
        integer(4)                   :: j
        integer(4)                   :: m
        integer(4)                   :: kount
        integer(4),    parameter     :: nmax=300
        real(kind=kr)                :: xm
        real(kind=kr)                :: xl
        real(kind=kr)                :: p1
        real(kind=kr)                :: p2
        real(kind=kr)                :: p3
        real(kind=kr)                :: pi
        real(kind=kr)                :: pp
        real(kind=kr)                :: z
        real(kind=kr)                :: z1
        real(kind=kr)                :: xtmp(nmax)  ! full set of abscissas
        real(kind=kr)                :: wtmp(nmax)  ! full set of weights
        real(kind=kr), parameter     :: eps=1.0e-30

        pi=4.0_kr*atan(1.0_kr)
        if (n > nmax) then
          write(0,'(a,1i6)') 'Gauss-Leg. integration problem --Increase PARAMETER: NMAX to at least:',n
          stop
        endif

        m=(n+1)/2
        xm=0.5_kr*(x2+x1)
        xl=0.5_kr*(x2-x1)
        do i=1,m
          z=cos(pi*(i-0.25_kr)/(n+0.5_kr))
      1   continue
          p1=1.0_kr
          p2=0.0_kr
          do j=1,n
            p3=p2
            p2=p1
            p1=((2.0_kr*j-1.0_kr)*z*p2-(j-1.0_kr)*p3)/j
          enddo
      !   p1 is now the desired Legendre polynomial. we next compute pp, its derivative,
      !   by a standard relation involving also p2, the polynomial of one lower order.
          pp=n*(z*p1-p2)/(z*z-1.0_kr)
          z1=z
          z=z1-p1/pp
          if (abs(z-z1) > eps) go to 1
          xtmp(i)=    xm-xl*z
          xtmp(n+1-i)=xm+xl*z
      !   the (n+1-i) terms are the symmetric counterparts
          wtmp(i)=2.0_kr*xl/((1.0_kr-z*z)*pp*pp)
          wtmp(n+1-i)=wtmp(i)
        enddo

      ! (half set and assumed symmetric)
        kount=0
        do i=1,n
          if (xtmp(i) >= 0.0_kr) then
            kount=kount+1
            x(kount)=xtmp(i)   ! abscissas
            w(kount)=wtmp(i)   ! weights
          endif
        enddo

      end subroutine gauleg

      subroutine solve_dd(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet,jp,jm,p,rho)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: bc(2)
        integer(4),    intent(in)    :: p
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q
        real(kind=kr), intent(in)    :: eps
        real(kind=kr), intent(in)    :: sigt
        real(kind=kr), intent(in)    :: sigs
        real(kind=kr), intent(in)    :: mu(n/2)
        real(kind=kr), intent(in)    :: w (n/2)
        real(kind=kr), intent(inout) :: phi (jmax)
        real(kind=kr), intent(inout) :: jnet(jmax+1)
        real(kind=kr), intent(inout) :: jp(jmax+1)
        real(kind=kr), intent(inout) :: jm(jmax+1)
        real(kind=kr), intent(inout) :: rho

        integer(4)                   :: sw
        integer(4)                   :: j
        integer(4)                   :: k
        integer(4)                   :: kount
        integer(4)                   :: m
        real(kind=kr)                :: tau
        real(kind=kr)                :: norm0
        real(kind=kr)                :: norm1
        real(kind=kr)                :: psi
        real(kind=kr)                :: psi_in
        real(kind=kr)                :: psi_bc(n/2)
        real(kind=kr), allocatable   :: c1(:,:)
        real(kind=kr), allocatable   :: c2(:,:)
        real(kind=kr), allocatable   :: phio(:)
        real(kind=kr), allocatable   :: phil(:)
        real(kind=kr), allocatable   :: s(:)
        real(kind=kr), allocatable   :: rhoi(:)

      ! pre-compute coeffs

        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        c1=0.0_kr
        c2=0.0_kr

        do m=1,n/2
          do j=1,jmax
            tau=sigt*h/mu(m)
            c1(j,m)=0.5_kr*tau
            c2(j,m)=0.5_kr*h/mu(m)
          enddo
        enddo

      ! solve problem

        allocate(phil(jmax))
        allocate(s(jmax))
        allocate(rhoi(kmax))
        phi =1.0_kr
        phil=0.0_kr
        rhoi=0.0_kr

        psi_in=0.0_kr
        psi_bc=0.0_kr
        kount=0

        allocate(phio(jmax))
        do k=1,kmax
          phio=phi
          do sw=1,sweeps
            do j=1,jmax
              s  (j)=0.5_kr*(sigs*phi(j)+q)
            enddo
            phi=0.0_kr
            jnet=0.0_kr
            jp=0.0_kr
            jm=0.0_kr
            !$omp parallel do private(j,psi_in,psi) reduction(+:phi,jp,jm)
            do m=1,n/2
              psi_in=psi_bc(m) ! left specular bc
              if (bc(1) == 0) psi_in=0.0_kr
              do j=1,jmax
                jp(j)   =jp(j)+psi_in*mu(m)*w(m)
                psi     =(s(j)*c2(j,m)+psi_in)/(1.0_kr+c1(j,m))
                phi(j)  =phi(j)+psi*w(m)
                psi_in  =2.0_kr*psi-psi_in
              enddo
              jp(jmax+1)=jp(jmax+1)+psi_in*mu(m)*w(m)
              if (bc(2) == 0) psi_in=0.0_kr
              jm(jmax+1)=jm(jmax+1)+psi_in*mu(m)*w(m)
              do j=jmax,1,-1
                psi     =(s(j)*c2(j,m)+psi_in)/(1.0_kr+c1(j,m))
                phi(j)  =phi(j)+psi*w(m)
                psi_in  =2.0_kr*psi-psi_in
                jm(j)   =jm(j)+psi_in*mu(m)*w(m)
              enddo
              psi_bc(m)=psi_in
            enddo
            !$omp end parallel do
            jnet=jp-jm
          enddo

          call pcmfd(sigt,sigs,h,p,jmax,bc,q,phi,phil,jp,jm)

          norm1=0.0_kr
          do j=1,jmax
            if (phi(j) > 1.0e+33 .or. phi(j) < 0.0_kr) then
              rho=10.0_kr
              return
            endif
            norm1=norm1+(phi(j)-phio(j))**2
          enddo
          norm1=sqrt(norm1)

          if (norm1 < 0.001_kr) then
            kount=kount+1
            rhoi(kount)=norm1/norm0
          endif
          norm0=norm1

          if (norm1 <= eps) exit 
          if (k > 100 .and. norm1 > 100.0_kr) exit

        enddo

        if (norm1 > eps) then
          rho=10.0_kr
        else
          rho=0.0_kr
          do j=1,kount
            rho=rho+rhoi(j)
          enddo
          rho=rho/kount
        endif

        deallocate(c1)
        deallocate(c2)
        deallocate(s)
        deallocate(rhoi)
        deallocate(phio)
        deallocate(phil)

      end subroutine solve_dd

      subroutine solve_sc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet,jp,jm,p,rho)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: bc(2)
        integer(4),    intent(in)    :: p
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q
        real(kind=kr), intent(in)    :: eps
        real(kind=kr), intent(in)    :: sigt
        real(kind=kr), intent(in)    :: sigs
        real(kind=kr), intent(in)    :: mu(n/2)
        real(kind=kr), intent(in)    :: w (n/2)
        real(kind=kr), intent(inout) :: phi (jmax)
        real(kind=kr), intent(inout) :: jnet(jmax+1)
        real(kind=kr), intent(inout) :: jp(jmax+1)
        real(kind=kr), intent(inout) :: jm(jmax+1)
        real(kind=kr), intent(inout) :: rho

        integer(4)                   :: sw
        integer(4)                   :: j
        integer(4)                   :: k
        integer(4)                   :: kount
        integer(4)                   :: m
        real(kind=kr)                :: tau
        real(kind=kr)                :: tau3
        real(kind=kr)                :: tau5
        real(kind=kr)                :: tau7
        real(kind=kr)                :: norm0
        real(kind=kr)                :: norm1
        real(kind=kr)                :: psi
        real(kind=kr)                :: psi_in
        real(kind=kr)                :: psi_bc(n/2)
        real(kind=kr), allocatable   :: alpha(:,:)
        real(kind=kr), allocatable   :: c1(:,:)
        real(kind=kr), allocatable   :: c2(:,:)
        real(kind=kr), allocatable   :: phio(:)
        real(kind=kr), allocatable   :: phil(:)
        real(kind=kr), allocatable   :: s(:)
        real(kind=kr), allocatable   :: rhoi(:)

      ! pre-compute coeffs

        allocate(alpha(jmax,n/2))
        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        alpha=0.0_kr
        c1   =0.0_kr
        c2   =0.0_kr

        do m=1,n/2
          do j=1,jmax
            tau=sigt*h/mu(m)
            if (tau < 0.01_kr) then
              tau3=tau *tau*tau
              tau5=tau3*tau*tau
              tau7=tau5*tau*tau
              alpha(j,m)=tau/6.0_kr-tau3/360.0_kr+tau5/15120.0_kr-tau7/604800.0_kr
            else
              alpha(j,m)=1.0_kr/tanh(tau/2.0_kr)-2.0_kr/tau
            endif
            c1(j,m)=0.5_kr*tau*(1.0_kr+alpha(j,m))
            c2(j,m)=0.5_kr*h  *(1.0_kr+alpha(j,m))/mu(m)
          enddo
        enddo

      ! solve problem

        allocate(phil(jmax))
        allocate(s(jmax))
        allocate(rhoi(kmax))
        phi =1.0_kr
        phil=0.0_kr
        rhoi=0.0_kr

        psi_in=0.0_kr
        psi_bc=0.0_kr
        kount=0

        allocate(phio(jmax))
        do k=1,kmax
          phio=phi
          do sw=1,sweeps
            do j=1,jmax
              s  (j)=0.5_kr*(sigs*phi  (j)+q)
            enddo
            phi=0.0_kr
            jnet=0.0_kr
            jp=0.0_kr
            jm=0.0_kr
            !$omp parallel do private(j,psi_in,psi) reduction(+:phi,jp,jm)
            do m=1,n/2
              psi_in=psi_bc(m) ! left specular bc
              if (bc(1) == 0) psi_in=0.0_kr
              do j=1,jmax
                jp(j)   =jp(j)+psi_in*mu(m)*w(m)
                psi     =(s(j)*c2(j,m)+psi_in)/(1.0_kr+c1(j,m))
                phi(j)  =phi(j)+psi*w(m)
                psi_in  =(2.0_kr*psi-(1.0_kr-alpha(j,m))*psi_in)/(1.0_kr+alpha(j,m))
              enddo
              jp(jmax+1)=jp(jmax+1)+psi_in*mu(m)*w(m)
              if (bc(2) == 0) psi_in=0.0_kr
              jm(jmax+1)=jm(jmax+1)+psi_in*mu(m)*w(m)
              do j=jmax,1,-1
                psi     =(s(j)*c2(j,m)+psi_in)/(1.0_kr+c1(j,m))
                phi(j)  =phi(j)+psi*w(m)
                psi_in  =(2.0_kr*psi-(1.0_kr-alpha(j,m))*psi_in)/(1.0_kr+alpha(j,m))
                jm(j)   =jm(j)+psi_in*mu(m)*w(m)
              enddo
              psi_bc(m)=psi_in
            enddo
            !$omp end parallel do
            jnet=jp-jm
          enddo

          call pcmfd(sigt,sigs,h,p,jmax,bc,q,phi,phil,jp,jm)

          norm1=0.0_kr
          do j=1,jmax
            if (phi(j) > 1.0e+33 .or. phi(j) < 0.0_kr) then
              rho=10.0_kr
              return
            endif
            norm1=norm1+(phi(j)-phio(j))**2
          enddo
          norm1=sqrt(norm1)

          if (norm1 < 0.001_kr) then
            kount=kount+1
            rhoi(kount)=norm1/norm0
          endif
          norm0=norm1

          if (norm1 <= eps) exit 
          if (k > 100 .and. norm1 > 100.0_kr) exit

        enddo

        if (norm1 > eps) then
          rho=10.0_kr
        else
          rho=0.0_kr
          do j=1,kount
            rho=rho+rhoi(j)
          enddo
          rho=rho/kount
        endif

        deallocate(alpha)
        deallocate(c1)
        deallocate(c2)
        deallocate(phio)
        deallocate(phil)
        deallocate(s)
        deallocate(rhoi)

      end subroutine solve_sc

      subroutine solve_ld(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet,jp,jm,p,rho)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: bc(2)
        integer(4),    intent(in)    :: p
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q
        real(kind=kr), intent(in)    :: eps
        real(kind=kr), intent(in)    :: sigt
        real(kind=kr), intent(in)    :: sigs
        real(kind=kr), intent(in)    :: mu(n/2)
        real(kind=kr), intent(in)    :: w (n/2)
        real(kind=kr), intent(inout) :: phi (jmax)
        real(kind=kr), intent(inout) :: jnet(jmax+1)
        real(kind=kr), intent(inout) :: jp(jmax+1)
        real(kind=kr), intent(inout) :: jm(jmax+1)
        real(kind=kr), intent(inout) :: rho

        integer(4)                   :: sw
        integer(4)                   :: j
        integer(4)                   :: k
        integer(4)                   :: kount
        integer(4)                   :: m
        real(kind=kr)                :: tau
        real(kind=kr)                :: norm0
        real(kind=kr)                :: norm1
        real(kind=kr)                :: psi
        real(kind=kr)                :: psil
        real(kind=kr)                :: psi_in
        real(kind=kr)                :: psi_out
        real(kind=kr)                :: psi_bc(n/2)
        real(kind=kr), allocatable   :: alpha(:,:)
        real(kind=kr), allocatable   :: c1(:,:)
        real(kind=kr), allocatable   :: c2(:,:)
        real(kind=kr), allocatable   :: phio(:)
        real(kind=kr), allocatable   :: phil(:)
        real(kind=kr), allocatable   :: philo(:)
        real(kind=kr), allocatable   :: s(:)
        real(kind=kr), allocatable   :: sl(:)
        real(kind=kr), allocatable   :: rhoi(:)

      ! pre-compute coeffs

        allocate(alpha(jmax,n/2))
        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        alpha=0.0_kr
        c1   =0.0_kr
        c2   =0.0_kr

        do m=1,n/2
          do j=1,jmax
            tau=sigt*h/mu(m)
            alpha(j,m)=1.0_kr/(1.0_kr+6.0_kr/tau)
            c1(j,m)   =       (2.0_kr/tau+alpha(j,m)-1.0_kr)
            c2(j,m)   =1.0_kr/(2.0_kr/tau+alpha(j,m)+1.0_kr)
          enddo
        enddo

      ! solve problem

        allocate(phil(jmax))
        allocate(s(jmax))
        allocate(sl(jmax))
        allocate(rhoi(kmax))
        phi =1.0_kr
        phil=0.0_kr
        rhoi=0.0_kr

        psi_in=0.0_kr
        psi_bc=0.0_kr
        kount=0

        allocate(phio(jmax))
        allocate(philo(jmax))
        do k=1,kmax
          phio=phi
          philo=phil
          do sw=1,sweeps
            do j=1,jmax
              s  (j)=0.5_kr*(sigs*phi  (j)+q)
              sl (j)=0.5_kr*(sigs*phil (j))
            enddo
            phi=0.0_kr
            phil=0.0_kr
            jnet=0.0_kr
            jp=0.0_kr
            jm=0.0_kr
            !$omp parallel do private(j,psi_in,psi_out,psi,psil) reduction(+:phi,phil,jp,jm)
            do m=1,n/2
              psi_in=psi_bc(m) ! left specular bc
              if (bc(1) == 0) psi_in=0.0_kr
              do j=1,jmax
                jp(j)  =jp(j)+psi_in*mu(m)*w(m)
                psi_out=c2(j,m)*(2.0_kr*(s(j)+alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
                psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr-alpha(j,m)*sl(j)/sigt
                psil   =psi_out-psi
                psi_in =psi_out
                phi(j) =phi(j)+psi*w(m)
                phil(j)=phil(j)+psil*w(m)
              enddo
              jp(jmax+1)=jp(jmax+1)+psi_in*mu(m)*w(m)
              if (bc(2) == 0) psi_in=0.0_kr
              jm(jmax+1)=jm(jmax+1)+psi_in*mu(m)*w(m)
              do j=jmax,1,-1
                psi_out=c2(j,m)*(2.0_kr*(s(j)-alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
                psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr+alpha(j,m)*sl(j)/sigt
                psil   =psi-psi_out
                psi_in =psi_out
                phi(j) =phi(j)+psi*w(m)
                phil(j)=phil(j)+psil*w(m)
                jm(j)  =jm(j)+psi_in*mu(m)*w(m)
              enddo
              psi_bc(m)=psi_in
            enddo
            !$omp end parallel do
            jnet=jp-jm
          enddo

          call pcmfd(sigt,sigs,h,p,jmax,bc,q,phi,phil,jp,jm)

          norm1=0.0_kr
          do j=1,jmax
            if (phi(j) > 1.0e+33 .or. phi(j) < 0.0_kr) then
              rho=10.0_kr
              return
            endif
            norm1=norm1+(phi(j)-phio(j))**2 ! +(phil(j)-philo(j))**2 is close to zero
          enddo
          norm1=sqrt(norm1)

          if (norm1 < 0.001_kr) then
            kount=kount+1
            rhoi(kount)=norm1/norm0
          endif
          norm0=norm1

          if (norm1 <= eps) exit 
          if (k > 100 .and. norm1 > 100.0_kr) exit

        enddo

        if (norm1 > eps) then
          rho=10.0_kr
        else
          rho=0.0_kr
          do j=1,kount
            rho=rho+rhoi(j)
          enddo
          rho=rho/kount
        endif

        deallocate(alpha)
        deallocate(c1)
        deallocate(c2)
        deallocate(phio)
        deallocate(phil)
        deallocate(philo)
        deallocate(s)
        deallocate(sl)
        deallocate(rhoi)

      end subroutine solve_ld

      subroutine solve_lc(n,jmax,kmax,h,q,eps,sigt,sigs,mu,w,bc,phi,jnet,jp,jm,p,rho)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: kmax
        integer(4),    intent(in)    :: bc(2)
        integer(4),    intent(in)    :: p
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q
        real(kind=kr), intent(in)    :: eps
        real(kind=kr), intent(in)    :: sigt
        real(kind=kr), intent(in)    :: sigs
        real(kind=kr), intent(in)    :: mu(n/2)
        real(kind=kr), intent(in)    :: w (n/2)
        real(kind=kr), intent(inout) :: phi (jmax)
        real(kind=kr), intent(inout) :: jnet(jmax+1)
        real(kind=kr), intent(inout) :: jp(jmax+1)
        real(kind=kr), intent(inout) :: jm(jmax+1)
        real(kind=kr), intent(inout) :: rho

        integer(4)                   :: sw
        integer(4)                   :: j
        integer(4)                   :: k
        integer(4)                   :: kount
        integer(4)                   :: m
        real(kind=kr)                :: tau
        real(kind=kr)                :: tau3
        real(kind=kr)                :: tau5
        real(kind=kr)                :: tau7
        real(kind=kr)                :: norm0
        real(kind=kr)                :: norm1
        real(kind=kr)                :: psi
        real(kind=kr)                :: psil
        real(kind=kr)                :: psi_in
        real(kind=kr)                :: psi_out
        real(kind=kr)                :: psi_bc(n/2)
        real(kind=kr), allocatable   :: alpha(:,:)
        real(kind=kr), allocatable   :: rbeta(:,:)
        real(kind=kr), allocatable   :: c1(:,:)
        real(kind=kr), allocatable   :: c2(:,:)
        real(kind=kr), allocatable   :: phio(:)
        real(kind=kr), allocatable   :: phil(:)
        real(kind=kr), allocatable   :: philo(:)
        real(kind=kr), allocatable   :: s(:)
        real(kind=kr), allocatable   :: sl(:)
        real(kind=kr), allocatable   :: rhoi(:)

      ! pre-compute coeffs

        allocate(alpha(jmax,n/2))
        allocate(rbeta(jmax,n/2))
        allocate(c1(jmax,n/2))
        allocate(c2(jmax,n/2))
        alpha=0.0_kr
        rbeta=0.0_kr
        c1   =0.0_kr
        c2   =0.0_kr

        do m=1,n/2
          do j=1,jmax
            tau=sigt*h/mu(m)
            if (tau < 0.01_kr) then
              tau3=tau *tau*tau
              tau5=tau3*tau*tau
              tau7=tau5*tau*tau
              alpha(j,m)=tau/6.0_kr-tau3/360.0_kr+tau5/15120.0_kr-tau7/604800.0_kr
            else
              alpha(j,m)=1.0_kr/tanh(tau/2.0_kr)-2.0_kr/tau
            endif
            rbeta(j,m)=1.0_kr/alpha(j,m)-6.0_kr/tau
            c1(j,m)=       (2.0_kr/tau+alpha(j,m)-1.0_kr)
            c2(j,m)=1.0_kr/(2.0_kr/tau+alpha(j,m)+1.0_kr)
          enddo
        enddo

      ! solve problem

        allocate(phil(jmax))
        allocate(s(jmax))
        allocate(sl(jmax))
        allocate(rhoi(kmax))
        phi =1.0_kr
        phil=0.0_kr
        rhoi=0.0_kr

        psi_in=0.0_kr
        psi_bc=0.0_kr
        kount=0

        allocate(phio(jmax))
        allocate(philo(jmax))
        do k=1,kmax
          phio=phi
          philo=phil
          do sw=1,sweeps
            do j=1,jmax
              s  (j)=0.5_kr*(sigs*phi  (j)+q)
              sl (j)=0.5_kr*(sigs*phil (j))
            enddo
            phi=0.0_kr
            phil=0.0_kr
            jnet=0.0_kr
            jp=0.0_kr
            jm=0.0_kr
            !$omp parallel do private(j,psi_in,psi_out,psi,psil) reduction(+:phi,phil,jp,jm)
            do m=1,n/2
              psi_in=psi_bc(m) ! left specular bc
              if (bc(1) == 0) psi_in=0.0_kr
              do j=1,jmax
                jp(j)  =jp(j)+psi_in*mu(m)*w(m)
                psi_out=c2(j,m)*(2.0_kr*(s(j)+alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
                psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr-alpha(j,m)*sl(j)/sigt
                psil   =((rbeta(j,m)+1.0_kr)*psi_out+(rbeta(j,m)-1.0_kr)*psi_in)/2.0_kr-rbeta(j,m)*psi
                psi_in =psi_out
                phi(j) =phi(j)+psi*w(m)
                phil(j)=phil(j)+psil*w(m)
              enddo
              jp(jmax+1)=jp(jmax+1)+psi_in*mu(m)*w(m)
              if (bc(2) == 0) psi_in=0.0_kr
              jm(jmax+1)=jm(jmax+1)+psi_in*mu(m)*w(m)
              do j=jmax,1,-1
                psi_out=c2(j,m)*(2.0_kr*(s(j)-alpha(j,m)*sl(j))/sigt+c1(j,m)*psi_in)
                psi    =((1.0_kr+alpha(j,m))*psi_out+(1.0_kr-alpha(j,m))*psi_in)/2.0_kr+alpha(j,m)*sl(j)/sigt
                psil   =-((rbeta(j,m)+1.0_kr)*psi_out+(rbeta(j,m)-1.0_kr)*psi_in)/2.0_kr+rbeta(j,m)*psi
                psi_in =psi_out
                phi(j) =phi(j)+psi*w(m)
                phil(j)=phil(j)+psil*w(m)
                jm(j)  =jm(j)+psi_in*mu(m)*w(m)
              enddo
              psi_bc(m)=psi_in
            enddo
            !$omp end parallel do
            jnet=jp-jm
          enddo

          call pcmfd(sigt,sigs,h,p,jmax,bc,q,phi,phil,jp,jm)

          norm1=0.0_kr
          do j=1,jmax
            if (phi(j) > 1.0e+33 .or. phi(j) < 0.0_kr) then
              rho=10.0_kr
              return
            endif
            norm1=norm1+(phi(j)-phio(j))**2 ! +(phil(j)-philo(j))**2 is close to zero
          enddo
          norm1=sqrt(norm1)

          if (norm1 < 0.001_kr) then
            kount=kount+1
            rhoi(kount)=norm1/norm0
          endif
          norm0=norm1

          if (norm1 <= eps) exit 
          if (k > 100 .and. norm1 > 100.0_kr) exit

        enddo

        if (norm1 > eps) then
          rho=10.0_kr
        else
          rho=0.0_kr
          do j=1,kount
            rho=rho+rhoi(j)
          enddo
          rho=rho/kount
        endif

        deallocate(alpha)
        deallocate(rbeta)
        deallocate(c1)
        deallocate(c2)
        deallocate(phio)
        deallocate(phil)
        deallocate(philo)
        deallocate(s)
        deallocate(sl)
        deallocate(rhoi)

      end subroutine solve_lc

      subroutine pcmfd(sigt,sigs,h,p,jmax,bc,q,phi,phil,jp,jm)

        use global

        implicit none

        integer(4),    intent(in)    :: p
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: bc(2)
        real(kind=kr), intent(in)    :: sigt
        real(kind=kr), intent(in)    :: sigs
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(in)    :: q
        real(kind=kr), intent(inout) :: phi (jmax)
        real(kind=kr), intent(inout) :: phil(jmax)
        real(kind=kr), intent(in)    :: jp(jmax+1)
        real(kind=kr), intent(in)    :: jm(jmax+1)

        integer(4)                   :: j
        integer(4)                   :: jj
        integer(4)                   :: n
        integer(4)                   :: nleft
        integer(4)                   :: nrite
        integer(4)                   :: nmax
        real(kind=kr)                :: dphi
        real(kind=kr)                :: dphil
        real(kind=kr)                :: dphir
        real(kind=kr)                :: dc
        real(kind=kr)                :: sumphi
        real(kind=kr)                :: siga
        real(kind=kr)                :: xj
        real(kind=kr), allocatable   :: phin(:)
        real(kind=kr), allocatable   :: phio(:)
        real(kind=kr), allocatable   :: jnp(:)
        real(kind=kr), allocatable   :: jnm(:)
        real(kind=kr), allocatable   :: s(:)
        real(kind=kr), allocatable   :: dd(:)
        real(kind=kr), allocatable   :: ccp(:)
        real(kind=kr), allocatable   :: ccm(:)
        real(kind=kr), allocatable   :: a(:)
        real(kind=kr), allocatable   :: b(:)
        real(kind=kr), allocatable   :: c(:)

        if (mod(jmax,p) /= 0) stop ' Fine mesh does not align with pCMFD.'
        if (bc(1) /= 0) stop ' Only vacuum BC supported on left edge.'
        if (bc(2) /= 1) stop ' Only reflective BC supported on right edge.'

        nmax=jmax/p

        allocate(phin(nmax))
        allocate(phio(nmax))
        allocate(jnp(nmax+1))
        allocate(jnm(nmax+1))
        allocate(s(nmax))
        allocate(dd(nmax+1))
        allocate(ccp(nmax+1))
        allocate(ccm(nmax+1))
        allocate(a(nmax))
        allocate(b(nmax))
        allocate(c(nmax))
        phin=0.0_kr
        jnp =0.0_kr
        jnm =0.0_kr
        s   =0.0_kr
        dd  =0.0_kr
        ccp =0.0_kr
        ccm =0.0_kr
        a   =0.0_kr
        b   =0.0_kr
        c   =0.0_kr

        j=1
        do n=1,nmax
          jnp(n)=jp(j)
          jnm(n)=jm(j)
          sumphi=0.0_kr
          do jj=1,p
            sumphi=sumphi+phi(j)
            j=j+1
          enddo
          phin(n)=sumphi/p
          s(n)   =p*h*q
        enddo
        jnp(nmax+1)=jp(jmax+1)
        jnm(nmax+1)=jm(jmax+1)

        phio  = phin
        dc    = 1.0_kr/(3.0_kr*p*h*sigt)
        dd    = dc
        dd(1) = dc/(0.5_kr+2.0_kr*dc)
        n     = nmax+1
        dd(n) = 0.0_kr
        ccm(1)=(2.0_kr*jnm(1)-dd(1)*phio(1))/(2.0_kr*phio(1))
        do n=2,nmax
          ccp(n)=(2.0_kr*jnp(n)+dd(n)*(phio(n)-phio(n-1)))/(2.0_kr*phio(n-1))
          ccm(n)=(2.0_kr*jnm(n)-dd(n)*(phio(n)-phio(n-1)))/(2.0_kr*phio(n))
        enddo
        siga=sigt-sigs
        b(1)=dd(2)+ccp(2)+0.5_kr*dd(1)+ccm(1)+p*h*siga
        c(1)=-(dd(2)+ccm(2))
        do n=2,nmax-1
          a(n)=-(dd(n)+ccp(n))
          b(n)=dd(n+1)+ccp(n+1)+dd(n)+ccm(n)+p*h*siga
          c(n)=-(dd(n+1)+ccm(n+1))
        enddo
        n   =nmax
        a(n)=-(dd(n)+ccp(n))
        b(n)=dd(n)+ccm(n)+p*h*siga

        phin=s
        call tdma(nmax,a,b,c,phin)

        if (linpro) then
          j=1
          do n=1,nmax
            nleft=n-1
            nrite=n+1
            if (n == 1) nleft=n
            if (n == nmax) nrite=n
            dphil=0.5_kr*((phin(nleft)/phio(nleft))+(phin(n)/phio(n)))
            dphir=0.5_kr*((phin(nrite)/phio(nrite))+(phin(n)/phio(n)))
            do jj=1,p
              xj=h/2.0_kr+(jj-1.0_kr)*h
              dphi=dphil+(xj/(p*h))*(dphir-dphil)
              phi(j) =phi(j) *dphi
              phil(j)=phil(j)*dphi
              j=j+1
            enddo
          enddo
        else
          j=1
          do n=1,nmax
            do jj=1,p
              phi(j) =(phin(n)/phio(n))*phi(j)
              phil(j)=(phin(n)/phio(n))*phil(j)
              j=j+1
            enddo
          enddo
        endif

        deallocate(phin)
        deallocate(phio)
        deallocate(jnp)
        deallocate(jnm)
        deallocate(s)
        deallocate(dd)
        deallocate(ccp)
        deallocate(ccm)
        deallocate(a)
        deallocate(b)
        deallocate(c)

      end subroutine pcmfd

      subroutine cmdsa(sigt,sigs,h,p,jmax,bc,phi,phi0)

        use global

        implicit none

        integer(4),    intent(in)    :: p
        integer(4),    intent(in)    :: jmax
        integer(4),    intent(in)    :: bc(2)
        real(kind=kr), intent(in)    :: sigt
        real(kind=kr), intent(in)    :: sigs
        real(kind=kr), intent(in)    :: h
        real(kind=kr), intent(inout) :: phi (jmax)
        real(kind=kr), intent(in)    :: phi0(jmax)

        integer(4)                   :: j
        integer(4)                   :: jj
        integer(4)                   :: n
        integer(4)                   :: nmax
        real(kind=kr)                :: dc
        real(kind=kr)                :: sumphi
        real(kind=kr)                :: sumphi0
        real(kind=kr)                :: siga
        real(kind=kr), allocatable   :: phin(:)
        real(kind=kr), allocatable   :: phino(:)
        real(kind=kr), allocatable   :: s(:)
        real(kind=kr), allocatable   :: a(:)
        real(kind=kr), allocatable   :: b(:)
        real(kind=kr), allocatable   :: c(:)

        if (mod(jmax,p) /= 0) stop ' Fine mesh does not align with CMDSA.'
        if (bc(1) /= 0) stop ' Only vacuum BC supported on left edge.'
        if (bc(2) /= 1) stop ' Only reflective BC supported on right edge.'

        nmax=jmax/p

        allocate(phin(nmax))
        allocate(phino(nmax))
        allocate(s(nmax))
        allocate(a(nmax))
        allocate(b(nmax))
        allocate(c(nmax))
        phin =0.0_kr
        phino=0.0_kr
        s    =0.0_kr
        a    =0.0_kr
        b    =0.0_kr
        c    =0.0_kr

        j=1
        do n=1,nmax
          sumphi =0.0_kr
          sumphi0=0.0_kr
          do jj=1,p
            sumphi =sumphi +phi (j)
            sumphi0=sumphi0+phi0(j)
            j=j+1
          enddo
          phin (n)=sumphi /p
          phino(n)=sumphi0/p
          s(n)    =sigs*p*h*(phin(n)-phino(n))
        enddo

        dc  =1.0_kr/(3.0_kr*p*h*sigt)
        siga=sigt-sigs
        b(1)=0.5_kr+dc+p*h*siga
        c(1)=-dc
        do n=2,nmax-1
          a(n)=-dc
          b(n)=2.0_kr*dc+p*h*siga
          c(n)=-dc
        enddo
        a(nmax)=-dc
        b(nmax)=dc+p*h*siga

        call tdma(nmax,a,b,c,s)

        j=1
        do n=1,nmax
          do jj=1,p
            phi(j) =phi(j) +s(n)
            j=j+1
          enddo
        enddo

        deallocate(phin)
        deallocate(phino)
        deallocate(s)
        deallocate(a)
        deallocate(b)
        deallocate(c)

      end subroutine cmdsa

      subroutine tdma(n,a,b,c,d)

        use global

        implicit none

        integer(4),    intent(in)    :: n
        real(kind=kr), intent(in)    :: a(n)
        real(kind=kr), intent(inout) :: b(n)
        real(kind=kr), intent(in)    :: c(n)
        real(kind=kr), intent(inout) :: d(n)

        integer(4)                   :: k
        integer(4)                   :: km1
        real(kind=kr)                :: xm

        ! simple solution for n = 1

        if (n == 1) then
          d(1) = d(1)/b(1)
          return
        endif

        ! forward substitution

        do k=2,n
          km1 =k-1
          xm  =a(k)/b(km1)
          b(k)=b(k)-xm*c(km1)
          d(k)=d(k)-xm*d(km1)
        enddo

        ! back substitition

        d(n)=d(n)/b(n)
        do k=n-1,1,-1
          d(k)=(d(k)-c(k)*d(k+1))/b(k)
        enddo

      end subroutine tdma

    end program main
