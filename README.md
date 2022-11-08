# space

	DOUBLE PRECISION v_t,x,z,Xtest,ZTEST,ttest,dtrecord
	DOUBLE PRECISION dt,dtmin,dttest
	EXTERNAL derivs,bsstep
	parameter (np=1000,tmax=1.0d4,dtmin=1.d-3,dttest=0.01d0)
	parameter (v_t=0.07d0,nrecord=60)
	DOUBLE PRECISION vp0(np,3),vp(np,3),var(6)
	DOUBLE PRECISION xp0(np),yp0(np),zp0(np),t0(np),xp(np),yp(np),zp(np)
	DOUBLE PRECISION vtmp0(3),vtmp(3),ring_theta,t,gasdev
	DOUBLE PRECISION en0(np),en(np),trecord(nrecord)
	DOUBLE PRECISION en_perp, en_parl,v_parl,e1,e2

	DOUBLE PRECISION vtest(3),bfld(2),efld(2)
	character*60 filename

	real(8) vxx,vyy,vzz,Tx,Ty,Tz,beta,vxxx,vyyy,vzzz,Tpara
        
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+ Parameter setting about the turblence of Alfven waves, where
C+ knum is the total number of waves used, kz is the wave number,
C+ psi is the wave phase, amp is the wave amplitude,theta is the
C+ wave propagation angle. Omega_i is ion normalized frequency. 
C+ Where, the frequency, the wave field strength and the speed are
C+ normalized with the proton gyrofrequency, the background magnetic
C+ field strength and the Alfven speed, respectively.
C+
C+Omega_i_L,  vphL  and kzL is the parameters of the left propagating waves
C+
C+ Proton: Omega_i=1.0d0
C+ He++:   Omega_i=2/4=0.5d0
C+ He+:    Omega_i=1/4=0.25d0
C+ O5+:    Omega_i=5/16=0.3125d0
C+ Mg9+:   Omega_i=9/24.3=0.3704d0  
C+   
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	DOUBLE PRECISION eps,pi,theta,omegamax,omegamin,Omega_i
	DOUBLE PRECISION Bindex,Bnormal,Btest,AL,AR
	parameter (omegamin=0.09d0,omegamax=0.13d0,Bindex=1.6667D0)
	parameter (eps=0.15d0,pi=3.1415925635897932385d0)
	DOUBLE PRECISION amp(200),domega_k,ampL(200)
	DOUBLE PRECISION omega_k(200),kz(200),vph(200),psi(200),randm(200)
	DOUBLE PRECISION omega_kL(200),kzL(200),vphL(200),psiL(200),randmL(200)
      integer knum
	common /wave/psi,amp,kz,omega_k,vph,theta,kzL,omega_kL,vphL,psiL,ampL
      common /knumber/ knum
	common /time/dt,Omega_i

	Omega_i=0.3125D0 ! O5+
      knum=40
	t=0.0d0
	dt=0.01d0
      AR=1.0d0! the amplitude of right propational waves
      AL=1.0d0! the amplitude of left propational waves
	theta=0.d0 !pi/4.0d0
	ring_theta=pi/4.0d0

	do i=1,10
		trecord(i)=20.D0*dble(i)
	enddo
	do i=11,nrecord-1
		trecord(i)=INT((tmax-trecord(10))/dble(nrecord-10))*dble(i-10)
     &               +trecord(10)
	    print*,trecord(i)
	enddo
	trecord(nrecord)=tmax-dt


C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+ Initialize the wave paprameters, where the wave phase is random.
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c	call seed(5000)
c	call random_seed(size=m)
c	print*,m
c	pause	

	call random_number (randm)
      call random_number (randmL)
      

	open(20,file='Bparameters.dat')
	write(20,'(1x,(g15.7),f10.4)') theta,eps
	do i=1, knum
          psi(i)=2.d0*pi*randm(i)
          psiL(i)=2.d0*pi*randmL(i)
      enddo
    


	domega_k=(dlog10(omegamax)-dlog10(omegamin))/dble(knum-1)


	do i=1,knum
          omega_k(i)=10**(dlog10(omegamin)+domega_k*dble(i-1))

          kz(i)=omega_k(i)

	    vph(i)=1.d0
	    write(*,*) i, omega_k(i), kz(i)
	enddo
	do i=1,knum
          omega_kL(i)=omega_k(i)

          kzL(i)=-omega_kL(i)

	    vphL(i)=omega_kL(i)/kzL(i)
	    write(*,*) i, omega_kL(i), kzL(i)
	enddo


	Bnormal=eps*(1.d0-Bindex)
     &		/(omegamax**(1.d0-Bindex)-omegamin**(1.d0-Bindex))
	write(*,*) omegamax,omegamin,Bnormal
	Btest=0.d0
	do i=2,knum-1
		amp(i)=dsqrt(Bnormal*0.5*(omega_k(i+1)-omega_k(i-1))
     &			*(omega_k(i)**(-Bindex)))
		Btest=Btest+amp(i)**2

	enddo
	amp(1)=dsqrt(Bnormal*(omega_k(2)-omega_k(1))
     &		*(omega_k(1)**(-Bindex)))
	amp(knum)=dsqrt(Bnormal*(omega_k(knum)-omega_k(knum-1))
     &				*(omega_k(knum)**(-Bindex)))
	Btest=Btest+amp(1)**2+amp(knum)**2
      Btest2=0.
	do i=1,knum
		amp(i)=dsqrt(eps/Btest)*amp(i)
                    
          ampL(i)=AL*amp(i)
          amp(i)=AR*amp(i)
          
	    Btest2=Btest2+amp(i)**2
		write(20,21) i,omega_k(i),kz(i),amp(i),vph(i),
     &			psi(i),omega_k(i)-kz(i)
	enddo
	write(*,*) 'eps, Btest=',eps,Btest,Btest2
21	format(1x,i3,6(g15.5))
	close(20)

	open(10,file='Bwave.dat')
	do z=-1.d3,1.d3,1d0

			call Bwave(x,z,0.d0,bfld,efld)
			write(10,20) z,bfld,dsqrt(bfld(1)**2+bfld(2)**2),
     &                     efld,dsqrt(efld(1)**2+efld(2)**2)

	enddo
20	format(1x,8(g15.5))
	close(10)
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



	open(2,file='single.dat')
      ifile=eps*100
	jfile=omegamin*100
      write(FILENAME,'(''thermal_eps='',i3.3,''_k='',i3.3,''.dat'')') 
     &ifile,jfile
      open(40,file=filename)







C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+ Giving the initial velocity and position of the test particles.+
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	open(1,file='v.dat')
	en_perp=0.d0
	en_parl=0.d0
	v_parl=0.d0
	iflag=0
	do i=1,np
		vtest(1)=gasdev(v_t)
		vtest(2)=gasdev(v_t)
		vtest(3)=gasdev(v_t)

		
		do j=1,3
			vp0(i,j)=vtest(j)
		enddo
		CALL RANDOM_NUMBER(XTEST)
		xp0(i)=XTEST*3.D3
		CALL RANDOM_NUMBER(ZTEST)
		zp0(i)=ZTEST*3.D3
		CALL RANDOM_NUMBER(tTEST)
		t0(i)=0.d0

		
		e1=vp0(i,1)**2+vp0(i,2)**2
		e2=vp0(i,3)**2
		en0(i)=e1+e2

		v_parl=v_parl+vp0(i,3)

		write(1,*),(vp0(i,j),j=1,3)
	enddo
	v_parl=v_parl/dble(np)
	do i=1,np
		e1=vp0(i,1)**2+vp0(i,2)**2
		e2=(vp0(i,3)-v_parl)**2
		en_perp=en_perp+e1
		en_parl=en_parl+e2
	enddo
	en_perp=en_perp/dble(np)
	en_parl=en_parl/dble(np)
	close(1)


	ifile=0

	write(FILENAME,'(''V_all_n'',I3.3,''.dat'')') ifile
	open(3,file=filename)
	write(3,61),t
	do i=1,np
		write(3,60),xp0(i),zp0(i),(vp0(i,j),j=1,3),
     &			(en(i)-en0(i))/en0(i)
	enddo
	write(*,*),'Writting File, number of particles:',np,NINT(t)
	close(3)

	open(94,file='V_average.dat')
	write(94,63),t,en_perp,en_parl,v_parl+1.d0
	
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+ Calculate the particle velocity and position at time t.      +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C	do isingle=1,10
C	    write(FILENAME,'(''V_single_n'',I3.3,''.dat'')') isingle
C	    open(isingle,file=filename)
C	enddo
C	print*,'dt=',dt
C	print*,'eps=',eps
C	print*,'Total number of test particles: ', np

	ifile=1
	do t=0.d0,tmax,dt
          
          write(2,69),t,vp0(1,1),vp0(1,2),vp0(1,3),xp0(1),yp0(1),zp0(1)
                        

          do i=1,np
              

			if (t0(i).le.t) then
				do j=1,3
					var(j)=vp0(i,j)
				enddo
				var(4)=xp0(i)
				var(5)=yp0(i)
				var(6)=zp0(i)
	            call odeint(var,6,t,t+dt,1.d-12,dttest,dtmin,nok,
     &				        nbad,derivs,bsstep)


				do j=1,3
					vp(i,j)=var(j)
				enddo
				xp(i)=var(4)
                yp(i)=var(5)
				zp(i)=var(6)
				do j=1,3
					vp0(i,j)=vp(i,j)
				enddo
				xp0(i)=xp(i)
                yp0(i)=yp(i)
			    zp0(i)=zp(i)
			else
				do j=1,3
					vp(i,j)=vp0(i,j)
				enddo				
				xp(i)=xp0(i)
                yp(i)=yp0(i)
				zp(i)=zp0(i)
			endif

		enddo






		if (dabs(INT(t/1.d0)*1.d0-t) .lt. dt) then



	    vxx=0
	    vyy=0
	    vzz=0
          vpara=0
	    do i=1,np
		vxx=vxx+vp(i,1)
		vyy=vyy+vp(i,2)
		vzz=vzz+vp(i,3)
          
	    
	    end do

          vxxx=0
          vyyy=0
	    vzzz=0
          vxx=vxx/np
		  vyy=vyy/np
		  vzz=vzz/np
	    do i=1,np
		    vxxx=vxxx+(vp(i,1)-vxx)**2
		    vyyy=vyyy+(vp(i,2)-vyy)**2
		    vzzz=vzzz+(vp(i,3)-vzz)**2

	    end do	  
          Tx=vxxx/np
	    Ty=vyyy/np
	    Tz=vzzz/np
          Tpara=(Tx+Ty)/2
          
          write(40,51),t,Tx,Ty,Tz,Tpara,vzz
          
	    vxx=0
	    vyy=0
	    vzz=0
	    vpara=0
	    vxxx=0
	    vyyy=0
	    vzzz=0
		
			

	    
	    endif
	    
51        format(f10.2,5(1x,f13.4))









		if (t .le. 1.D2) then
			dtrecord=0.5D0
		else
			dtrecord=5.D0
		endif
		if (dabs(INT(t/dtrecord)*dtrecord-t) .lt. dt) then
			en_perp=0.d0
			en_parl=0.d0
			v_parl=0.d0
			do i=1,np
				v_parl=v_parl+vp(i,3)
			enddo
			v_parl=v_parl/dble(np)
			do i=1,np
				en_perp=en_perp+vp(i,1)**2+vp(i,2)**2
				en_parl=en_parl+(vp(i,3)-v_parl)**2
			enddo
			en_perp=en_perp/dble(np)
			en_parl=en_parl/dble(np)
			write(94,63),t,en_perp,en_parl,v_parl+1.d0
		endif

		if (dabs(t-trecord(ifile)).lt. dt) then
			write(FILENAME,'(''V_all_n'',I3.3,''.dat'')') ifile
			open(93,file=filename)
			write(93,61),t
			mm=0
			do i=1,np
				if (t0(i).le.t) then
					en(i)=vp(i,1)**2+vp(i,2)**2+vp(i,3)**2
					write(93,60),xp(i),zp(i),(vp(i,j),j=1,3),
     &					(en(i)-en0(i))/en0(i)
					mm=mm+1
				endif		
			enddo

			write(*,*),'Writting File, Particle Number:',mm,NINT(t)
			close(93)
			ifile=ifile+1
		endif
	enddo
	do isingle=1,10
	    close(isingle)
	enddo
	close(94)
      close(2)
60	format(1x,6(g13.4))
61	format(1x,g13.4)
62	format(1x,I7)
63    format(1x,4(g13.4))
69	format(7(1x,g18.8))
	END



C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+  Subroutine for calculating the wave magnetic field and electric field      +
C+  in plasma frame, where Bw(1), Ew(1) is the x component and Bw(2), Ew(2)    +
C+  is the y component. x and z is the position.                               +
C+  Waves are the Alfven-cyclotron waves with left-hand circular polarization  +
C+  propagating parallel to the ambient magnetic field line.                   +
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 	subroutine Bwave(x,z,t,Bw,Ew)
!	parameter (knum=41)
	DOUBLE PRECISION Bw(2),Ew(2),x,z,t,psi(200),amp(200),ampL(200)
	DOUBLE PRECISION omega_k(200),kz(200),vph(200),theta
	DOUBLE PRECISION omega_kL(200),kzL(200),vphL(200),psiL(200)
      integer knum
      common /wave/psi,amp,kz,omega_k,vph,theta,kzL,omega_kL,vphL,psiL
     &             ,ampL
      common /knumber/ knum
      
	Bw(1)=0.d0
	Bw(2)=0.d0
	Ew(1)=0.d0
	Ew(2)=0.d0


	do i=1,knum

		Bw(2)=Bw(2)+amp(i)*dsin(omega_k(i)*t-kz(i)*z-psi(i))
		Ew(1)=Ew(1)+vph(i)*amp(i)
     &			*dsin(omega_k(i)*t-kz(i)*z-psi(i))

	enddo

	do i=1,knum

		Bw(2)=Bw(2)+ampL(i)*dsin(omega_kL(i)*t-kzL(i)*z-psiL(i))
		Ew(1)=Ew(1)+vphL(i)*ampL(i)
     &			*dsin(omega_kL(i)*t-kzL(i)*z-psiL(i))

	enddo

	return
	end

		

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+  Subroutine for calculating Vx, Vy and Vz at time (n+1) from         +
C+  their values at time (n) and the corresponding magnetic field.      +
C+  An coordinate transformation scheme is used in this subroutine.         +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


	subroutine velocity(V2,V1,Bw)
	DOUBLE PRECISION V2(3),V1(3),Bw(2),dt,B,VA1,dot
	DOUBLE PRECISION Vxp,Vyp,phi,Omega_i
	DOUBLE PRECISION V1_parr(3),V1_perp(3),Vparr1,Vperp1,V2_perp(3)
	DOUBLE PRECISION UnitB(3),Unitxp(3),Unityp(3),Unitzp(3)
	DOUBLE PRECISION Unitmp1(3),Unitperp(3)
	DOUBLE PRECISION pi
	parameter (pi=3.1415925635897932385d0)
	common /time/ dt,Omega_i

C 计算磁场B方向的单位矢量。

	B=dsqrt(Bw(1)**2+Bw(2)**2+1.d0)
	UnitB(1)=Bw(1)/B
	UnitB(2)=Bw(2)/B
	UnitB(3)=1.d0/B

C 计算速度V1的大小及其所在方向的单位矢量。
c	VA1=dsqrt(V1(1)**2+V1(2)**2+V1(3)**2)
c	do i=1,3
c		UnitV1(i)=V1(i)/VA1
c	enddo

C 计算x'方向的单位矢量,也即在速度V1垂直于磁场分量方向上的单位矢量.

	call cross(UnitB,V1,Unitmp1)
	call cross(Unitmp1,UnitB,Unitperp)
	Vperp1=dsqrt(dot(Unitperp,Unitperp))
	do i=1,3
		Unitxp(i)=Unitperp(i)/Vperp1
	enddo	
        
c	print*,'Unitperp=',Unitperp
c	print*,'Unitxp=',Unitxp,dot(Unitxp,Unitxp)
c	pause

C 计算z'方向的单位向量,也即平行于磁场方向的单位向量.

	do i=1,3
		Unitzp(i)=UnitB(i)
	enddo

C 计算y'方向的单位向量

	call cross(Unitzp,Unitxp,Unityp)

C 计算V1平行于磁场方向的速度分量及其大小.

	Vparr1=dot(V1,UnitB)
	do i=1,3
		V1_parr(i)=Vparr1*UnitB(i)
	enddo

C 计算V1垂直于磁场方向的速度大小

	Vperp1=dsqrt(dot(V1,V1)-Vparr1**2)
        
C 计算垂直于磁场方向的速度在dt时间旋转的角度,其中负号对应正的离子.

	phi=-dt*B*0.25d0 
        
C 机算在t+dt时刻,垂直于磁场方向的速度在x'y'z'坐标系中的各分量.

	Vxp=Vperp1*dcos(phi)
	Vyp=Vperp1*dsin(phi)

C 将垂直于磁场方向的速度丛x'y'z'坐标变换到xyz坐标系.

	do i=1,3
		V2_perp(i)=Unitxp(i)*Vxp+Unityp(i)*Vyp
		V2(i)=V2_perp(i)+V1_parr(i)
	enddo
	
	return
	end
        
C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+  Subroutine for calculating Vx, Vy and Vz at time (n+1) from         +
C+  their values at time (n) and the corresponding wave field.          +
C+  The Boris algorithm is used in this subroutine.                     +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	subroutine velocity_Boris(V2,V1,Bw,Ew)
	DOUBLE PRECISION V2(3),V1(3),Bw(2),Ew(2),dt,s,C,Omega_i
	DOUBLE PRECISION Vtmp(3),Vtmp1(3),Vtmp2(3),Vtmp3(3)
	DOUBLE PRECISION UnitA(3),UnitC(3),Uphi(3)
	DOUBLE PRECISION pi
	parameter (pi=3.1415925635897932385d0)
	common /time/ dt,Omega_i

C 计算电场A矢量。

	UnitA(1)=Omega_i*Ew(1)
	UnitA(2)=Omega_i*Ew(2)
	UnitA(3)=0.d0

C 计算磁场C矢量。

	UnitC(1)=Omega_i*Bw(1)
	UnitC(2)=Omega_i*Bw(2)
	UnitC(3)=Omega_i*1.d0
	C=dsqrt(UnitC(1)**2+UnitC(2)**2+UnitC(3)**2)

C 计算电场的前半步贡献

	do i=1,3
		Vtmp1(i)=V1(i)+0.5d0*dt*UnitA(i)
	enddo
        
C 计算磁场的第一次旋转

	call cross(Vtmp1,UnitC,Vtmp)
	do i=1,3
		Vtmp2(i)=Vtmp1(i)+0.5d0*dt*Vtmp(i)
	enddo
        
C 第二次旋转

	s=dt/(1.d0+0.25d0*C*C*dt*dt)
	do i=1,3
		Uphi(i)=s*UnitC(i)
	enddo
	call cross(Vtmp2,Uphi,Vtmp)
	do i=1,3
		Vtmp3(i)=Vtmp1(i)+Vtmp(i)
	enddo
        
C 计算电场的后半步贡献

      do i=1,3
	    V2(i)=Vtmp3(i)+0.5d0*dt*UnitA(i)
      enddo 
	return
	end
        
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



C+ Returns a normally distributed deviate with zero mean and variance. +
C+ V_t/sqrt(2), where V_t is the thermal speed of a gas.               +
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	function gasdev(v_t)
	DOUBLE PRECISION gasdev
	DOUBLE PRECISION v_t
	integer iset
	DOUBLE PRECISION fac,gset,rsq,v1,v2
	save iset,gset
	data iset/0/
	if (iset.eq.0) then

1		call random_number(v1)
		call random_number(v2)
		v1=2.d0*v1-1.d0
		v2=2.d0*v2-1.d0
		rsq=v1**2+v2**2
		if(rsq.ge.1.d0.or.rsq.eq.0.d0) goto 1
		fac=dsqrt(-2.d0*dlog(rsq)/rsq)
		gset=v1*fac*v_t/dsqrt(2.d0)
		gasdev=v2*fac*v_t/dsqrt(2.d0)
		iset=1
	else
		gasdev=gset
		iset=0
	endif
	return
	end
	

C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+  Calculate the cross of two vector.                                   +
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	subroutine cross(vect1,vect2,vect3)
	DOUBLE PRECISION vect1(3),vect2(3),vect3(3)
	vect3(1)=vect1(2)*vect2(3)-vect1(3)*vect2(2)
	vect3(2)=vect2(1)*vect1(3)-vect1(1)*vect2(3)
	vect3(3)=vect1(1)*vect2(2)-vect1(2)*vect2(1)
	return
	end


C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+  Calculate the dot of two vector.                                    +
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	function dot(vect1,vect2)
	DOUBLE PRECISION dot,vect1(3),vect2(3)
	dot=vect1(1)*vect2(1)+vect1(2)*vect2(2)+vect1(3)*vect2(3)
	return
	end


C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+   Computes the derivatives.                                         +
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	SUBROUTINE derivs(t,v,dvdt)
	DOUBLE PRECISION t,v(6),dvdt(6)
	DOUBLE PRECISION Bfld(2),Efld(2),x,z
	DOUBLE PRECISION dt,Omega_i
      integer knum
	EXTERNAL Bwave
	common /time/ dt,Omega_i
      common /knumber/ knum
	x=v(4)
	z=v(6)
	call Bwave(x,z,t,Bfld,Efld)
	if (t .le. 2.d2) then
		do j=1,2
			Bfld(j)=Bfld(j)*dexp(-(t-2.d2)**2/50.d0/50.d0)
			Efld(j)=Efld(j)*dexp(-(t-2.d2)**2/50.d0/50.d0)
		enddo
      endif
      
	dvdt(1)=Omega_i*(Efld(1)+v(2)-v(3)*Bfld(2))
	dvdt(2)=Omega_i*(Efld(2)+v(3)*Bfld(1)-v(1))
	dvdt(3)=Omega_i*(v(1)*Bfld(2)-v(2)*Bfld(1))
	dvdt(4)=v(1)
	dvdt(5)=v(2)
	dvdt(6)=v(3)
	return
	END



      SUBROUTINE odeint(ystart,nvar,x1,x2,eps,h1,hmin,nok,nbad,derivs,
     *bsstep)
      INTEGER nbad,nok,nvar,KMAXX,MAXSTP,NMAX
      DOUBLE PRECISION eps,h1,hmin,x1,x2,ystart(nvar),TINY
      EXTERNAL derivs,bsstep
      PARAMETER (MAXSTP=10000,NMAX=50,KMAXX=200,TINY=1.d-30)
      INTEGER i,kmax,kount,nstp
      DOUBLE PRECISION dxsav,h,hdid,hnext,x,xsav,dydx(NMAX),xp(KMAXX),
     *y(NMAX),yp(NMAX,KMAXX),yscal(NMAX)
      COMMON /path/ kmax,kount,dxsav,xp,yp
      x=x1
      h=dsign(h1,x2-x1)
      nok=0
      nbad=0
      kount=0
      do 11 i=1,nvar
        y(i)=ystart(i)
11    continue
      if (kmax.gt.0) xsav=x-2.d0*dxsav
      do 16 nstp=1,MAXSTP
        call derivs(x,y,dydx)
        do 12 i=1,nvar
          yscal(i)=dabs(y(i))+dabs(h*dydx(i))+TINY
12      continue
        if(kmax.gt.0)then
          if(dabs(x-xsav).gt.dabs(dxsav)) then
            if(kount.lt.kmax-1)then
              kount=kount+1
              xp(kount)=x
              do 13 i=1,nvar
                yp(i,kount)=y(i)
13            continue
              xsav=x
            endif
          endif
        endif
        if((x+h-x2)*(x+h-x1).gt.0.) h=x2-x
        call bsstep(y,dydx,nvar,x,h,eps,yscal,hdid,hnext,derivs)
        if(hdid.eq.h)then
          nok=nok+1
        else
          nbad=nbad+1
        endif
        if((x-x2)*(x2-x1).ge.0.)then
          do 14 i=1,nvar
            ystart(i)=y(i)
14        continue
          if(kmax.ne.0)then
            kount=kount+1
            xp(kount)=x
            do 15 i=1,nvar
              yp(i,kount)=y(i)
15          continue
          endif
          return
        endif
c        if(dabs(hnext).lt.hmin) pause
c     *'stepsize smaller than minimum in odeint'
        if(dabs(hnext).lt.hmin) write(*,*)
     *'stepsize smaller than minimum in odeint'
        h=hnext
16    continue
      pause 'too many steps in odeint'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.



      SUBROUTINE bsstep(y,dydx,nv,x,htry,eps,yscal,hdid,hnext,derivs)
      INTEGER nv,NMAX,KMAXX,IMAX
      DOUBLE PRECISION eps,hdid,hnext,htry,x,dydx(nv),y(nv),yscal(nv),
     *SAFE1,SAFE2,REDMAX,REDMIN,TINY,SCALMX
      PARAMETER (NMAX=50,KMAXX=8,IMAX=KMAXX+1,SAFE1=.25d0,SAFE2=.7d0,
     *REDMAX=1.d-5,REDMIN=.7d0,TINY=1.d-30,SCALMX=.1d0)
CU    USES derivs,mmid,pzextr
      INTEGER i,iq,k,kk,km,kmax,kopt,nseq(IMAX)
      DOUBLE PRECISION eps1,epsold,errmax,fact,h,red,scale,work,wrkmin,
     *xest,xnew,a(IMAX),alf(KMAXX,KMAXX),err(KMAXX),yerr(NMAX),
     *ysav(NMAX),yseq(NMAX)
      LOGICAL first,reduct
      SAVE a,alf,epsold,first,kmax,kopt,nseq,xnew
      EXTERNAL derivs
      DATA first/.true./,epsold/-1.d0/
      DATA nseq /2,4,6,8,10,12,14,16,18/
      if(eps.ne.epsold)then
        hnext=-1.d29
        xnew=-1.d29
        eps1=SAFE1*eps
        a(1)=nseq(1)+1
        do 11 k=1,KMAXX
          a(k+1)=a(k)+nseq(k+1)
11      continue
        do 13 iq=2,KMAXX
          do 12 k=1,iq-1
            alf(k,iq)=eps1**((a(k+1)-a(iq+1))/((a(iq+1)-a(1)+1.d0)*(2*k+
     *1)))
12        continue
13      continue
        epsold=eps
        do 14 kopt=2,KMAXX-1
          if(a(kopt+1).gt.a(kopt)*alf(kopt-1,kopt))goto 1
14      continue
1       kmax=kopt
      endif
      h=htry
      do 15 i=1,nv
        ysav(i)=y(i)
15    continue
      if(h.ne.hnext.or.x.ne.xnew)then
        first=.true.
        kopt=kmax
      endif
      reduct=.false.
2     do 17 k=1,kmax
        xnew=x+h
        if(xnew.eq.x)pause 'step size underflow in bsstep'
        call mmid(ysav,dydx,nv,x,h,nseq(k),yseq,derivs)
        xest=(h/nseq(k))**2
        call pzextr(k,xest,yseq,y,yerr,nv)
        if(k.ne.1)then
          errmax=TINY
          do 16 i=1,nv
            errmax=max(errmax,dabs(yerr(i)/yscal(i)))
16        continue
          errmax=errmax/eps
          km=k-1
          err(km)=(errmax/SAFE1)**(1./(2*km+1))
        endif
        if(k.ne.1.and.(k.ge.kopt-1.or.first))then
          if(errmax.lt.1.)goto 4
          if(k.eq.kmax.or.k.eq.kopt+1)then
            red=SAFE2/err(km)
            goto 3
          else if(k.eq.kopt)then
            if(alf(kopt-1,kopt).lt.err(km))then
              red=1.d0/err(km)
              goto 3
            endif
          else if(kopt.eq.kmax)then
            if(alf(km,kmax-1).lt.err(km))then
              red=alf(km,kmax-1)*SAFE2/err(km)
              goto 3
            endif
          else if(alf(km,kopt).lt.err(km))then
            red=alf(km,kopt-1)/err(km)
            goto 3
          endif
        endif
17    continue
3     red=min(red,REDMIN)
      red=max(red,REDMAX)
      h=h*red
      reduct=.true.
      goto 2
4     x=xnew
      hdid=h
      first=.false.
      wrkmin=1.d35
      do 18 kk=1,km
        fact=max(err(kk),SCALMX)
        work=fact*a(kk+1)
        if(work.lt.wrkmin)then
          scale=fact
          wrkmin=work
          kopt=kk+1
        endif
18    continue
      hnext=h/scale
      if(kopt.ge.k.and.kopt.ne.kmax.and..not.reduct)then
        fact=max(scale/alf(kopt-1,kopt),SCALMX)
        if(a(kopt+1)*fact.le.wrkmin)then
          hnext=h/fact
          kopt=kopt+1
        endif
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.



      SUBROUTINE mmid(y,dydx,nvar,xs,htot,nstep,yout,derivs)
      INTEGER nstep,nvar,NMAX
      DOUBLE PRECISION htot,xs,dydx(nvar),y(nvar),yout(nvar)
      EXTERNAL derivs
      PARAMETER (NMAX=50)
      INTEGER i,n
      DOUBLE PRECISION h,h2,swap,x,ym(NMAX),yn(NMAX)
      h=htot/nstep
      do 11 i=1,nvar
        ym(i)=y(i)
        yn(i)=y(i)+h*dydx(i)
11    continue
      x=xs+h
      call derivs(x,yn,yout)
      h2=2.d0*h
      do 13 n=2,nstep
        do 12 i=1,nvar
          swap=ym(i)+h2*yout(i)
          ym(i)=yn(i)
          yn(i)=swap
12      continue
        x=x+h
        call derivs(x,yn,yout)
13    continue
      do 14 i=1,nvar
        yout(i)=0.5d0*(ym(i)+yn(i)+h*yout(i))
14    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.


      SUBROUTINE pzextr(iest,xest,yest,yz,dy,nv)
      INTEGER iest,nv,IMAX,NMAX
      DOUBLE PRECISION xest,dy(nv),yest(nv),yz(nv)
      PARAMETER (IMAX=13,NMAX=50)
      INTEGER j,k1
      DOUBLE PRECISION delta,f1,f2,q,d(NMAX),qcol(NMAX,IMAX),x(IMAX)
      SAVE qcol,x
      x(iest)=xest
      do 11 j=1,nv
        dy(j)=yest(j)
        yz(j)=yest(j)
11    continue
      if(iest.eq.1) then
        do 12 j=1,nv
          qcol(j,1)=yest(j)
12      continue
      else
        do 13 j=1,nv
          d(j)=yest(j)
13      continue
        do 15 k1=1,iest-1
          delta=1.d0/(x(iest-k1)-xest)
          f1=xest*delta
          f2=x(iest-k1)*delta
          do 14 j=1,nv
            q=qcol(j,k1)
            qcol(j,k1)=dy(j)
            delta=d(j)-q
            dy(j)=f1*delta
            d(j)=f2*delta
            yz(j)=yz(j)+dy(j)
14        continue
15      continue
        do 16 j=1,nv
          qcol(j,iest)=dy(j)
16      continue
      endif
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software ,4-#.


