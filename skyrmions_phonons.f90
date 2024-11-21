  !----------------------!
  MODULE systemvariables
  !----------------------!
  
  INTEGER :: N,Ns,Nc,Neql,Nmes,Nt,Ntemp,Nhz,Ncopy,NcopySz
  PARAMETER  (N=30,Nc=N*N,Ns=N*N,Neql=13*Ns,Nmes=4*Neql,Nt=5,Ntemp=60,Nhz=1,Ncopy=1,NcopySz=5,Nlines=100)
  INTEGER, allocatable :: NN1(:,:)
  ! S(j,M): global matrix.  j-th component of spin at the site M model 1
  REAL*8, allocatable :: S(:,:)   
  REAL*8 ::  PI= 4.0d0*DATAN(1.0d0) 
  REAL*8 :: sqrt3=DSQRT(3.0d0)
  REAL*8 :: sqrt3on2=DSQRT(3.0d0)*0.5d0
  COMPLEX*16 ::im=(0.0d0,1.0d0)
  
  END MODULE systemvariables
  !--------------------------!



PROGRAM Metrolopis
    USE systemvariables
    IMPLICIT NONE
    INTEGER :: I,J,K,KK,jh,I2,I3,jt,jn
    INTEGER :: Idum,Iter,It,Ih,Kx,Ky,L,M,M1,qx,qy,Nover, Iloop
    REAL*8 :: T,hz,aux,auxx,q,r
    REAL*8 :: J1,E0,DM1,D,A2,alpha,beta,kappa
    REAL*8 :: deltaS,wB,dE,z,p,H2,Spr,Haux(3),H(3),Sp(3)

    REAL*8 :: Etot,Eav,E2av,C,Stot,S2av,Mag,Mav,M2av,ChiM,Chi,ChiAv,Chirtot,Chir1,Chir2

    REAL*8 :: Idumv(NLines)
    REAL(8) :: Hu,Hup                           ! Hamiltonian terms
    REAL(8), DIMENSION(Ns, 2) :: u  ! u_i for each site, with 2 components (x, y)
    ! Declare vectors for u(i) and u(j) (2D)
    REAL, DIMENSION(2) :: u_i, u_j, u_diff_x, u_diff_y


    CHARACTER(len=300):: NameFile1,NameFile2,NameFile3,NameFile4
    INTEGER :: position
    COMMON Idum
    !----------------
    ALLOCATE(NN1(4,Ns)) ! first neighbors
    ALLOCATE(S(3,Ns))
    !-----------------------------------------
    !   +/- Integer <2147483648=2^31 to ativate generator
    !   loop in Idum to get random seed 
    !-----------------------------------------
    !  Open idums file to load random seed
    !-----------------------------------------
    !NameFile4='idums.txt'
    OPEN(44,FILE='idums', STATUS='OLD', ACTION='READ')

	DO i=1,Nlines
        READ(44,*) Idumv(i)
        WRITE(*,*) Idumv(i)
	ENDDO

    CLOSE(44)

    !-----------------------------------------
    DO Iloop=1, 1
        WRITE(*,*) Iloop,Idum
        Idum=Idumv(Iloop)
        Pi=acos(-1.0d0)
        !-----------------------------------------
        ! subroutine to create the lattice
        !-----------------------------------------
        CALL geometryS()
        !-----------------------------------------
        J1=-1.0d0
        DM1=1.0d0
        A2=0.1d0
        alpha=0.1
        beta=0.1
        kappa=0.5
        !hz=0.5d0
        !-----------------------------------------
        ! loop in magnetic field
    
        DO aux= 1,1
    
           !hz=aux*0.1d0
            hz=VVVd0    
        !-----------------------------------------    
              
        !-----------------------------------------
        ! The beginning
        !-----------------------------------------
        DO Iter=1,Ncopy
                !-----------------------------------------
                ! subroutine to create the file of data
                WRITE(*,*) Iter, hz, J1, DM1, A2, Iloop, NameFile1, NameFile2, NameFile3
                !-----------------------------------------
                CALL CreateFile(Iloop,hz,J1,DM1,alpha,beta,NameFile1,NameFile2,NameFile3)
                !-----------------------------------------
                CALL InitSpin()  ! call initial spin configuration
                !-----------------------------------------
                ! Bucle in temperature
                !-----------------------------------------
                DO jt=1,Ntemp
                    !---------------------------------------
                    T=2.5d0*(0.91d0)**jt
                    !========
                    Nover=5 !
                    IF(T.le.1.5d0) Nover=13
                    deltaS=T
                    IF(T.ge.2.0d0)deltaS=2.0d0
                    !---------------------------------------
                    ! Thermal equilibration ----------------
                    !---------------------------------------
                    ! First write the displacements vectors u_i and it's neighbors
                    !----------------------------------------
                    ! Loop to calculate the u components for each site
                    !----------------------------------------
                    DO i = 1, Ns
                        ! X-component of u_i
                        u(i, 1) = (-J * alpha / kappa) * (S(1,i) * S(1,NN1(1,i)) + S(2,i) * S(2,NN1(1,i)) + S(3,i) * S(3,NN1(1,i)))  + &  
                            (J * alpha / kappa) * (S(1,i) * S(1,NN1(3,i)) + S(2,i) * S(2,NN1(3,i)) + S(3,i) * S(3,NN1(3,i))) -  &  
                            (D * beta / kappa) * (S(2,i) * S(3,NN1(1,i)) - S(2,NN1(1,i)) * S(3,i) + S(1,NN1(1,i)) * S(3,i) - S(1,i) * S(3,NN1(1,i))) + &
                            (D * beta / kappa) * (S(2,i) * S(3,NN1(3,i)) - S(2,NN1(3,i)) * S(3,i) + S(1,NN1(3,i)) * S(3,i) - S(1,i) * S(3,NN1(3,i)))

                        ! Y-component of u_i 
                        u(i, 2) = (-1)*(J * alpha / kappa) * (S(1,i) * S(1,NN1(2,i)) + S(2,i) * S(2,NN1(2,i)) + S(3,i) * S(3,NN1(2,i))) +   &  
                            (J * alpha / kappa) * (S(1,i) * S(1,NN1(4,i)) + S(2,i) * S(2,NN1(4,i)) + S(3,i) * S(3,NN1(4,i))) -   &  
                            (D * beta / kappa) * (S(2,i) * S(3,NN1(2,i)) - S(2,NN1(2,i)) * S(3,i) + S(1,NN1(2,i)) * S(3,i) - S(1,i) * S(3,NN1(2,i))) +  &  
                            (D * beta / kappa) * (S(2,i) * S(3,NN1(4,i)) - S(2,NN1(4,i)) * S(3,i) + S(1,NN1(4,i)) * S(3,i) - S(1,i) * S(3,NN1(4,i)))  
                            
                    ENDDO

                    WRITE(*,*) 'u vectors are OK, ux=', u(1,1)

                   ! Initialize the Hamiltonian Hu, contribution of terms depending on vectors to 0
            Hu = 0.0d0

            ! Loop over all sites
			DO i = 1, Ns

                ! First neighbor (NN1)
                jn = NN1(1, i)
                Hu = Hu + (J * alpha) * u(i, 1) * (S(1, i) * S(1, jn) + S(2, i) * S(2, jn) + S(3, i) * S(3, jn)) + &
                             (D * beta)   * u(i, 1) * ( (S(2, i) * S(3, jn) - S(3, i) * S(2, jn)) + &
                                                     (S(3, i) * S(1, jn) - S(1, i) * S(3, jn)) )

                ! Second neighbor (NN2)
                jn = NN1(2, i)
                Hu = Hu + (J * alpha) * u(i, 2) * (S(1, i) * S(1, jn) + S(2, i) * S(2, jn) + S(3, i) * S(3, jn)) + &
                              (D * beta)  * u(i, 2) * ( (S(2, i) * S(3, jn) - S(3, i) * S(2, jn)) + &
                                                     (S(3, i) * S(1, jn) - S(1, i) * S(3, jn)) )

                ! Third neighbor (NN3)
                jn = NN1(3, i)
                Hu = Hu - (J * alpha) * u(i, 1) * (S(1, i) * S(1, jn) + S(2, i) * S(2, jn) + S(3, i) * S(3, jn)) - &
                             (D * beta) * u(i, 1) * ( (S(2, i) * S(3, jn) - S(3, i) * S(2, jn)) + &
                                                     (S(3, i) * S(1, jn) - S(1, i) * S(3, jn)) )

                ! Fourth neighbor (NN4)
                jn = NN1(4, i)
                Hu = Hu - (J * alpha) * u(i, 2) * (S(1, i) * S(1, jn) + S(2, i) * S(2, jn) + S(3, i) * S(3, jn)) - &
                              (D * beta) * u(i, 2) * ( (S(2, i) * S(3, jn) - S(3, i) * S(2, jn)) + &
                                                     (S(3, i) * S(1, jn) - S(1, i) * S(3, jn)) )

                ! Elastic energy term
                Hu = Hu + 0.5d0 * kappa * (u(i, 1)**2 + u(i, 2)**2)


			END DO

            WRITE(*,*) 'Hamiltonian dep on u Hu=', Hu

                    DO K=1,Neql
                        DO L=1,Ns
                            !Potential energy term
                            !H=(-1/2 * kappa) * (J**2 * alpha**2) * (S(1, M) * S(1, k) + S(2, M) * S(2, k) + S(3, M) * S(3, k))**2 + &
                            !(-J * alpha * beta / kappa) * (S(1, M) * S(1, k) + S(2, M) * S(2, k) + S(3, M) * S(3, k)) * &
                            !    (DMI(1) * (S(2, M) * S(3, l) - S(2, l) * S(3, M)) + &
                            !     DMI(2) * (S(1, l) * S(3, M) - S(1, M) * S(3, l)) + &
                            !     DMI(3) * (S(1, M) * S(2, l) - S(1, l) * S(2, M))) + &
                            !(-1/2 * kappa) * (D**2 * beta**2) * &
                            !    ((S(2, M) * S(3, l) - S(2, l) * S(3, M))**2 + &
                            !     (S(1, l) * S(3, M) - S(1, M) * S(3, l))**2 + &
                            !     (S(1, M) * S(2, l) - S(1, l) * S(2, M))**2)
                            CALL random(p)
                            M=1+INT(ABS(p)*Ns) !choose a random site
                            H(1)=J1*SUM(S(1,NN1(:,M)))+DM1*(-S(3,NN1(2,M))+S(3,NN1(4,M)) )
                            H(2)=J1*SUM(S(2,NN1(:,M)))+DM1*(S(3,NN1(1,M))-S(3,NN1(3,M)) )
                            H(3)=J1*SUM(S(3,NN1(:,M)))+DM1*(S(1,NN1(2,M))-S(1,NN1(4,M))-S(2,NN1(1,M))+S(2,NN1(3,M)) )-hz
                            !New spin ====
                            !===============
                            CALL NewSpin(M,Sp,deltaS)
                            ! Loop over all sites i to calculate Hu', (Hu evaluated in Sp)
            			    DO i = 1, Ns
                                ! Loop over neighbors of site i and calculate vector u(Sp')
                                ! X-component of u_i
                                u(i, 1) = (-J * alpha / kappa) * (Sp(1) * S(1,NN1(1,i)) + Sp(2) * S(2,NN1(1,i)) + Sp(3) * S(3,NN1(1,i))) + &  
                                          (J * alpha / kappa) * (Sp(1) * S(1,NN1(3,i)) + Sp(2) * S(2,NN1(3,i)) + Sp(3) * S(3,NN1(3,i))) + &  
                                          (-D * beta / kappa) * (Sp(2) * S(3,NN1(1,i)) - S(2,NN1(1,i)) * Sp(3) + S(1,NN1(1,i)) * Sp(3) - Sp(1) * S(3,NN1(1,i))) + &  
                                          (D * beta / kappa) * (Sp(2) * S(3,NN1(3,i)) - S(2,NN1(3,i)) * Sp(3) + S(1,NN1(3,i)) * Sp(3) - Sp(1) * S(3,NN1(3,i)))
                                ! Y-component of u_i
                                u(i, 2) = (-1) * (J * alpha / kappa) * (Sp(1) * S(1,NN1(2,i)) + Sp(2) * S(2,NN1(2,i)) + Sp(3) * S(3,NN1(2,i))) &  
                                          + (J * alpha / kappa) * (Sp(1) * S(1,NN1(4,i)) + Sp(2) * S(2,NN1(4,i)) + Sp(3) * S(3,NN1(4,i))) &  
                                          - (D * beta / kappa) * (Sp(2) * S(3,NN1(2,i)) - S(2,NN1(2,i)) * Sp(3) + S(1,NN1(2,i)) * Sp(3) - Sp(1) * S(3,NN1(2,i))) &  
                                          + (D * beta / kappa) * (Sp(2) * S(3,NN1(4,i)) - S(2,NN1(4,i)) * Sp(3) + S(1,NN1(4,i)) * Sp(3) - Sp(1) * S(3,NN1(4,i)))


                                ! Compute u_i - u_j in the desired format
                                ! First compute the differences in x and y components
                                u_diff_x = u(i, 1) - u(j, 1)
                                u_diff_y = u(i, 2) - u(j, 2)

                               ! First neighbor (NN1)
                                jn = NN1(1, i)
                                Hup = Hup + (J * alpha) * u(i, 1) * (Sp(1) * S(1, jn) + Sp(2) * S(2, jn) + Sp(3) * S(3, jn)) + &
                                             (D * beta)   * u(i, 1) * ( (Sp(2) * S(3, jn) - Sp(3) * S(2, jn)) + &
                                                                     (Sp(3) * S(1, jn) - Sp(1) * S(3, jn)) )

                                ! Second neighbor (NN2)
                                jn = NN1(2, i)
                                Hup = Hup + (J * alpha) * u(i, 2) * (Sp(1) * S(1, jn) + Sp(2) * S(2, jn) + Sp(3) * S(3, jn)) + &
                                              (D * beta)  * u(i, 2) * ( (Sp(2) * S(3, jn) - Sp(3) * S(2, jn)) + &
                                                                     (Sp(3) * S(1, jn) - Sp(1) * S(3, jn)) )

                                ! Third neighbor (NN3)
                                jn = NN1(3, i)
                                Hup = Hup - (J * alpha) * u(i, 1) * (Sp(1) * S(1, jn) + Sp(2) * S(2, jn) + Sp(3) * S(3, jn)) - &
                                             (D * beta) * u(i, 1) * ( (Sp(2) * S(3, jn) - Sp(3) * S(2, jn)) + &
                                                                     (Sp(3) * S(1, jn) - Sp(1) * S(3, jn)) )

                                ! Fourth neighbor (NN4)
                                jn = NN1(4, i)
                                Hup = Hup - (J * alpha) * u(i, 2) * (Sp(1) * S(1, jn) + Sp(2) * S(2, jn) + Sp(3) * S(3, jn)) - &
                                              (D * beta) * u(i, 2) * ( (Sp(2) * S(3, jn) - Sp(3) * S(2, jn)) + &
                                                                     (Sp(3) * S(1, jn) - Sp(1) * S(3, jn)) )

                                ! Elastic energy term
                                Hup = Hup + 0.5d0 * kappa * (u(i, 1)**2 + u(i, 2)**2)

            			     END DO

                        WRITE(*,*) 'Hamiltonian dep on u(Sp) Hup=', Hup


                            dE=dot_product(Sp(:)-S(:,M),H(:))    - (Hup-Hu)
    
                            IF(dE.le.0.0)THEN
                                S(:,M)=Sp(:)
                            ELSE
                                wB=EXP(-dE/T)
                                CALL random(p)
                                IF(ABS(p).le.wB) THEN
                                    S(:,M)=Sp(:)
                                ELSE
                                    Haux(1)= H(2)
                                    Haux(2)=-H(1)
                                    Haux(3)=0.0d0
                                    H2=dot_product(Haux,Haux)
                                    Spr=2.d0*dot_product(S(:,M),Haux(:))/H2
                                    S(:,M)=S(:,M)-Spr*Haux(:)
                                ENDIF!IF(ABS(p).le.wB) THEN
                            ENDIF!IF(dE.le.0.0d0)THEN
                        ENDDO !DO L=1,Ns
                    ENDDO 
                    !-EEND->DO K=1,Neql
                    !---------------------------------------
                    ! Normalization of spins (just in case) 
                    !---------------------------------------
                        DO M=1,Ns
                                z=sqrt(dot_product(S(:,M),S(:,M)))
                                S(:,M)=S(:,M)/z
                        ENDDO 
                        !EEND->DO M=1,Ns
                        Eav=0.d0 ! energy
                        Mav=0.d0 ! magnetization along the field	
                        E2av=0.d0 ! energy
                        M2av=0.d0 ! magnetization along the field
                        ChiAv=0.0d0
                  	    Chir1=0
 	    	            Chir2=0
 	  	                Chirtot=0

                        !---------------------------------------------------------------------------------------
                        ! Measurements -------------------------------------------------------------------------
                        !---------------------------------------------------------------------------------------
                        DO K=1,Nmes
                                !---------------------------------------------------------------------------------------
                                ! Overrelaxation -----------------------------------------------------------------------
                                !---------------------------------------------------------------------------------------
                                CALL random(p)
                                M1=1+INT(ABS(p)*Ns)!choose a random site
                                DO Kx=1,Nover
                                    DO L=0,Ns-1
                                            M=MODULO(M1+L,Ns)
                                            IF(M.eq.0) M=Ns
                                            !===============
                                            H(1)=J1*SUM(S(1,NN1(:,M)))+DM1*(-S(3,NN1(2,M))+S(3,NN1(4,M)) )
                                            H(2)=J1*SUM(S(2,NN1(:,M)))+DM1*(S(3,NN1(1,M))-S(3,NN1(3,M)) )
                                            H(3)=J1*SUM(S(3,NN1(:,M)))+DM1*(S(1,NN1(2,M))-S(1,NN1(4,M))-S(2,NN1(1,M))+S(2,NN1(3,M)) )-hz
                                            !===============
                                            !H(1)=Ja(M,1)*S(1,NN1(1,M))+ &
                                            !Ja(M,2)*S(1,NN1(2,M))+ &
                                            !Ja(NN1(3,M),1)*S(1,NN1(3,M))+ &
                                            !Ja(NN1(4,M),2)*S(1,NN1(4,M)) + &
                                            !(DMIa(M,2)*-S(3,NN1(2,M)))+(DMIa(NN1(4,M),2)*S(3,NN1(4,M)))
                                          !
                                            !H(2)=Ja(M,1)*S(2,NN1(1,M))+ &
                                            !Ja(M,2)*S(2,NN1(2,M))+ &
                                            !Ja(NN1(3,M),1)*S(2,NN1(3,M))+ &
                                            !Ja(NN1(4,M),2)*S(2,NN1(4,M)) + &
                                            !(DMIa(M,1)*S(3,NN1(1,M)))-(DMIa(NN1(3,M),1)*S(3,NN1(3,M)))
            !
                                            !H(3)=Ja(M,1)*S(3,NN1(1,M))+ &
                                            !Ja(M,2)*S(3,NN1(2,M))+ &
                                            !Ja(NN1(3,M),1)*S(3,NN1(3,M))+ &
                                            !Ja(NN1(4,M),2)*S(3,NN1(4,M)) + &
                                            !(DMIa(M,2)*S(1,NN1(2,M)))-(DMIa(NN1(4,M),2)*S(1,NN1(4,M))) + &
                                            !(DMIa(M,1)*-S(2,NN1(1,M)))+(DMIa(NN1(3,M),1)*S(2,NN1(3,M))) - hz
                                            !===============
                                            Haux(1)= H(2)
                                            Haux(2)=-H(1)
                                            Haux(3)=0.0d0
                                            H2=dot_product(Haux,Haux)
                                            Spr=2.d0*dot_product(S(:,M),Haux(:))/H2
                                            S(:,M)=S(:,M)-Spr*Haux(:)
                                    ENDDO !DO L=0,Ns-1
                                ENDDO !DO Kx=1,Nover
                                !-----------------------------------
                                ! Metropolis step ------------------
                                !-----------------------------------
                                DO L=1,Ns*Nt
                                        CALL random(p)
                                        M=1+int(abs(p)*Ns)!choose a random site

                                        H(1)=J1*SUM(S(1,NN1(:,M)))+DM1*(-S(3,NN1(2,M))+S(3,NN1(4,M)) )
                                        H(2)=J1*SUM(S(2,NN1(:,M)))+DM1*(S(3,NN1(1,M))-S(3,NN1(3,M)) )
                                        H(3)=J1*SUM(S(3,NN1(:,M)))+DM1*(S(1,NN1(2,M))-S(1,NN1(4,M))-S(2,NN1(1,M))+S(2,NN1(3,M)) )-hz
                                        !  New spin ====
                                        !===============
                                        CALL NewSpin(M,Sp,deltaS)
                                        dE=dot_product(Sp(:)-S(:,M),H(:))    - (Hup-Hu)
                                        if(dE.le.0.) then
                                                S(:,M)=Sp(:)
                                        else
                                                wB=exp(-dE/T)
                                                CALL random(p)
                                                IF(ABS(p).le.wB) THEN
                                                        S(:,M)=Sp(:)
                                                ELSE
                                                        Haux(1)= H(2)
                                                        Haux(2)=-H(1)
                                                        Haux(3)=0.0d0
                                                        H2=dot_product(Haux,Haux)
                                                        Spr=2.d0*dot_product(S(:,M),Haux(:))/H2
                                                        S(:,M)=S(:,M)-Spr*Haux(:)
                                                ENDIF!IF(ABS(p).le.wB) THEN
                                        ENDIF
                                ENDDO!DO L=0,Ns-1
                                ! ----------------------------------
                                ! Magnetization along the field ----
                                Mag=SUM(S(3,:))
                                Mav=Mav+Mag/float(Nmes)
                                ! ----------------------------------
                                ! Energy ---------------------------
                                ! ----------------------------------
                                Etot=0.0d0
                                DO M=1,Ns

                                    H(1)=J1*SUM(S(1,NN1(:,M)))+DM1*(-S(3,NN1(2,M))+S(3,NN1(4,M)) )
                                    H(2)=J1*SUM(S(2,NN1(:,M)))+DM1*(S(3,NN1(1,M))-S(3,NN1(3,M)) )
                                    H(3)=J1*SUM(S(3,NN1(:,M)))+DM1*(S(1,NN1(2,M))-S(1,NN1(4,M))-S(2,NN1(1,M))+S(2,NN1(3,M)) )-hz




                                    Etot=Etot+dot_product(S(:,M),H(:)) + Hu
                                ENDDO !DO M=1,Ns

                                Etot=0.5d0*Etot
                                Eav=Eav+Etot/float(Nmes)
                                E2av=E2av+Etot**2/float(Nmes)
                                !--------------Files Sz
                                !------------------------------------------------			
	                       !scalar producs and chirality
	   		              DO M=1,Ns
	          	 	       !Quiralidad del skyrmion
          	 	 	       Chir1=Chir1 + (S(1,M)*S(2,NN1(1,M))*S(3,NN1(2,M)) - S(1,M)*S(3,NN1(1,M))*S(2,NN1(2,M))  &
           		      	 	 + S(2,M)*S(3,NN1(1,M))*S(1,NN1(2,M)) - S(2,M)*S(1,NN1(1,M))*S(3,NN1(2,M)) &
           		      		  + S(3,M)*S(1,NN1(1,M))*S(2,NN1(2,M)) - S(3,M)*S(2,NN1(1,M))*S(1,NN1(2,M))      )
           	        
           		      		  Chir2=Chir2 + (S(1,M)*S(2,NN1(3,M))*S(3,NN1(4,M)) - S(1,M)*S(3,NN1(3,M))*S(2,NN1(4,M))  &
           	      		 	 + S(2,M)*S(3,NN1(3,M))*S(1,NN1(4,M)) - S(2,M)*S(1,NN1(3,M))*S(3,NN1(4,M)) &
           	     			   + S(3,M)*S(1,NN1(3,M))*S(2,NN1(4,M)) - S(3,M)*S(2,NN1(3,M))*S(1,NN1(4,M))      )
           	        
                   
            	       		 Chirtot = (Chir1 + Chir2)/(8*PI*Nmes)
          	  	      		 !WRITE(*,*) Chirtot
               
            	        ENDDO!end loop NxN
    
	   
	   	
                        ENDDO!END->DO K=1,Nmes
                        !---------------------------------------------------------------------------------------
                        ! END Measurements -------------------------------------------------------------------------
                        !---------------------------------------------------------------------------------------
!                         ChiM=(M2av-Mav*Mav)/(T*dfloat(Ns))
                        C=(E2av-Eav*Eav)/(T*T*dfloat(Ns))
                        Mav=Mav/float(Ns)
                        Eav=Eav/DFLOAT(Ns)
        			    
    
                      
                        ! --------------------------------------
                        ! End of measurement  ------------------
                        ! --------------------------------------
                        !save data
                        OPEN(UNIT=99,FILE=trim(NameFile1),STATUS="OLD",ACTION="WRITE",POSITION="APPEND")
                        WRITE(99,'(A)') "Magnetic Field (Hz), Temperature (T), Energy (Eav), Magnetization (Mav), Specific Heat (C), Chirality (Chirtot)"
                        WRITE(99,'(30f20.10)')hz,T,Eav,Mav,C,Chirtot
                        WRITE(* ,'(30f15.10)')hz,T,Eav,Mav,C,Chirtot
                        CLOSE(99)
                        ! --------------------------------------
                        OPEN(UNIT=98,FILE=trim(NameFile2),STATUS="OLD",ACTION="WRITE",POSITION="APPEND")
                        WRITE(98,'(f20.7)') T
                        DO I2=1,Ns
                            WRITE(98,'(3f20.7)')S(1,I2),S(2,I2),S(3,I2)
                        ENDDO
                        CLOSE(98)





                        
                        
                        ! --------------------------------------
                        ! End of magnetic temperature Loop  ---
                        ! --------------------------------------
                        
	   
    
                ENDDO 
                !ENDDo DO jt=1,Ntemp

         
            ! --------------------------------------
            ! End of magnetic field Loop  ---
            ! --------------------------------------
        ENDDO

    ENDDO
        ! --------------------------------------
        ! End of Idum loop  ---
        ! --------------------------------------    

    ENDDO 

    !DO Iter=1,Ncopy
  
1999 continue  
  !----------------free memory
  
  deallocate(NN1)
  deallocate(S)
  !----------------
END PROGRAM Metrolopis




! =========+=========+=========+=========+=========+=========+=========+=========+==
! PROGRAM: GEOMETRYS
! TYPE	 : subroutine
! PURPOSE: J1-J2  on a HONEYCOMB lattice
! =========+=========+=========+=========+=========+=========+=========+=========+==
  SUBROUTINE geometryS()
  USE systemvariables
  IMPLICIT NONE
  INTEGER::Kx,Kxp,Kxm,Ky,Kyp,Kym,M,Kxpp,Kxmm,Kypp,Kymm
  INTEGER,EXTERNAL :: position
  
  DO Ky=1,N
    Kyp=MODULO(Ky,N)+1
    Kypp=MODULO(MODULO(Ky,N)+1,N)+1
    Kym=Ky-1
    Kymm=Ky-2
    IF(Ky.eq.1)Kym=N
    IF(Ky.eq.1)Kymm=N-1
    IF(Ky.eq.2)Kymm=N
      
    DO Kx=1,N
        Kxp=MODULO(Kx,N)+1
        Kxpp=MODULO(MODULO(Kx,N)+1,N)+1
        Kxm=Kx-1
        Kxmm=Kx-2
        IF(Kx.eq.1)Kxm=N
        IF(Kx.eq.1)Kxmm=N-1
        IF(Kx.eq.2)Kxmm=N
        !===========================
        M=position(Kx,Ky)
        NN1(1,M)=position(Kxp,Ky)
        NN1(2,M)=position(Kx,Kyp)
        NN1(3,M)=position(Kxm,Ky)
        NN1(4,M)=position(Kx,Kym)
        !===========================
    ENDDO ! Kx
  ENDDO !Ky
  END SUBROUTINE geometryS





! =========+=========+=========+=========+=========+=========+=========+=========+==
! PROGRAM: position
! TYPE	 : function
! PURPOSE: lists the unit cells of the lattice as a linear chain
! =========+=========+=========+=========+=========+=========+=========+=========+==
  INTEGER FUNCTION position(Kx,Ky)
  USE systemvariables
  IMPLICIT NONE
  INTEGER :: Kx,Ky,M
    position=1+(Kx-1)+N*(Ky-1)
  END FUNCTION position


! =========+=========+=========+=========+=========+=========+=========+=========+==
! PROGRAM: New spin
! TYPE	 : subroutine
! PURPOSE: computes components of spin 
! version: 31/05/2011
! =========+=========+=========+=========+=========+=========+=========+=========+==
  SUBROUTINE NewSpin(M,Sp,deltaS)
  USE systemvariables
  IMPLICIT NONE
  INTEGER ::M
  REAL*8 :: Sp(3),deltaS
  REAL*8 :: x,y,z,xl,yl,zl,zt,csn,sns
222  CONTINUE
    CALL random(x)
    CALL random(y)      
    z=x*x+y*y      
    IF(z.gt.1.d0) goto 222      
    zl=1.d0-deltaS*z      
    zt=sqrt(deltaS*(1.d0+zl))      
    xl=x*zt     
    yl=y*zt      
    zt=sqrt(S(1,M)**2+S(2,M)**2)    
    csn=S(1,M)/zt      
    sns=S(2,M)/zt       
    Sp(1)=xl*csn*S(3,M)-yl*sns+zl*S(1,M)     
    Sp(2)=xl*sns*S(3,M)+yl*csn+zl*S(2,M)     
    Sp(3)=-xl*zt+zl*S(3,M)     
  END SUBROUTINE NewSpin


! =========+=========+=========+=========+=========+=========+=========+=========+==
! PROGRAM: RANDOM /FROM MARSAGLIA/
! TYPE	 : subroutine
! PURPOSE: GENERATE RANDOM NUMBERS IN [-1,1]
! COMMENT: initialize 'Idum' with (+/-)integer <2147483648=2^31
! version: 01/04/2002
! =========+=========+=========+=========+=========+=========+=========+=========+==

  SUBROUTINE random(rand1)
  implicit none
  integer*4 :: Irand,Ix,Idum
  real*8 :: delta,rand1
  parameter (delta=4.65661287307739D-10)
  common Idum
        Ix=Idum
        Idum=xor(Idum,lshift(Idum,13))
        Idum=xor(Idum,rshift(Idum,17))
        Idum=xor(Idum,lshift(Idum,5))
        Irand=Ix+Idum
        rand1=float(Irand)*delta
  END SUBROUTINE random

! =========+=========+=========+=========+=========+=========+=========+=========+==
! PROGRAM: spin inicial
! TYPE	 : subroutine
! PURPOSE: computes components of spin inicial
! version: 31/05/2011
! =========+=========+=========+=========+=========+=========+=========+=========+==
  SUBROUTINE InitSpin()
  USE systemvariables
  IMPLICIT NONE
  INTEGER :: M
  REAL*8::x,y,z,zt
  
    DO M=1,Ns
        S(1,M)=0.1d0
        S(2,M)=0.05d0
        S(3,M)=DSQRT(1.0d0-S(1,M)**2-S(2,M)**2)
    ENDDO !M    
    WRITE(*,*)'-----------------'
  END SUBROUTINE InitSpin

! =========+=========+=========+=========+=========+=========+=========+=========+==
! PROGRAM: CreateFile
! TYPE	 : subroutine
! PURPOSE: creates the output file for differents parameters
! =========+=========+=========+=========+=========+=========+=========+=========+==
  SUBROUTINE CreateFile(Iloop,hz,J1,DM1,alpha,beta,NameFile1,NameFile2,NameFile3)
  USE systemvariables
  IMPLICIT NONE
  INTEGER :: Iloop,Nimpurezas
  REAL*8 :: J1,DM1,hz,alpha,beta
  CHARACTER(len=300):: NameFile1,NameFile2, NameFile3
  CHARACTER(len=14):: Aux1,Aux2,Aux3,Aux4,Aux5,Aux6,Aux7,Aux8
  
        !======================
        !parameters to write the name of the files
        WRITE(Aux1,*) N
        WRITE(Aux2,*)Iloop
        WRITE(Aux3,'(F8.4)')J1
        WRITE(Aux4,'(F8.4)')DM1
        WRITE(Aux5,'(F8.4)')hz
        WRITE(Aux6,'(F8.4)')alpha
        WRITE(Aux7,'(F8.4)')beta
        
        
        Aux1=adjustl(Aux1)
        Aux2=adjustl(Aux2)
        Aux3=adjustl(Aux3)
        Aux4=adjustl(Aux4)
        Aux5=adjustl(Aux5)
        Aux6=adjustl(Aux6)
        Aux7=adjustl(Aux7)


        !======================

        !  h_ext, M, |Sz|
        ! NameFile  is the name with all of the information about, J2,Delta (D), It, Hz
        NameFile1='dataM/Phonons_mag.L_'//trim(Aux1)//'_'
        NameFile1=trim(NameFile1)//'Ncopy_'//trim(Aux2)//'_'
        NameFile1=trim(NameFile1)//'J1_'//trim(Aux3)//'_'
        NameFile1=trim(NameFile1)//'DM1_'//trim(Aux4)//'_'
        NameFile1=trim(NameFile1)//'hz_'//trim(Aux5)//'_'
        NameFile1=trim(NameFile1)//'alpha_'//trim(Aux6)//'_'
        NameFile1=trim(NameFile1)//'beta_'//trim(Aux7)//'_.dat'
        WRITE(*,'(A)')trim(NameFile1)
        OPEN(99,FILE=trim(NameFile1),STATUS="REPLACE",ACTION="WRITE")
        CLOSE(99)
        
        NameFile2='dataS/S_Phonons.L_'//trim(Aux1)//'_'
        NameFile2=trim(NameFile2)//'Ncopy_'//trim(Aux2)//'_'
        NameFile2=trim(NameFile2)//'J1_'//trim(Aux3)//'_'
        NameFile2=trim(NameFile2)//'DM1_'//trim(Aux4)//'_'
        NameFile2=trim(NameFile2)//'hz_'//trim(Aux5)//'_'
        NameFile2=trim(NameFile2)//'alpha_'//trim(Aux6)//'_'
        NameFile2=trim(NameFile2)//'beta_'//trim(Aux7)//'_.dat'
        WRITE(*,'(A)')trim(NameFile2)
        OPEN(98,FILE=trim(NameFile2),STATUS="REPLACE",ACTION="WRITE")
        CLOSE(98)


  END SUBROUTINE CreateFile
