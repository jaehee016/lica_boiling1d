!C***********************************************************************
!C          MULTIPLE LEVEL-SET AND VOLUME OF FLUID COUPLING METHOD 
!C                           IN REFINED GRID
!C
!C
!C      13.08.28 LEVEL-SET METHOD IS IMPLEMENTED FOR BUBBLY FLOW .
!C      14.09.25 LES WITH VELOCITY RECREATION PART IS IMPLENETED.
!C      15.03.18 DIRICHELT, NEUMANN, PERIODIC BOUNDARY CONDITION IS 
!C               POSSIBLE IN ALL DIRECTION.
!C      15.06.01 MCLS (MASS CONSERVING LEVEL-SET) AND GLOBALL CORRECTION
!C              (ZHANG ET AL., 2010) IS IMPLEMENTED FOR VOLUME CONSERVATION.
!C      15.06.12 MULTIPLE LEVEL-SET METHOD IS IMPLEMENTED PREVENTING 
!C               NUMERICAL COALESCENCE
!C
!C                                                  CODED BY KIYOUNG KIM


!C***********************************************************************
       subroutine Jacobi_bub(a,x,abserr,n)
!""http://ww2.odu.edu/~agodunov/computing/progra"ms/book2/Ch07/Jacobi.f90
!===========================================================
! Evaluate eigenvalues and eigenvectors
! of a real symmetric matrix a(n,n): a*x = lambda*x 
! method: Jacoby method for symmetric matrices 
! Alex G. (December 2009)
!-----------------------------------------------------------
! input ...
! a(n,n) - array of coefficients for matrix A
! n      - number of equations
! abserr - abs tolerance [sum of (off-diagonal elements)^2]
! output ...
! a(i,i) - eigenvalues
! x(i,j) - eigenvectors
! comments ...
!===========================================================
       implicit none
       integer i, j, k, n, ITER_JACOBI
       REAL a(n,n),x(n,n)
       REAL abserr, b2, bar
       REAL beta, coeff, c, s, cs, sc

! initialize x(i,j)=0, x(i,i)=1
! *** the array operation x=0.0 is specific for Fortran 90/95
       x = 0.0
       do i=1,n
        x(i,i) = 1.0
       end do

! find the sum of all off-diagonal elements (squared)
       b2 = 0.0
       do i=1,n
         do j=1,n
           if (i.ne.j) b2 = b2 + a(i,j)**2
         end do
       end do

       IF (B2 .LE. ABSERR) RETURN

       ! average for off-diagonal elements /2
       bar = 0.5*b2/float(n*n)

       DO ITER_JACOBI=1,500 !ITERATION LIMIT
         do i=1,n-1
           do j=i+1,n
             if (a(j,i)**2 <= bar) cycle  ! do not touch small elements
             b2 = b2 - 2.0*a(j,i)**2
             bar = 0.5*b2/float(n*n)
       ! calculate coefficient c and s for Givens matrix
             beta = (a(j,j)-a(i,i))/(2.0*a(j,i))
             coeff = 0.5*beta/sqrt(1.0+beta**2)
             s = sqrt(max(0.5+coeff,0.0))
             c = sqrt(max(0.5-coeff,0.0))
       ! recalculate rows i and j
             do k=1,n
               cs =  c*a(i,k)+s*a(j,k)
               sc = -s*a(i,k)+c*a(j,k)
               a(i,k) = cs
               a(j,k) = sc
             end do
       ! new matrix a_{k+1} from a_{k}, and eigenvectors 
             do k=1,n
               cs =  c*a(k,i)+s*a(k,j)
               sc = -s*a(k,i)+c*a(k,j)
               a(k,i) = cs
               a(k,j) = sc
               cs =  c*x(k,i)+s*x(k,j)
               sc = -s*x(k,i)+c*x(k,j)
               x(k,i) = cs
               x(k,j) = sc
             end do
           end do
         end do
         IF (b2.LE.abserr) RETURN

       end do
       return
       end
       
       
       
!*******************************************************************
      FUNCTION ASSIGN_LVS(XX,YY,ZZ,LLVS,PI)
!*******************************************************************

      use flow_geom_var
      USE LVS_VAR
      real*8 r_tem
      
      assign_lvs=zz-0.0607!0.0609!0.0431!0.4998d0  !0.5002!
      
      ASSIGN_LVS=ASSIGN_LVS+1.E-12 !small disturbance for not aligning with grid.
      RETURN
      END      
      
!c*******************************************************************
      SUBROUTINE LVSINIT(U,V,W,PSI_CN,VOL_TOT_ORI)
!c*******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)

      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)
      REAL VOL_TOT_ORI(MLVS)

      INTEGER BUBCOL(MLVS)
      INTEGER BUBCOL2(MLVS)
      
      
      INTEGER I,J,K,N,NN,LLVS
      INTEGER LVS_START,LVS_END,NUMS,NUMS_INITIAL,NUMS_CRI
      REAL XX,YY,ZZ,ALPHI_TMP,DVOL,CRI
      INTEGER ISKIP,JSKIP,KSKIP
      REAL PI
      REAL ASSIGN_LVS
      ! REAL ALPHI_EXT2(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      
      ALLOCATE( ALPHI_EXT(-2:M1F+3,-2:M2F+3,-2:M3F+3,1))

      IF ( ITRACKING .EQ. 1 ) THEN
!c       DO LLVS=1,NLVS
!c        CALL BAND_GENERATION(LLVS)
!c        CALL ALPHI_BC(LLVS)
!C        CALL REINITIALIZATION(N_MAX-3,N_MAX-1,N_MAX,LLVS)
!C        CALL ALPHI_BC(LLVS)
!c       ENDDO
      DO LLVS=1,NLVS
      
      IF (LVSON(LLVS) .EQ. 0) GOTO 774 ! 2017-08-27
      
!$OMP PARALLEL DO private(I,J)
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       ALPHI_EXT(I,J,K,1)=0.
      ENDDO
      ENDDO
      ENDDO

       DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         ALPHI_EXT(I,J,K,1)=ALPHI(NN,N,LLVS)
        ENDDO
       ENDDO
       CALL ALPHI_BC(LLVS)

       CALL CAL_PSIF(1,N_MAX,LLVS,U,V,W,PSI_CN,BUBCOL,bubcol2
     &, vol_tot_ori)
     
       ! if (llvs. eq. 3) then
       ! write(*,*) 'first lvs written'      
       ! OPEN(146,FILE='N#alphi_contour_AFTER_BANDG.DAT')
       ! WRITE(146,*) 'VARIABLES="X","Y","Z","alphi"'
      ! WRITE(146,*) 'ZONE I=',N1FM/2,',J=',N2FM/2,',K=',N3FM/2,',F=POINT'
      ! DO K=1,N3FM,2
      ! DO J=1,N2FM,2
      ! DO I=1,N1FM,2
        ! WRITE(146,147) XPF(I),YPF(J),ZPF(K),alphi_ext(i,j,k,1)
      ! ENDDO
      ! ENDDO
      ! ENDDO          
       ! CLOSE(146)
 ! 147  FORMAT(4F15.8)       
       ! endif     
      write(*,*)'LVS&volume=',llvs,vol_tot_ori(llvs)
     
     
     
     
  774 CONTINUE !2017-08-27     
      ENDDO

       WRITE(*,*) 'LVS FIELD IS READ!!'

!**********************************************************************
!***********MAKE NEW LEVEL SET****************************************
      ELSE IF ( ITRACKING .EQ. 2 .or. ITRACKING .EQ. 12 ) THEN

!       IF (ITRACKING .EQ. 2) THEN
        LVS_START=1
        LVS_END=NLVS !MLVS 2017-10-08
!       ELSE
!
!       DO LLVS=1,NLVS
!        CALL CAL_PSIF(1,N_MAX,LLVS,U,V,W)
!       ENDDO
!         LVS_START=NLVS+1
!        LVS_END=MLVS
!       ENDIF 

      PI=ACOS(-1.)

      DO LLVS=LVS_START,LVS_END
      
      IF (LVSON(LLVS) .EQ. 0) GOTO 775 ! 2017-08-27
      
!$OMP PARALLEL DO private(I,J)
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       ALPHI_EXT(I,J,K,1)=0.
      ENDDO
      ENDDO
      ENDDO

!CCCCCC------------------------FIRST-STEP LVSINIT---------------------CCC
      IF (N3F .EQ. 2 ) THEN
       CRI=1.2*MAX(SDXF,SDYF)
      ELSE
       CRI=1.2*MAX(SDXF,SDYF,SDZF)  
      ENDIF

      !THIS IS FOR EQUALLY DISTRIBUTE THE INITIAL LEVEL-SET TO EACH NUMBER OF BAND.
      NUMS_INITIAL=0
!$OMP PARALLEL DO private(I,J,XX,YY,ZZ,ALPHI_TMP)
!$OMP&reduction(+:NUMS_INITIAL)
      DO 1000 K=1,N3FM
      DO 1000 J=1,N2FM
      DO 1000 I=1,N1FM 
        XX=XPF(I)
        YY=YPF(J)
        ZZ=ZPF(K)
        ALPHI_TMP=ASSIGN_LVS(XX,YY,ZZ,LLVS,PI)
        IF( ABS(ALPHI_TMP) .LE. CRI ) NUMS_INITIAL=NUMS_INITIAL+1
 1000  CONTINUE
          NUMS_CRI=NUMS_INITIAL/(N_MAX-1)  !-1 IS FOR SAFTY.

      NUMS_INITIAL=0
      N=1
C!$OMP PARALLEL DO private(I,J,XX,YY,ZZ,ALPHI_TMP)
C!$OMP&reduction(+:NUMS_INITIAL)
      DO 1001 K=1,N3FM
      DO 1001 J=1,N2FM
      DO 1001 I=1,N1FM 
        XX=XPF(I)
        YY=YPF(J)
        ZZ=ZPF(K)
        ALPHI_TMP=ASSIGN_LVS(XX,YY,ZZ,LLVS,PI)
        IF( ABS(ALPHI_TMP) .LE. CRI ) THEN
          NUMS_INITIAL=NUMS_INITIAL+1
          I_B(NUMS_INITIAL,N,LLVS)=I
          J_B(NUMS_INITIAL,N,LLVS)=J
          K_B(NUMS_INITIAL,N,LLVS)=K
          ALPHI_EXT(I,J,K,1)=ALPHI_TMP
          IF (NUMS_INITIAL .EQ. NUMS_CRI) THEN
            NUMA(N,LLVS)=NUMS_CRI
            NUMS_INITIAL=0
            N=N+1
          ENDIF
        ENDIF
 1001  CONTINUE
          NUMA(N,LLVS)=NUMS_INITIAL !REMAIN
!CCCCCC------------------------FIRST-STEP LVSINIT---------------------CCC

!CCCCCC----------------------------BAND_INIT--------------------------CCC
      CALL BAND_GENERATION(LLVS)
            ! DO N=1,11!(change)
        ! DO NN=1,NUMA(N,LLVS)
         ! I=I_B(NN,N,LLVS)
         ! J=J_B(NN,N,LLVS)
         ! K=K_B(NN,N,LLVS)
         ! XX=XPF(I)
         ! YY=YPF(J)
         ! ZZ=ZPF(K)
        ! alphi_ext2(i,j,k)=ALPHI_EXT(I,J,K,1)
        ! ENDDO
       ! ENDDO
       ! if ( llvs .eq. 3 ) then 
       ! OPEN(146,FILE='N#alphi_contour_AFTER_BANDG.DAT')
       ! WRITE(146,*) 'VARIABLES="X","Y","Z","alphi"'
      ! WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
      ! DO K=1,N3FM
      ! DO J=1,N2FM
      ! DO I=1,N1FM
        ! WRITE(146,147) XPF(I),YPF(J),ZPF(K),alphi_ext2(i,j,k)
      ! ENDDO
      ! ENDDO
      ! ENDDO          
       ! CLOSE(146)
 ! 147  FORMAT(4F15.8)
       ! endif
       DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K,XX,YY,ZZ)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         XX=XPF(I)
         YY=YPF(J)
         ZZ=ZPF(K)
         ALPHI_EXT(I,J,K,1)=ASSIGN_LVS(XX,YY,ZZ,LLVS,PI)
        ENDDO
       ENDDO

!!!!2017-08-19_ to see how alphi changes according to band#...
!     DO N=1,n_max
!       DO NN=1,NUMA(N,LLVS)
!        I=I_B(NN,N,LLVS)
!        J=J_B(NN,N,LLVS)
!        K=K_B(NN,N,LLVS)
!        XX=XPF(I)
!        YY=YPF(J)
!        ZZ=ZPF(K)
!       alphi_ext2(i,j,k)=ALPHI_EXT(I,J,K,1)
!       ENDDO
!      ENDDO
!      
!      if ( llvs .eq. NLVS ) then 
!      OPEN(146,FILE='N#alphi_contour_AFTER_BANDG.DAT')
!      WRITE(146,*) 'VARIABLES="X","Y","Z","alphi"'
!     WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!     DO K=1,N3FM
!     DO J=1,N2FM
!     DO I=1,N1FM
!       WRITE(146,147) XPF(I),YPF(J),ZPF(K),alphi_ext2(i,j,k)
!     ENDDO
!     ENDDO
!     ENDDO          
!      CLOSE(146)
!147  FORMAT(4F15.8)
!      endif
 
       ! 동일한 routine 이 band_generation 에 있음
       NUMA_MAX=0
       DO N=1,N_MAX
         NUMA_MAX=max(NUMA_MAX,NUMA(N,LLVS))
       ENDDO
        WRITE(*,*) 'NUMA_MAX/MF_BAND=',FLOAT(NUMA_MAX)/FLOAT(MF_BAND) ! MF_BAND 는 총 그리드 수의 어떤 PORTION이 되도록 설정
       IF (NUMA_MAX .GE. MF_BAND) THEN                                ! 활성화 된 CELL 의 개수가 위의 PORTION보다 크면 문제. 
        WRITE(*,*) 'MF_BAND IS NOT ENOUGH!!'
        WRITE(*,*) 'SIMULATION IS STOPPED!!'
        STOP
       ENDIF
!CCCCCC----------------------------BAND_INIT--------------------------CCC

       CALL CAL_PSIF(1,N_MAX,LLVS,U,V,W,PSI_CN,BUBCOL,bubcol2
     &, vol_tot_ori)
!C-----CALCULATE_ORIGINAL_TOTAL_MASS OF BUBBLE
        CALL CAL_PSIF_MTP(1,N_MAX,LLVS,PSIF_MTP)
       IF (N3FM .EQ. 1) THEN
        DVOL=SDXF*SDYF
       ELSE
        DVOL=SDXF*SDYF*SDZF
       ENDIF

        VOL_TOT_ORI(LLVS)=0.
C!$OMP PARALLEL DO private(I,J)
C!$OMP&reduction(+:VOL_TOT_ORI)
        DO K=1,N3FM
        DO J=1,N2FM
        DO I=1,N1FM
         VOL_TOT_ORI(LLVS)=VOL_TOT_ORI(LLVS)+(1.-PSIF_MTP(I,J,K))*DVOL
        ENDDO
        ENDDO
        ENDDO
      WRITE(*,*) 'TOTAL BUBBLE VOL = ',LLVS,VOL_TOT_ORI(LLVS)

       DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
          ALPHI(NN,N,LLVS)=ALPHI_EXT(I,J,K,1)
        ENDDO
       ENDDO

 775  CONTINUE    
       ENDDO  !DO L=1,NLVS

       WRITE(*,*) 'LEVEL-SET FIELD IS GENERATED!!!'

!C=====SAVE
       ISKIP=2
       JSKIP=2
       KSKIP=2
       IF (N3FM .EQ. 1) THEN
       OPEN(146,FILE='0PSIF_LVSINIT.DAT')
       WRITE(146,*) 'VARIABLES="X","Y","PSIF"'
      WRITE(146,*) 'ZONE I=',N1FM/ISKIP,',J=',N2FM/JSKIP,'F=POINT'
      DO J=1,N2FM,JSKIP
      DO I=1,N1FM,ISKIP
        WRITE(146,150) XPF(I),YPF(J),PSIF(i,j,1)
      ENDDO
      ENDDO
       CLOSE(146)
 150  FORMAT(3F15.8) 
       ELSE
       OPEN(146,FILE='0PSIF_LVSINIT.DAT')
       WRITE(146,*) 'VARIABLES="X","Y","Z","PSIF"'
      WRITE(146,*) 'ZONE I=',N1FM/ISKIP,',J=',N2FM/JSKIP,',K='
     &,N3FM/KSKIP,',F=POINT'
      DO K=1,N3FM,KSKIP
      DO J=1,N2FM,JSKIP
      DO I=1,N1FM,ISKIP
        WRITE(146,151) XPF(I),YPF(J),ZPF(K),PSIF(i,j,k)
      ENDDO
      ENDDO
      ENDDO
       CLOSE(146)
 151  FORMAT(4F15.8) 
       ENDIF !IF (N3FM .EQ. 1) THEN

       ENDIF
       
      

!C=====SAVE
       OPEN(146,FILE='0alphi_LVSINIT.DAT')
       WRITE(146,*) 'VARIABLES="X","Y","Z","alphi"'
      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
      DO K=1,N3FM
      DO J=1,N2FM
      DO I=1,N1FM
        WRITE(146,152) XPF(I),YPF(J),ZPF(K),alphi_ext(i,j,k,1)
      ENDDO
      ENDDO
      ENDDO
       CLOSE(146)
 152  FORMAT(4F15.8) 
!      stop
 
       DEALLOCATE(ALPHI_EXT)
       
       RETURN
       END
       
       
!*******************************************************************
       SUBROUTINE ALPHI_BC(LLVS)
!******************************************************************* 
      USE PARAM_VAR

      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR

      IMPLICIT NONE
      
      INTEGER LLVS
      
      INTEGER I,J,K

      !X-DIR
      IF (IPX .NE. 1) THEN      
!$OMP PARALLEL DO private(J)
      DO K=1,N3F
      DO J=1,N2F
       IF (ABS(ALPHI_EXT(1,J,K,1)) .LT. 1.E-8) THEN
        ALPHI_EXT(0,J,K,1)=0.
        ALPHI_EXT(-1,J,K,1)=0.
        ALPHI_EXT(-2,J,K,1)=0.
       ELSE IF (ALPHI_EXT(1,J,K,1) .GT. 0.) THEN
        ALPHI_EXT(0,J,K,1)=ALPHI_EXT(1,J,K,1)+SDXF
        ALPHI_EXT(-1,J,K,1)=ALPHI_EXT(0,J,K,1)+SDXF
        ALPHI_EXT(-2,J,K,1)=ALPHI_EXT(-1,J,K,1)+SDXF
       ELSE IF (ALPHI_EXT(1,J,K,1) .LT. 0.) THEN
        ALPHI_EXT(0,J,K,1)=ALPHI_EXT(1,J,K,1)-SDXF
        ALPHI_EXT(-1,J,K,1)=ALPHI_EXT(0,J,K,1)-SDXF
        ALPHI_EXT(-2,J,K,1)=ALPHI_EXT(-1,J,K,1)-SDXF
       ENDIF

       IF (ABS(ALPHI_EXT(N1FM,J,K,1)) .LT. 1.E-8) THEN
       ALPHI_EXT(N1F,J,K,1)=0.
       ALPHI_EXT(N1F+1,J,K,1)=0.
       ALPHI_EXT(N1F+2,J,K,1)=0.
       ALPHI_EXT(N1F+3,J,K,1)=0.
       ELSE IF (ALPHI_EXT(N1FM,J,K,1) .GT. 0.) THEN
       ALPHI_EXT(N1F,J,K,1)=ALPHI_EXT(N1FM,J,K,1)+SDXF
       ALPHI_EXT(N1F+1,J,K,1)=ALPHI_EXT(N1F,J,K,1)+SDXF
       ALPHI_EXT(N1F+2,J,K,1)=ALPHI_EXT(N1F+1,J,K,1)+SDXF
       ALPHI_EXT(N1F+3,J,K,1)=ALPHI_EXT(N1F+2,J,K,1)+SDXF
       ELSE IF (ALPHI_EXT(N1FM,J,K,1) .LT. 0.) THEN
       ALPHI_EXT(N1F,J,K,1)=ALPHI_EXT(N1FM,J,K,1)-SDXF
       ALPHI_EXT(N1F+1,J,K,1)=ALPHI_EXT(N1F,J,K,1)-SDXF
       ALPHI_EXT(N1F+2,J,K,1)=ALPHI_EXT(N1F+1,J,K,1)-SDXF
       ALPHI_EXT(N1F+3,J,K,1)=ALPHI_EXT(N1F+2,J,K,1)-SDXF
       ENDIF
      ENDDO
      ENDDO
      ENDIF

      !Y-DIR
      IF (IPY .NE. 1) THEN      
!$OMP PARALLEL DO private(I)
      DO K=1,N3F
      DO I=1,N1F
       IF (ABS(ALPHI_EXT(I,1,K,1)) .LT. 1.E-8) THEN
       ALPHI_EXT(I,0,K,1)=0.
       ALPHI_EXT(I,-1,K,1)=0.
       ALPHI_EXT(I,-2,K,1)=0.
       ELSE IF (ALPHI_EXT(I,1,K,1) .GT. 0.) THEN
       ALPHI_EXT(I,0,K,1)=ALPHI_EXT(I,1,K,1)+SDYF
       ALPHI_EXT(I,-1,K,1)=ALPHI_EXT(I,0,K,1)+SDYF
       ALPHI_EXT(I,-2,K,1)=ALPHI_EXT(I,-1,K,1)+SDYF
       ELSE IF (ALPHI_EXT(I,1,K,1) .LT. 0.) THEN
       ALPHI_EXT(I,0,K,1)=ALPHI_EXT(I,1,K,1)-SDYF
       ALPHI_EXT(I,-1,K,1)=ALPHI_EXT(I,0,K,1)-SDYF
       ALPHI_EXT(I,-2,K,1)=ALPHI_EXT(I,-1,K,1)-SDYF
       ENDIF

       IF (ABS(ALPHI_EXT(I,N2FM,K,1)) .LT. 1.E-8) THEN
       ALPHI_EXT(I,N2F,K,1)=0.
       ALPHI_EXT(I,N2F+1,K,1)=0.
       ALPHI_EXT(I,N2F+2,K,1)=0.
       ALPHI_EXT(I,N2F+3,K,1)=0.
       ELSE IF (ALPHI_EXT(I,N2FM,K,1) .GT. 0.) THEN
       ALPHI_EXT(I,N2F,K,1)=ALPHI_EXT(I,N2FM,K,1)+SDYF
       ALPHI_EXT(I,N2F+1,K,1)=ALPHI_EXT(I,N2F,K,1)+SDYF
       ALPHI_EXT(I,N2F+2,K,1)=ALPHI_EXT(I,N2F+1,K,1)+SDYF
       ALPHI_EXT(I,N2F+3,K,1)=ALPHI_EXT(I,N2F+2,K,1)+SDYF
       ELSE IF (ALPHI_EXT(I,N2FM,K,1) .LT. 0.) THEN
       ALPHI_EXT(I,N2F,K,1)=ALPHI_EXT(I,N2FM,K,1)-SDYF
       ALPHI_EXT(I,N2F+1,K,1)=ALPHI_EXT(I,N2F,K,1)-SDYF
       ALPHI_EXT(I,N2F+2,K,1)=ALPHI_EXT(I,N2F+1,K,1)-SDYF
       ALPHI_EXT(I,N2F+3,K,1)=ALPHI_EXT(I,N2F+2,K,1)-SDYF
       ENDIF
      ENDDO
      ENDDO
      ENDIF

      !Z-DIR
      IF (IPZ .NE. 1) THEN      
!$OMP PARALLEL DO private(I)
      DO J=1,N2F
      DO I=1,N1F
       IF (ABS(ALPHI_EXT(I,J,1,1)) .LT. 1.E-8) THEN
       ALPHI_EXT(I,J,0,1)=0.
       ALPHI_EXT(I,J,-1,1)=0.
       ALPHI_EXT(I,J,-2,1)=0.
       ELSE IF (ALPHI_EXT(I,J,1,1) .GT. 0.) THEN
       ALPHI_EXT(I,J,0,1)=ALPHI_EXT(I,J,1,1)+SDZF
       ALPHI_EXT(I,J,-1,1)=ALPHI_EXT(I,J,0,1)+SDZF
       ALPHI_EXT(I,J,-2,1)=ALPHI_EXT(I,J,-1,1)+SDZF
       ELSE IF (ALPHI_EXT(I,J,1,1) .LT. 0.) THEN
       ALPHI_EXT(I,J,0,1)=ALPHI_EXT(I,J,1,1)-SDZF
       ALPHI_EXT(I,J,-1,1)=ALPHI_EXT(I,J,0,1)-SDZF
       ALPHI_EXT(I,J,-2,1)=ALPHI_EXT(I,J,-1,1)-SDZF
       ENDIF

       IF (ABS(ALPHI_EXT(I,J,N3FM,1)) .LT. 1.E-8) THEN
       ALPHI_EXT(I,J,N3F,1)=0.
       ALPHI_EXT(I,J,N3F+1,1)=0.
       ALPHI_EXT(I,J,N3F+2,1)=0.
       ALPHI_EXT(I,J,N3F+3,1)=0.
       ELSE IF (ALPHI_EXT(I,J,N3FM,1) .GT. 0.) THEN
       ALPHI_EXT(I,J,N3F,1)=ALPHI_EXT(I,J,N3FM,1)+SDZF
       ALPHI_EXT(I,J,N3F+1,1)=ALPHI_EXT(I,J,N3F,1)+SDZF
       ALPHI_EXT(I,J,N3F+2,1)=ALPHI_EXT(I,J,N3F+1,1)+SDZF
       ALPHI_EXT(I,J,N3F+3,1)=ALPHI_EXT(I,J,N3F+2,1)+SDZF
       ELSE IF (ALPHI_EXT(I,J,N3FM,1) .LT. 0.) THEN
       ALPHI_EXT(I,J,N3F,1)=ALPHI_EXT(I,J,N3FM,1)-SDZF
       ALPHI_EXT(I,J,N3F+1,1)=ALPHI_EXT(I,J,N3F,1)-SDZF
       ALPHI_EXT(I,J,N3F+2,1)=ALPHI_EXT(I,J,N3F+1,1)-SDZF
       ALPHI_EXT(I,J,N3F+3,1)=ALPHI_EXT(I,J,N3F+2,1)-SDZF
       ENDIF

      ENDDO
      ENDDO
      ENDIF

      
!Periodic Condition_ADDED 2018-08-31-- CAN IGNORE WHEN DOING BUBBLY FLOW  
      !X-DIR
      IF (IPX .EQ. 1) THEN      
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F    !START FROM ZERO BECAUSE OF EDGE OF DOMAIN.(NEED FOR INTER_PSI)
      DO J=0,N2F
      DO I=0,N1F_BD
        ALPHI_EXT(-I,J,K,1)=ALPHI_EXT(N1FM-I,J,K,1)
        ALPHI_EXT(N1F+I,J,K,1)=ALPHI_EXT(1+I,J,K,1)
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      !Y-DIR
      IF (IPY .NE. 1) THEN      
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F
      DO J=0,N2F_BD
      DO I=0,N1F
        ALPHI_EXT(I,-J,K,1)=ALPHI_EXT(I,N2FM-J,K,1)
        ALPHI_EXT(I,N2F+J,K,1)=ALPHI_EXT(I,1+J,K,1)
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      !Z-DIR
      IF (IPZ .NE. 1) THEN      
      DO K=0,N3F_BD
!$OMP PARALLEL DO private(I)
      DO J=0,N2F
      DO I=0,N1F
        ALPHI_EXT(I,J,-K,1)=ALPHI_EXT(I,J,N3FM-K,1)
        ALPHI_EXT(I,J,N3F+K,1)=ALPHI_EXT(I,J,1+K,1)
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      
      
      
      
      
       RETURN
       END 


!C*******************************************************************
      SUBROUTINE REIN_TRIGGER(NN1,LLVS,DALPHIMAX,DALPHIMIN)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE

      INTEGER NN1,LLVS
      REAL DALPHIMAX,DALPHIMIN
      
      INTEGER I,J,K,N,NN
      REAL DALPHIMAX_OMP,DALPHIMIN_OMP
      REAL DALPHIX,DALPHIY,DALPHIZ,DALPHI_SUM

        DALPHIMAX=0.
        DALPHIMIN=1.

       DO N=1,NN1                ! T-BAND ONLY (HERRMANN 2008)
        DALPHIMAX_OMP=0.
        DALPHIMIN_OMP=1.
!$OMP PARALLEL DO private(I,J,K,DALPHIX,DALPHIY,DALPHIZ,DALPHI_SUM)
!$OMP&reduction(MAX:DALPHIMAX_OMP)
!$OMP&reduction(MIN:DALPHIMIN_OMP)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
        DALPHIX=0.5*( ALPHI_EXT(IPF(I),J,K,1)
     &                               -ALPHI_EXT(IMF(I),J,K,1) )*SSDXF
        DALPHIY=0.5*( ALPHI_EXT(I,JPF(J),K,1)
     &                               -ALPHI_EXT(I,JMF(J),K,1) )*SSDYF
        DALPHIZ=0.5*( ALPHI_EXT(I,J,KPF(K),1)
     &                               -ALPHI_EXT(I,J,KMF(K),1) )*SSDZF
        DALPHI_SUM=DALPHIX**2+DALPHIY**2+DALPHIZ**2 

        DALPHIMAX_OMP=MAX(DALPHI_SUM,DALPHIMAX_OMP)
        DALPHIMIN_OMP=MIN(DALPHI_SUM,DALPHIMIN_OMP)
       ENDDO
         DALPHIMAX=MAX(DALPHIMAX,DALPHIMAX_OMP)
         DALPHIMIN=MIN(DALPHIMIN,DALPHIMIN_OMP)
       ENDDO
         DALPHIMAX=SQRT(DALPHIMAX)
         DALPHIMIN=SQRT(DALPHIMIN)

!         WRITE(*,900) DALPHIMAX,DALPHIMIN
! 900   FORMAT('DPHIMAX=',F10.5,'  DPHIMIN=',F10.5)

      RETURN
      END    


!C*******************************************************************
      SUBROUTINE TRANSPORT(U,V,W,PSI_CN,VOL_TOT_ORI)
!C*******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      use phase_boiling
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)

      REAL VOL_TOT_ORI(MLVS)

      INTEGER BUBCOL(MLVS)
      INTEGER BUBCOL2(MLVS)
      INTEGER L,LL,L2  !2017-08-16
      REAL TCON
      INTEGER LCOUNT              
      
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: ALPHI_RK3
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: VOF
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: UF,VF,WF
      REAL, DIMENSION (:,:,:), ALLOCATABLE :: UT,VT,WT
      
       
      INTEGER NN1,NN2,I,J,K,N,NN,LLVS,ISUBCYCLE,IMASS
      REAL DTLVS,TIME_TRAN,CFLM,CFLM_OMP,CFLX,CFLY,CFLZ,CFLI,UUVVWW
      REAL DALPHIMAX,DALPHIMIN
      
      real*8 dvol
      
      ALLOCATE( ALPHI_EXT(-2:M1F+3,-2:M2F+3,-2:M3F+3,1))

      ALLOCATE(UT(-1:M1+2,-1:M2+2,-1:M3+2),VT(-1:M1+2,-1:M2+2,-1:M3+2),
     &                                     WT(-1:M1+2,-1:M2+2,-1:M3+2))
      

       NN1=N_MAX-3
       NN2=N_MAX-1

            
      !여기서 TRANSPORT 전에 처리하자.
      !ALL TO 1 LVS
      !CALL LVSALLTOONE
      
       DO 776 LLVS=1,NLVS !MULTIPLE_LEVEL_SET_APPROACH
      IF (LVSON(LLVS) .EQ. 0) GOTO 776 ! 2017-08-28

!$OMP PARALLEL DO private(I,J,K)
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       ALPHI_EXT(I,J,K,1)=0.
      ENDDO
      ENDDO
      ENDDO

       DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         ALPHI_EXT(I,J,K,1)=ALPHI(NN,N,LLVS)
        ENDDO
       ENDDO
       CALL ALPHI_BC(LLVS)

      ALLOCATE(UF(M1F,M2F,M3F),VF(M1F,M2F,M3F),WF(M1F,M2F,M3F))
      ALLOCATE(ALPHI_RK3(-2:M1F+3,-2:M2F+3,-2:M3F+3))
      
        IF (LLVS .EQ. 1) THEN
          CALL VEL_INTP(2,NN2,LLVS,UF,VF,WF,U,V,W,UT,VT,WT)
        ELSE
          CALL VEL_INTP(0,NN2,LLVS,UF,VF,WF,U,V,W,UT,VT,WT)
        ENDIF
        
        !2017-12-22
        ! ! lag_prt transport temp
        ! call lagtrans(u,v,w)      
        
        
      !SUBCYCLING : Simulating Two-Phase Flows Using the 
      !                    refned Level Set Grid Method - HERRMANN(2006)
        TIME_TRAN=0.
      DO ISUBCYCLE=1,10    !SELECTED discretionally

!C-----CAL. TRANSPORT DT
      CFLM=0.
       DO N=1,NN1
      CFLM_OMP=0.
!$OMP PARALLEL DO private(I,J,K,CFLX,CFLY,CFLZ,CFLI)
!$OMP&reduction(MAX:CFLM_OMP)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         CFLX=ABS(0.5*(UF(IMF(I),J,K)+UF(I,J,K)))*SSDXF
         CFLY=ABS(0.5*(VF(I,JMF(J),K)+VF(I,J,K)))*SSDYF
         CFLZ=ABS(0.5*(WF(I,J,KMF(K))+WF(I,J,K)))*SSDZF
         CFLI=CFLX+CFLY+CFLZ
         CFLM_OMP=MAX(CFLM_OMP,CFLI)
        ENDDO
         CFLM=MAX(CFLM,CFLM_OMP)
       ENDDO

      IF ( CFLM .EQ. 0. ) THEN
        UUVVWW=0.
!$OMP PARALLEL DO private(NN,I,J,K)
!$OMP&reduction(+:UUVVWW)
         DO N=1,NN1
           DO NN=1,NUMA(N,LLVS)
            I=I_B(NN,N,LLVS)
            J=J_B(NN,N,LLVS)
            K=K_B(NN,N,LLVS)
              UUVVWW=UUVVWW+(UF(I,J,K)+VF(I,J,K)+WF(I,J,K))
         ENDDO
        ENDDO        
        
         IF ( UUVVWW .EQ. 0. ) THEN 
             WRITE(*,*) 'ALL VELOCITY IS ZERO!!', llvs
             if (mlvs .gt. 1) lvson(llvs)=0
             if (iphs .ne. 1) GOTO 3000 !SKIP THE TRANSPORT
             dtlvs=dtconst/10d0
         Else
         DTLVS=1./CFLM
         endif
      else
      DTLVS=1./CFLM
      ENDIF



       IF( TIME_TRAN+DTLVS .GE. DTCONST) DTLVS=DTCONST-TIME_TRAN
        TIME_TRAN=TIME_TRAN+DTLVS

!C-----3RD_ORDER TVD RK
       DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)    
         ALPHI_RK3(I,J,K)=ALPHI_EXT(I,J,K,1)
        ENDDO
       ENDDO


      CALL TRANSPORT_RK3(DTLVS,NN1,LLVS,UF,VF,WF) 
      CALL TRANSPORT_RK3(DTLVS,NN1,LLVS,UF,VF,WF)
      
       DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)             
         ALPHI_EXT(I,J,K,1)=0.75*ALPHI_RK3(I,J,K)
     &                                      +0.25*ALPHI_EXT(I,J,K,1)  
         ENDDO
        ENDDO

      CALL TRANSPORT_RK3(DTLVS,NN1,LLVS,UF,VF,WF)

       DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)                
         ALPHI_EXT(I,J,K,1)=1./3.*ALPHI_RK3(I,J,K)
     &                                     +2./3.*ALPHI_EXT(I,J,K,1)
         ENDDO
        ENDDO

       !IBM & LVS INTERACTION
!       CALL WALL_INTERACTION(NN1,LLVS) !AFTER SUBROUTIN GLOBAL_MASS_CORRECTION
!-------------------------------------1dpart----------------------------
!$OMP PARALLEL DO private(J)
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
      ALPHI_EXT(I,J,K,1)=ALPHI_EXT(2,2,K,1)
      ENDDO
      ENDDO
      ENDDO
!-----------------------------------------------------------------------
    
!c-----RE-INITIALIZATION_TRIGGER_START
       CALL REIN_TRIGGER(NN1,LLVS,DALPHIMAX,DALPHIMIN)
  
       IF( (DALPHIMAX .GT. 2.) .OR. (DALPHIMIN .LT. 1.E-4) ) THEN
         IMASS=1
       ELSE
         IMASS=0
       ENDIF
       
       
       IF( BUBCOL(LLVS) .EQ. 1 .or. iphs .eq. 1) THEN
       call reinitialization(nn1,nn2,nn2,llvs)
       !imass=1
       IF (MOD(M,50) .EQ. 0) IMASS=1  !EVERY 50-STEP
       ENDIF

       !RE-INITIALIZATION CRITERIA (HERRMANN 2008)
       IF( (DALPHIMAX .GT. 2.) .OR. (DALPHIMIN .LT.  1.E-4)
     &                               .OR. IMASS .EQ. 1 ) THEN  ! BUBCOL(LLVS) .EQ. 1 ) THEN  
       !WHEN BUB_COLLISION, DO THE REINITIALIZATION EXCEPT MASS_CORRECTION.
        WRITE(*,4000) LLVS,IMASS
 4000  FORMAT('!!RE-INITIALIZATION, LLVS=',I5,',IMASS=',I5)

          IF (MGLOBAL .EQ. 1 ) THEN !GLOBAL MASS CORRECTION
          IF (IMASS .EQ. 1) THEN
            write(*,*) 'at transport'
          	CALL GLOBAL_MASS_CORRECTION(NN1,LLVS,VOL_TOT_ORI,dtlvs)
          ENDIF
          ENDIF

        CALL ALPHI_BC(LLVS)
        CALL BAND_GENERATION(LLVS)
        CALL ALPHI_BC(LLVS)
        CALL REINITIALIZATION(NN1,NN2,NN2,LLVS)
       ENDIF   !RE-INITIALIZATION CRITERIA (HERRMANN 2008)
        CALL ALPHI_BC(LLVS)

!        write(*,*) 'solve reinitialization step per every step!'
!c-----RE-INITIALIZATION_TRIGGER_END
!              WRITE(*,*) 'NO REIN_TRIGGER!!!!!!!!!!!!'
!C-------------------------MASS_CONSERVING_LEVEL_SET--------------------C
! !LEVEL-SET CORRECTION WITH N+1 LEVEL-SET FUNC
       ! IF(MCLS .EQ. 1) THEN
         ! CALL CAL_MCLS(2,NN1,MGLOBAL,LLVS,DTLVS,VOF,UF,VF,WF)
      ! DEALLOCATE(VOF)
       ! ENDIF
!C-------------------------MASS_CONSERVING_LEVEL_SET--------------------C

       IF ( TIME_TRAN .EQ. DTCONST ) GOTO 1000
       IF( (DALPHIMAX .GT. 2.) .OR. (DALPHIMIN .LT. 1.E-4) ) THEN  
        CALL VEL_INTP(0,NN2,LLVS,UF,VF,WF,U,V,W,UT,VT,WT) !VELOCITY FOR NEW BAND
       ENDIF

      ENDDO

        WRITE(*,*) 'TRANSPORT SUBCYCLE ITERATION IS EXCEEDED!!',LLVS
        !STOP

 1000  CONTINUE

 3000  CONTINUE
 
!-----before deallocating uf,vf,wf... we have to assess whether this 
!-----level set will need more resolution.
!      IF (MSUB.EQ.3) CALL THIN_CRI(UF,VF,wf,LLVS)
 
      DEALLOCATE(UF,VF,WF)
      DEALLOCATE(ALPHI_RK3)
       !write(*,*) 'working on lvs=',llvs ! 2017-10-06, comment after check
  
       CALL CAL_PSIF(0,NN1,LLVS,U,V,W,PSI_CN,BUBCOL,bubcol2
     &,vol_tot_ori)
     
       DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
          ALPHI(NN,N,LLVS)=ALPHI_EXT(I,J,K,1)
        ENDDO
       ENDDO
       
       
      if (imass .eq. 0) then 
      IF (IPHS .EQ. 1) THEN 
      write(*,*) 'after transport'
      call global_mass_correction(NN1,LLVS,VOL_TOT_ORI,dtconst) 
      ENDIF  
      endif             
       
  776 CONTINUE

  
  
!C=====SAVE
!       OPEN(146,FILE='0alphi_transport2.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","alphi"'
!      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!      DO K=1,N3FM
!      DO J=1,N2FM
!      DO I=1,N1FM
!        WRITE(146,152) XPF(I),YPF(J),ZPF(K),alphi_ext(i,j,k,1)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
! 152  FORMAT(4F15.8) 
!      if (m .eq. 10) stop  

  
  
!        !2017-12-22
!        ! ! lag_prt transport temp
!        ! 2018-08-02_call lagtrans(u,v,w,VOL_TOT_ORI) 
!        ! 2018-08-02_CALL LAGLVS(U,V,W,VOL_TOT_ORI) ! MOVED TO HERE DUE TO VOL_TOT_ORI 2018-01-18
!       !CONTACT TIME COUNTER 2017-08-16---------------------------------            
!       DO L=2,NLVS
!       L2=L-1
!       DO LL=1,L2
!       
!      IF ((CONTACT_OLD(L,LL) .EQ. 0).AND.(CONTACT(L,LL) .EQ. 0)) THEN
!       CONTACT_TIME(L,LL)=0.
!       CONTACT_START(L,LL)=0.
!       !write(*,*) 'no contact'
!       
!      ELSEIF ((CONTACT_OLD(L,LL) .EQ. 0).AND.(CONTACT(L,LL) .EQ. 1))THEN
!       CONTACT_START(L,LL)=TIME        
!       CONTACT_TIME(L,LL)=TIME-CONTACT_START(L,LL)   
!       CALL CAL_TCON(L,LL,vol_tot_ori)
!       WRITE(*,800) L,LL,TIME,T_CON(L,LL)
!       
!      ELSEIF ((CONTACT_OLD(L,LL) .EQ. 1).AND.(CONTACT(L,LL) .EQ. 0))THEN
!       CONTACT_TIME(L,LL)=TIME-CONTACT_START(L,LL) ! CONTACT 끗
!       CONTACT_START(L,LL)=0.
!       WRITE(*,801) L,LL,TIME 
!       WRITE(*,802) L,LL,CONTACT_TIME(L,LL)    
!       CONTACT_TIME(L,LL)=0. ! INITIALIZE CONTACT_TIME       
!      ELSEIF ((CONTACT_OLD(L,LL) .EQ. 1).AND.(CONTACT(L,LL) .EQ. 1))THEN 
!      CONTACT_TIME(L,LL)=TIME-CONTACT_START(L,LL)  
!      WRITE(*,803) L,LL,CONTACT_START(L,LL), TIME  
!      !CONTACT_START(L,LL)=CONTACT_START(LL)
!      ENDIF
!
!      CONTACT_OLD(L,LL)=CONTACT(L,LL) ! 2017-08-18
!      CONTACT_OLD(LL,L)=CONTACT(LL,L) ! 2017-08-20
!      CONTACT(LL,L)=0                 ! 2017-08-20
!      CONTACT(l,ll)=0                 ! 2017-08-18
!       IF ((CONTACT_TIME(L,LL).GE.T_CON(L,LL)).and.(msub .eq. 3)) then ! no contact_time_delay
!       ! WRITE(*,*) 'No contact time delay'
!       CALL CAL_COALESCENCE(L,LL,PSI_CN,VOL_TOT_ORI)       ! LL is smaller lvs# then L
!       write(*,*) 'calculated coalescence'
!       CONTACT_TIME(L,LL)=0.
!       CONTACT_START(L,LL)=0.
!       t_con(l,ll)=100.                           
!       endif
!       ENDDO
!       ENDDO              
!       !CONTACT TIME COUNTER 2017-08-16-------------------------------- 
! 800  FORMAT('-',I5,' AND',I5,' START CONTACT AT', F12.5,' TCON=',F12.5)
! 801  FORMAT('-',I5,' AND',I5,' END   CONTACT AT', F12.5)
! 802  FORMAT('-',I5,' AND',I5,' END   CONTACT, dt_C=', F12.5)     
! 803  FORMAT('-',I5,' AND',I5,' KEEP  CONTACT FROM', F12.5,' TO',F12.5)     
!       !CONTACT TIME COUNTER 2017-08-20--------------------------------      
       
!         LCOUNT=0
!       DO I=1,MLVS
!       IF (LVSON(I) .EQ. 1) THEN
!       LCOUNT=LCOUNT+1
!       ENDIF
!       ENDDO
!       WRITE(*,*) 'NUMBER OF FINAL LVSs,nlvs', LCOUNT, nlvs
        
        
        ! if (msub.eq. 3) then
        ! OPEN(146,FILE='0alp_contourA3.DAT')
        ! WRITE(146,*) 'VARIABLES="X","Y","Z","alphi"'
       ! WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
       ! DO K=1,N3FM
       ! DO J=1,N2FM
       ! DO I=1,N1FM
         ! WRITE(146,147) XPF(I),YPF(J),ZPF(K),alphi_ext(i,j,k,1)
       ! ENDDO
       ! ENDDO
       ! ENDDO
        ! CLOSE(146)
       ! endif
        DEALLOCATE(UT,VT,WT)
        DEALLOCATE(ALPHI_EXT)

        ! IF (IPHS .EQ. 1) DEALLOCATE(VOL_TOT_IPHS)
!       if (msub .eq. 1) then
!       OPEN(146,FILE='0psif_contourA1.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","psif","vof","alphi"'
!      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!      DO K=1,N3FM
!      DO J=1,N2FM
!      DO I=1,N1FM
!        WRITE(146,147) XPF(I),YPF(J),ZPF(K),psif(i,j,k),vof(i,j,k)
!     &,alphi(i,j,k)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
!       else if (msub .eq. 2) then
!      OPEN(146,FILE='0psif_contourA2.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","psif","vof","alphi"'
!      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!      DO K=1,N3FM
!      DO J=1,N2FM
!      DO I=1,N1FM
!        WRITE(146,147) XPF(I),YPF(J),ZPF(K),psif(i,j,k),vof(i,j,k)
!     &,alphi(i,j,k)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
!       else if (msub.eq. 3) then
!       OPEN(146,FILE='0psif_contourA3.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","psif","vof","alphi"'
!      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!      DO K=1,N3FM
!      DO J=1,N2FM
!      DO I=1,N1FM
!        WRITE(146,147) XPF(I),YPF(J),ZPF(K),psif(i,j,k),vof(i,j,k)
!     &,alphi(i,j,k)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
!      endif
!
 147  FORMAT(4F15.8) 

      RETURN
      END   

!C*******************************************************************
      SUBROUTINE TRANSPORT_RK3(DTLVS,NN1,LLVS,UF,VF,WF)
!C*******************************************************************
!C=====5th-WENO
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      !BOILING
      USE TWO_PHASE_PROPERTY
      USE PHASE_BOILING
      IMPLICIT NONE
      
       REAL ALPHI_NEW(MF_BAND,M_MAX)
       REAL DPHI_X(-2:M1F+3,-2:M2F+3,-2:M3F+3)
     &      ,DPHI_Y(-2:M1F+3,-2:M2F+3,-2:M3F+3)
     &      ,DPHI_Z(-2:M1F+3,-2:M2F+3,-2:M3F+3)
       REAL UF(M1F,M2F,M3F),VF(M1F,M2F,M3F),WF(M1F,M2F,M3F) 
       
      INTEGER NN1,LLVS
      REAL DTLVS
      
      INTEGER N,NN,I,J,K,IP1,IP2,IM1,IM2,IM3
      INTEGER JP1,JP2,JM1,JM2,JM3,KP1,KP2,KM1,KM2,KM3
      REAL CRI,UU1,UU2,UU3,V1,V2,V3,V4,V5
      REAL WENO_X,WENO_Y,WENO_Z,UX_MAX,UY_MAX,UZ_MAX
      REAL ADLPHIXM,ADLPHIXP,ALPHI_X,ADLPHIYM,ADLPHIYP,ALPHI_Y
      REAL ADLPHIZM,ADLPHIZP,ALPHI_Z
      REAL FUNCBODY

!C======================TO AVOID REPEATED CAL._START====================C
      !DPHI_X,DPHI_Y,DPHI_Z IS DEFINED AT CELL-INTERFACE
       DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         DPHI_X(I,J,K)=(ALPHI_EXT(IPF(I),J,K,1)
     &                                    -ALPHI_EXT(I,J,K,1))*SSDXF
         DPHI_Y(I,J,K)=(ALPHI_EXT(I,JPF(J),K,1)
     &                                    -ALPHI_EXT(I,J,K,1))*SSDYF        
         DPHI_Z(I,J,K)=(ALPHI_EXT(I,J,KPF(K),1)
     &                                    -ALPHI_EXT(I,J,K,1))*SSDZF
         ENDDO
       ENDDO

      N=NN1
!$OMP PARALLEL DO private(I,J,K,IP1,IP2,IM1,IM2,IM3)
!$OMP&private(JP1,JP2,JM1,JM2,JM3,KP1,KP2,KM1,KM2,KM3)
       DO NN=1,NUMA(N,LLVS)
        I=I_B(NN,N,LLVS)
        J=J_B(NN,N,LLVS)
        K=K_B(NN,N,LLVS)
        IM1=IMF(I)        !!IF INTERFACE TOUCHES THE BOUNDARY, 
        IM2=IMF(IM1)        !!BOUNDARY CONDITION IS NEEDED HERE!!
        IM3=IMF(IM2)
        IP1=IPF(I)
        IP2=IPF(IP1)
      DPHI_X(IM3,J,K)=(ALPHI_EXT(IM2,J,K,1)
     &                                   -ALPHI_EXT(IM3,J,K,1))*SSDXF            
      DPHI_X(IM2,J,K)=(ALPHI_EXT(IM1,J,K,1)
     &                                   -ALPHI_EXT(IM2,J,K,1))*SSDXF
      DPHI_X(IM1,J,K)=(ALPHI_EXT(I,J,K,1)
     &                                   -ALPHI_EXT(IM1,J,K,1))*SSDXF
      DPHI_X(IP1,J,K)=(ALPHI_EXT(IP2,J,K,1)
     &                                   -ALPHI_EXT(IP1,J,K,1))*SSDXF
      DPHI_X(IP2,J,K)=(ALPHI_EXT(IPF(IP2),J,K,1)
     &                                   -ALPHI_EXT(IP2,J,K,1))*SSDXF
        JM1=JMF(J)
        JM2=JMF(JM1)
        JM3=JMF(JM2)
        JP1=JPF(J)
        JP2=JPF(JP1)
      DPHI_Y(I,JM3,K)=(ALPHI_EXT(I,JM2,K,1)
     &                                   -ALPHI_EXT(I,JM3,K,1))*SSDYF
      DPHI_Y(I,JM2,K)=(ALPHI_EXT(I,JM1,K,1)
     &                                   -ALPHI_EXT(I,JM2,K,1))*SSDYF
      DPHI_Y(I,JM1,K)=(ALPHI_EXT(I,J,K,1)
     &                                   -ALPHI_EXT(I,JM1,K,1))*SSDYF
      DPHI_Y(I,JP1,K)=(ALPHI_EXT(I,JP2,K,1)
     &                                   -ALPHI_EXT(I,JP1,K,1))*SSDYF
      DPHI_Y(I,JP2,K)=(ALPHI_EXT(I,JPF(JP2),K,1)
     &                                   -ALPHI_EXT(I,JP2,K,1))*SSDYF
        KM1=KMF(K)
        KM2=KMF(KM1)
        KM3=KMF(KM2)
        KP1=KPF(K)
        KP2=KPF(KP1)
      DPHI_Z(I,J,KM3)=(ALPHI_EXT(I,J,KM2,1)
     &                                   -ALPHI_EXT(I,J,KM3,1))*SSDZF
      DPHI_Z(I,J,KM2)=(ALPHI_EXT(I,J,KM1,1)
     &                                   -ALPHI_EXT(I,J,KM2,1))*SSDZF
      DPHI_Z(I,J,KM1)=(ALPHI_EXT(I,J,K,1)
     &                                   -ALPHI_EXT(I,J,KM1,1))*SSDZF
      DPHI_Z(I,J,KP1)=(ALPHI_EXT(I,J,KP2,1)
     &                                   -ALPHI_EXT(I,J,KP1,1))*SSDZF
      DPHI_Z(I,J,KP2)=(ALPHI_EXT(I,J,KPF(KP2),1)
     &                                   -ALPHI_EXT(I,J,KP2,1))*SSDZF
       ENDDO
!C======================TO AVOID REPEATED CAL._END======================C

!C=======================SOLVING_TRANSPORT_EQN_START====================C

       CRI=2.*(SDXF*SDYF*SDZF)**(1./3.)

      !Roe scheme with local lax-friedrich entropy correction
       DO N=1,NN1
!$OMP PARALLEL DO 
!$OMP&private(I,J,K,UU1,UU2,UU3,V1,V2,V3,V4,V5,IP1,IP2,IM1,IM2,IM3)
!$OMP&private(JP1,JP2,JM1,JM2,JM3,KP1,KP2,KM1,KM2,KM3)
!$OMP&private(WENO_X,WENO_Y,WENO_Z,UX_MAX,UY_MAX,UZ_MAX)
!$OMP&private(ADLPHIXM,ADLPHIXP,ALPHI_X,ADLPHIYM,ADLPHIYP,ALPHI_Y)
!$OMP&private(ADLPHIZM,ADLPHIZP,ALPHI_Z)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)

      IF (FUNCBODY(XPF(I),YPF(J),ZPF(K)) .LE. CRI) THEN !FIRST-ORDER UPWIND NEAR THE WALL

       IP1=IPF(I)
       IM1=IMF(I)
       JP1=JPF(J)
       JM1=JMF(J)
       KP1=KPF(K)
       KM1=KMF(K)

      !X-DIR
      UU1=UF(I,J,K)
      UU2=UF(IP1,J,K)
      UU3=0.5*(UU1+UU2)
      IF( (UU1 .GE. 0.).AND.(UU2 .GE. 0.).AND.(UU3 .GE. 0.) )THEN
        ADLPHIXM=( ALPHI_EXT(I,J,K,1)-ALPHI_EXT(IM1,J,K,1) )*SSDXF
        ALPHI_X=UU3*ADLPHIXM
      ELSE IF((UU1.LE.0.).AND.(UU2.LE.0.).AND.(UU3 .LE. 0.))THEN
        ADLPHIXP=( ALPHI_EXT(IP1,J,K,1)-ALPHI_EXT(I,J,K,1) )*SSDXF
        ALPHI_X=UU3*ADLPHIXP
      ELSE
        ADLPHIXM=( ALPHI_EXT(I,J,K,1)-ALPHI_EXT(IM1,J,K,1) )*SSDXF
        ADLPHIXP=( ALPHI_EXT(IP1,J,K,1)-ALPHI_EXT(I,J,K,1) )*SSDXF
        UX_MAX=MAX(ABS(UU1),ABS(UU2),ABS(UU3))
        ALPHI_X=0.5*(  UU3*(ADLPHIXP+ADLPHIXM)
     &                              -UX_MAX*(ADLPHIXP-ADLPHIXM)  ) 
      ENDIF

      !Y-DIR
      UU1=VF(I,J,K)
      UU2=VF(I,JP1,K)
      UU3=0.5*(UU1+UU2)
      IF( (UU1 .GE. 0.).AND.(UU2 .GE. 0.).AND.(UU3 .GE. 0.) )THEN
        ADLPHIYM=( ALPHI_EXT(I,J,K,1)-ALPHI_EXT(I,JM1,K,1) )*SSDYF
        ALPHI_Y=UU3*ADLPHIYM
      ELSE IF((UU1.LE.0.).AND.(UU2.LE.0.).AND.(UU3 .LE. 0.) )THEN
        ADLPHIYP=( ALPHI_EXT(I,JP1,K,1)-ALPHI_EXT(I,J,K,1) )*SSDYF
        ALPHI_Y=UU3*ADLPHIYP
      ELSE
        ADLPHIYM=( ALPHI_EXT(I,J,K,1)-ALPHI_EXT(I,JM1,K,1) )*SSDYF
        ADLPHIYP=( ALPHI_EXT(I,JP1,K,1)-ALPHI_EXT(I,J,K,1) )*SSDYF
        UY_MAX=MAX(ABS(UU1),ABS(UU2),ABS(UU3))
        ALPHI_Y=0.5*(UU3*(ADLPHIYP+ADLPHIYM) 
     &                                 -UY_MAX*(ADLPHIYP-ADLPHIYM) )
      ENDIF

      !Z-DIR
      UU1=WF(I,J,K)
      UU2=WF(I,J,KP1)
      UU3=0.5*(UU1+UU2)
      IF( (UU1 .GE. 0.).AND.(UU2 .GE. 0.).AND.(UU3 .GE. 0.) )THEN
        ADLPHIZM=( ALPHI_EXT(I,J,K,1)-ALPHI_EXT(I,J,KM1,1) )*SSDZF
         ALPHI_Z=UU3*ADLPHIZM
      ELSE IF((UU1.LE.0.).AND.(UU2.LE.0.).AND.(UU3 .LE. 0.) )THEN
        ADLPHIZP=( ALPHI_EXT(I,J,KP1,1)-ALPHI_EXT(I,J,K,1) )*SSDZF
        ALPHI_Z=UU3*ADLPHIZP
      ELSE
        ADLPHIZM=( ALPHI_EXT(I,J,K,1)-ALPHI_EXT(I,J,KM1,1) )*SSDZF
        ADLPHIZP=( ALPHI_EXT(I,J,KP1,1)-ALPHI_EXT(I,J,K,1) )*SSDZF  
        UZ_MAX=MAX(ABS(UU1),ABS(UU2),ABS(UU3))
        ALPHI_Z=0.5*(  UU3*(ADLPHIZP+ADLPHIZM)
     &                                -UZ_MAX*(ADLPHIZP-ADLPHIZM) )
      ENDIF

       ELSE !5TH-ORDER WENO

      !X-DIR
        IM1=IMF(I)        !!IF INTERFACE TOUCHES THE BOUNDARY, 
        IM2=IMF(IM1)        !!BOUNDARY CONDITION IS NEEDED HERE!!
        IM3=IMF(IM2)
        IP1=IPF(I)
        IP2=IPF(IP1)

      UU1=UF(I,J,K)
      UU2=UF(IP1,J,K)
      UU3=0.5*(UU1+UU2)
      IF( (UU1 .GE. 0.).AND.(UU2 .GE. 0.).AND.(UU3 .GE. 0.) )THEN
         V1=DPHI_X(IM3,J,K)
         V2=DPHI_X(IM2,J,K)
         V3=DPHI_X(IM1,J,K)
         V4=DPHI_X(I,J,K)
         V5=DPHI_X(IP1,J,K)
          CALL HJWENO(V2-V1,V3-V2,V4-V3,V5-V4,WENO_X)
         ADLPHIXM=1./12.*(-V2+7.*V3+7.*V4-V5)-WENO_X
         ALPHI_X=UU3*ADLPHIXM
      ELSE IF((UU1.LE.0.).AND.(UU2.LE.0.).AND.(UU3 .LE. 0.))THEN
         V1=DPHI_X(IP2,J,K)
         V2=DPHI_X(IP1,J,K)
         V3=DPHI_X(I,J,K)
         V4=DPHI_X(IM1,J,K)
         V5=DPHI_X(IM2,J,K)
          CALL HJWENO(V1-V2,V2-V3,V3-V4,V4-V5,WENO_X)
         ADLPHIXP=1./12.*(-V5+7.*V4+7.*V3-V2)+WENO_X
         ALPHI_X=UU3*ADLPHIXP
      ELSE
         V1=DPHI_X(IM3,J,K)
         V2=DPHI_X(IM2,J,K)
         V3=DPHI_X(IM1,J,K)
         V4=DPHI_X(I,J,K)
         V5=DPHI_X(IP1,J,K)
          CALL HJWENO(V2-V1,V3-V2,V4-V3,V5-V4,WENO_X)
         ADLPHIXM=1./12.*(-V2+7.*V3+7.*V4-V5)-WENO_X
         V1=DPHI_X(IP2,J,K)
         V2=DPHI_X(IP1,J,K)
         V3=DPHI_X(I,J,K)
         V4=DPHI_X(IM1,J,K)
         V5=DPHI_X(IM2,J,K)
          CALL HJWENO(V1-V2,V2-V3,V3-V4,V4-V5,WENO_X)
         ADLPHIXP=1./12.*(-V5+7.*V4+7.*V3-V2)+WENO_X
         UX_MAX=MAX(ABS(UU1),ABS(UU2),ABS(UU3))
         ALPHI_X=0.5*(  UU3*(ADLPHIXP+ADLPHIXM)
     &                              -UX_MAX*(ADLPHIXP-ADLPHIXM)  ) 
      ENDIF

      !Y-DIR
        JM1=JMF(J)
        JM2=JMF(JM1)
        JM3=JMF(JM2)
        JP1=JPF(J)
        JP2=JPF(JP1)
      UU1=VF(I,J,K)
      UU2=VF(I,JP1,K)
      UU3=0.5*(UU1+UU2)
      IF( (UU1 .GE. 0.).AND.(UU2 .GE. 0.).AND.(UU3 .GE. 0.) )THEN
        V1=DPHI_Y(I,JM3,K)
        V2=DPHI_Y(I,JM2,K)
        V3=DPHI_Y(I,JM1,K)
        V4=DPHI_Y(I,J,K)
        V5=DPHI_Y(I,JP1,K)
         CALL HJWENO(V2-V1,V3-V2,V4-V3,V5-V4,WENO_Y)
        ADLPHIYM=1./12.*(-V2+7.*V3+7.*V4-V5)-WENO_Y 
        ALPHI_Y=UU3*ADLPHIYM
      ELSE IF((UU1.LE.0.).AND.(UU2.LE.0.).AND.(UU3 .LE. 0.) )THEN
        V1=DPHI_Y(I,JP2,K)
        V2=DPHI_Y(I,JP1,K)
        V3=DPHI_Y(I,J,K)  
        V4=DPHI_Y(I,JM1,K)
        V5=DPHI_Y(I,JM2,K)   
         CALL HJWENO(V1-V2,V2-V3,V3-V4,V4-V5,WENO_Y)
        ADLPHIYP=1./12.*(-V5+7.*V4+7.*V3-V2)+WENO_Y
        ALPHI_Y=UU3*ADLPHIYP    
      ELSE
        V1=DPHI_Y(I,JM3,K)
        V2=DPHI_Y(I,JM2,K)
        V3=DPHI_Y(I,JM1,K)
        V4=DPHI_Y(I,J,K)
        V5=DPHI_Y(I,JP1,K)
         CALL HJWENO(V2-V1,V3-V2,V4-V3,V5-V4,WENO_Y)
        ADLPHIYM=1./12.*(-V2+7.*V3+7.*V4-V5)-WENO_Y 
        V1=DPHI_Y(I,JP2,K)
        V2=DPHI_Y(I,JP1,K)
        V3=DPHI_Y(I,J,K)  
        V4=DPHI_Y(I,JM1,K)
        V5=DPHI_Y(I,JM2,K)  
         CALL HJWENO(V1-V2,V2-V3,V3-V4,V4-V5,WENO_Y)
        ADLPHIYP=1./12.*(-V5+7.*V4+7.*V3-V2)+WENO_Y
        UY_MAX=MAX(ABS(UU1),ABS(UU2),ABS(UU3))
        ALPHI_Y=0.5*(UU3*(ADLPHIYP+ADLPHIYM) 
     &                                 -UY_MAX*(ADLPHIYP-ADLPHIYM) )
      ENDIF

      !Z-DIR
        KM1=KMF(K)
        KM2=KMF(KM1)
        KM3=KMF(KM2)
        KP1=KPF(K)
        KP2=KPF(KP1)
      UU1=WF(I,J,K)
      UU2=WF(I,J,KP1)
      UU3=0.5*(UU1+UU2)
      IF( (UU1 .GE. 0.).AND.(UU2 .GE. 0.).AND.(UU3 .GE. 0.) )THEN
        V1=DPHI_Z(I,J,KM3)
        V2=DPHI_Z(I,J,KM2)
        V3=DPHI_Z(I,J,KM1)
        V4=DPHI_Z(I,J,K)
        V5=DPHI_Z(I,J,KP1)
         CALL HJWENO(V2-V1,V3-V2,V4-V3,V5-V4,WENO_Z)
         ADLPHIZM=1./12.*(-V2+7.*V3+7.*V4-V5)-WENO_Z     
         ALPHI_Z=UU3*ADLPHIZM
      ELSE IF((UU1.LE.0.).AND.(UU2.LE.0.).AND.(UU3 .LE. 0.) )THEN
        V1=DPHI_Z(I,J,KP2)
        V2=DPHI_Z(I,J,KP1)
        V3=DPHI_Z(I,J,K)
        V4=DPHI_Z(I,J,KM1)
        V5=DPHI_Z(I,J,KM2)
         CALL HJWENO(V1-V2,V2-V3,V3-V4,V4-V5,WENO_Z)
        ADLPHIZP=1./12.*(-V5+7.*V4+7.*V3-V2)+WENO_Z  
        ALPHI_Z=UU3*ADLPHIZP
      ELSE
        V1=DPHI_Z(I,J,KM3)
        V2=DPHI_Z(I,J,KM2)
        V3=DPHI_Z(I,J,KM1)
        V4=DPHI_Z(I,J,K)
        V5=DPHI_Z(I,J,KP1)
         CALL HJWENO(V2-V1,V3-V2,V4-V3,V5-V4,WENO_Z)
        ADLPHIZM=1./12.*(-V2+7.*V3+7.*V4-V5)-WENO_Z     
        V1=DPHI_Z(I,J,KP2)
        V2=DPHI_Z(I,J,KP1)
        V3=DPHI_Z(I,J,K)
        V4=DPHI_Z(I,J,KM1)
        V5=DPHI_Z(I,J,KM2)
         CALL HJWENO(V1-V2,V2-V3,V3-V4,V4-V5,WENO_Z)
        ADLPHIZP=1./12.*(-V5+7.*V4+7.*V3-V2)+WENO_Z      
        UZ_MAX=MAX(ABS(UU1),ABS(UU2),ABS(UU3))
        ALPHI_Z=0.5*(  UU3*(ADLPHIZP+ADLPHIZM)
     &                                -UZ_MAX*(ADLPHIZP-ADLPHIZM) )
      ENDIF
      
      ENDIF !IF (FUNCBODY(XPF(I),YPF(J),ZPF(K)) .LE. CRI) THEN

!C***** TIME ADVANCEMENT USING EE
      ALPHI_NEW(NN,N)=ALPHI_EXT(I,J,K,1)
     &                                  -DTLVS*(ALPHI_X+ALPHI_Y+ALPHI_Z)

       ENDDO
       ENDDO
       
       
      IF (IPHS .EQ. 1) THEN
      DO N=1,NN1
      DO NN=1,NUMA(N,LLVS)            
      if (alphi(nn,n,llvs) .gt. 0d0) then 
      ALPHI_NEW(NN,N)=ALPHI_NEW(NN,N)-DTLVS*MFLVS(NN,N,LLVS)/DENR
      else   
      ALPHI_NEW(NN,N)=ALPHI_NEW(NN,N)-DTLVS*MFLVS(NN,N,LLVS)
      endif
      ENDDO
      ENDDO
      ENDIF       

       
!C=======================SOLVING_TRANSPORT_EQN_END======================C

       DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
          I=I_B(NN,N,LLVS)
          J=J_B(NN,N,LLVS)
          K=K_B(NN,N,LLVS)
           ALPHI_EXT(I,J,K,1)=ALPHI_NEW(NN,N)   
       ENDDO
       ENDDO

       RETURN
       END       
       
       

!C*******************************************************************
      SUBROUTINE WALL_INTERACTION(NN1,LLVS)
!C******************************************************************* 
      !IB WALL AND LEVEL-SET INTERACTION
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING

!      IF (N3FM .EQ. 1) THEN
!       HEIGHT_CRI=SQRT(SDXF**2+SDYF**2)
!      ELSE
!       HEIGHT_CRI=SQRT(SDXF**2+SDYF**2+SDZF**2)
!      ENDIF
!
!       HEIGHT_SPARE=0.5*HEIGHT_CRI
!
!       DO N=1,NN1
!!$OMP PARALLEL DO private(I,J,K,BODY)
!        DO NN=1,NUMA(N,LLVS)
!         I=I_B(NN,N,LLVS)
!         J=J_B(NN,N,LLVS)
!         K=K_B(NN,N,LLVS)  
!         BODY=FUNCBODY(XPF(I),YPF(J),ZPF(K))
!
!         IF ( ABS(BODY) .LE. HEIGHT_CRI ) THEN         !LEVEL-SET, IB WALL BOUNDARY CONDITION
!
!          IF( BODY .LT. 0.) THEN
!           IF (ALPHI_EXT(I,J,K,LLVS) .LE. -BODY ) THEN
!            ALPHI_EXT(I,J,K,LLVS)=-BODY+HEIGHT_SPARE
!           ENDIF
!          ELSE IF( BODY .GT. 0.) THEN
!           IF ( -ALPHI_EXT(I,J,K,LLVS) .GE. BODY ) THEN
!            ALPHI_EXT(I,J,K,LLVS)=-BODY-HEIGHT_SPARE
!           ENDIF
!          ENDIF
!
!         ENDIF !IF ( ABS(BODY) .LE. HEIGHT_CRI ) THEN
!
!        ENDDO
!       ENDDO






c-----contact line (not finished)
!      IF (N3FM .EQ. 1) THEN
!       HEIGHT_CRI=SQRT(SDXF**2+SDYF**2)
!      ELSE
!       HEIGHT_CRI=SQRT(SDXF**2+SDYF**2+SDZF**2)
!      ENDIF
!
!       DO N=1,NN1
!!$OMP PARALLEL DO private(I,J,K,BODY1,BODY2,BODY3,BODY4)
!!$OMP&private(ANORI,ANOR1,ANOR2,L,II,JJ,I2,J2,I3,J3,X2,Y2,X3,Y3)
!!$OMP&private(DPX,DPY,DPDN,DENO,B1,B2,B3)
!        DO NN=1,NUMA(N,LLVS)
!         I=I_B(NN,N,LLVS)
!         J=J_B(NN,N,LLVS)
!         K=K_B(NN,N,LLVS)  
!
!!         IF ( ABS(FUNCBODY(XPF(I),YPF(J),ZPF(K))) .LE. HEIGHT_CRI ) THEN
!         IF ( FUNCBODY(XPF(I),YPF(J),ZPF(K)) .LE. HEIGHT_CRI ) THEN
!           
!         BODY1=FUNCBODY(XPF(IPF(I)),YPF(J),ZPF(K))
!         BODY2=FUNCBODY(XPF(IMF(I)),YPF(J),ZPF(K))
!         BODY3=FUNCBODY(XPF(I),YPF(JPF(J)),ZPF(K))
!         BODY4=FUNCBODY(XPF(I),YPF(JMF(J)),ZPF(K))
!
!         ANOR1=0.5*(BODY1-BODY2)*SSDXF
!         ANOR2=0.5*(BODY3-BODY4)*SSDYF
!         ANORI=1./SQRT(ANOR1**2+ANOR2**2)
!         ANOR1=ANOR1*ANORI
!         ANOR2=ANOR2*ANORI
!
!         II=INT(SIGN(1.,ANOR1))
!         JJ=INT(SIGN(1.,ANOR2))
!         
!         I2=I+II
!         J2=J+JJ
!         X2=XPF(I2)
!         Y2=YPF(J2)
!
!         I3=I2+II
!         J3=J2+JJ
!         X3=XPF(I3)
!         Y3=YPF(J3)
!
!         DPX=ALPHI_EXT(I3+1,J3,K,LLVS)-ALPHI_EXT(I3-1,J3,K,LLVS)
!         DPY=ALPHI_EXT(I3,J3+1,K,LLVS)-ALPHI_EXT(I3,J3-1,K,LLVS)
!         DPDN=DPX*ANOR1+DPY*ANOR2
!
!        DENO=1./(ANOR1*(Y3-Y2)-ANOR2*(X3-X2))
!
!        B1=((-X2*Y3+X3*Y2)*DPDN
!     &              +(ANOR1*Y3-ANOR2*X3)*ALPHI_EXT(I2,J2,K,LLVS)
!     &              -(-ANOR1*Y2+ANOR2*X2)*ALPHI_EXT(I3,J3,K,LLVS))*DENO
!        B2=((Y3-Y2)*DPDN+ANOR2*ALPHI_EXT(I2,J2,K,LLVS)
!     &                  -ANOR2*ALPHI_EXT(I3,J3,K,LLVS))*DENO
!        B3=(-(X3-X2)*DPDN-ANOR1*ALPHI_EXT(I2,J2,K,LLVS)
!     &                   +ANOR1*ALPHI_EXT(I3,J3,K,LLVS))*DENO
!     
!         ALPHI_EXT(I,J,K,LLVS)=B1+B2*XPF(I)+B3*YPF(J)
!
!
!         WRITE(*,*) ALPHI_EXT(I2,J2,K,LLVS),B1+B2*XPF(I2)+B3*YPF(J2)
!         
!         ENDIF  !IF ( ABS(BODY) .LE. HEIGHT_CRI ) THEN
!
!        ENDDO
!       ENDDO


       RETURN
       END 
    


!C*******************************************************************
      SUBROUTINE BAND_GENERATION(LLVS)
!C******************************************************************* 
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      ! use flow_var
      
      IMPLICIT NONE

      INTEGER MASK_S(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      
      INTEGER LLVS,NN2
      
      INTEGER IX_TOUCH1,IX_TOUCH2,IY_TOUCH1,IY_TOUCH2
     &,IZ_TOUCH1,IZ_TOUCH2
      INTEGER IX_TOUCH_N,IY_TOUCH_N,IZ_TOUCH_N
      INTEGER I,J,K,II,JJ,KK,I2,J2,K2,N,NN,NM,NUMS
      REAL DX,DY,DZ
     
!C      MASK_S = 0    !NOT ACTIVATE CELL
!C      MASK_S = 1    !ACTIVATED CELL & NEAR INTERFACE
      nn2=n_max-1
      
       IX_TOUCH1=0
       IX_TOUCH2=0
       IY_TOUCH1=0
       IY_TOUCH2=0
       IZ_TOUCH1=0
       IZ_TOUCH2=0
       IX_TOUCH_N=N_MAX+1
       IY_TOUCH_N=N_MAX+1
       IZ_TOUCH_N=N_MAX+1

!      NUMA_OLD=0
!      DO N=1,N_MAX
!        NUMA_OLD=NUMA_OLD+NUMA(N,LLVS)
!      ENDDO

!$OMP PARALLEL DO private(I,J)
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
      MASK_S(I,J,K)=0
      ENDDO
      ENDDO
      ENDDO

!C===============================MAKE S-TUBE============================C
      NUMS=0
       DO N=1,NN2
!c$OMP PARALLEL DO private(I,J,K,NUMS)  !ORDER
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
!AT THE TUBE BOUNDARY ALPHI1*ALPHI2 CAN BE ZERO WHICH IS NOT DESIRABLE.
!ASSUME ALPHI NEVER BE ZERO WITHIN COMPUTATION BECAUSE OF NUMERICAL ERROR.
!SO I USE 'LESS THAN' NOT 'LESS EQUAL' TO ZERO
       IF(    (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(IPF(I),J,K,1) .LT. 0.)    
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(IMF(I),J,K,1) .LT. 0.)
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,JPF(J),K,1) .LT. 0.)
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,JMF(J),K,1) .LT. 0.)     
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,J,KPF(K),1) .LT. 0.)     
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,J,KMF(K),1) .LT. 0.)
     &                                                            ) THEN
           NUMS=NUMS+1
           I_B(NUMS,1,LLVS)=I
           J_B(NUMS,1,LLVS)=J
           K_B(NUMS,1,LLVS)=K
           MASK_S(I,J,K)=1

       ELSE
           MASK_S(I,J,K)=2
       ENDIF

       ENDDO
      ENDDO
           NUMA(1,LLVS)=NUMS
           
           
           
!C===============================MAKE S-TUBE============================C

        DO 1003 N=2,N_MAX
          NUMS=0
!C==========================makes a I_B,J_B,K_B=========================C
!c!$OMP PARALLEL DO private(I,J,K)  !MASK_S IS SHOULD BE SYNCRONIZED EVERY STEP.. HIGH COST. ORDER IS IMPORTANT..
          NM=N-1
         DO 1004 NN=1,NUMA(NM,LLVS)  !AT S-TUBE
          I=I_B(NN,NM,LLVS)
          J=J_B(NN,NM,LLVS)
          K=K_B(NN,NM,LLVS)  

         DO K2=K-1,K+1
            KK=K2
            IF (IPZ .EQ. 1) THEN
             IF (KK .EQ. 0) KK=N3FM
             IF (KK .EQ. N3F) KK=1  
            ELSE
             IF(K2 .EQ. 1) THEN        !IF INTERFACE TOUCHES THE BOUNDARY
              IZ_TOUCH1=IZ_TOUCH1+1
              IZ_TOUCH_N=MIN(IZ_TOUCH_N,N)
              goto 1005
             ELSE IF(K2 .EQ. N3F) THEN
              IZ_TOUCH2=IZ_TOUCH2+1
              IZ_TOUCH_N=MIN(IZ_TOUCH_N,N)
              goto 1005
             ENDIF
            ENDIF
         DO J2=J-1,J+1
            JJ=J2
            IF (IPY .EQ. 1) THEN
             IF (JJ .EQ. 0) JJ=N2FM
             IF (JJ .EQ. N2F) JJ=1  
            ELSE
             IF(J2 .EQ. 1) THEN 
              IY_TOUCH1=IY_TOUCH1+1
              IY_TOUCH_N=MIN(IY_TOUCH_N,N)
              goto 1006
             ELSE IF(J2 .EQ. N2F) THEN
              IY_TOUCH2=IY_TOUCH2+1
              IY_TOUCH_N=MIN(IY_TOUCH_N,N)
              goto 1006
             ENDIF
            ENDIF
          DO I2=I-1,I+1
            II=I2
            IF (IPX .EQ. 1) THEN
             IF (II .EQ. 0) II=N1FM
             IF (II .EQ. N1F) II=1  
            ELSE
             IF(I2 .EQ. 1) THEN 
               IX_TOUCH1=IX_TOUCH1+1
               IX_TOUCH_N=MIN(IX_TOUCH_N,N)
              goto 1007
             ELSE IF (I2 .EQ. N1FM) THEN
               IX_TOUCH2=IX_TOUCH2+1
               IX_TOUCH_N=MIN(IX_TOUCH_N,N)
              goto 1007
             ENDIF
            ENDIF

!       IF ( MASK_S(II,JJ,KK) .EQ. 1 ) ->JUST GO THROUGH

       IF ( MASK_S(II,JJ,KK) .NE. 1 ) THEN

        !IF ( MASK_S(II,JJ,KK) .EQ. 2 )  -> ALREADY IN I_B,J_B,K_B & HAVE ALPHI VALUE 

!C-----NEW_CELL ASSIGN
      IF( (MASK_S(II,JJ,KK).EQ. 0) .OR. (N .GE. NN2) ) THEN

      !ALPHI ASSIGN
       !this if statesment can have problem, when z-dir has no difference.
       !because only k2=0 is activated, but k2=n3f is no acticated.
       !that is direction is diffrent only at the k2=0
       !if there is difference along z-dir, it may not have problem in my thingking.
       ! z방향으로 평평하면 참조하는 cell에서의 거리랑 수직방향으로의 거리가 일치하지 않아서 문제가 된다. 
       ! --> re-initialization과 꼭 함께 사용해야한다. 
      IF ( K2 .EQ. 0) THEN
        IF (J2 .EQ. 0) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ENDIF
        ELSE IF (J2 .EQ. N2F) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ENDIF
        ELSE
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ENDIF
        ENDIF
       ELSE IF (K2 .EQ. N3F) THEN
         IF (J2 .EQ. 0) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ENDIF
        ELSE IF (J2 .EQ. N2F) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ENDIF
        ELSE
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ENDIF
        ENDIF
        ELSE
        IF (J2 .EQ. 0) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ENDIF
        ELSE IF (J2 .EQ. N2F) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ENDIF
        ELSE
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-ZPF(K)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-ZPF(K)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-ZPF(K)
          ENDIF
        ENDIF
      ENDIF

       IF (N3FM .EQ. 1) DZ=0. !2D CASE
       
       ! z방향 평행할 떼 (z방향만 이렇게 처리되어 있으니까, x-y 방향도 이렇게 해야하면 따로 건드려야함)
       if (alphi_ext(imf(i),j,k,1) .eq. alphi_ext(ipf(i),j,k,1) .and.
     &     alphi_exT(i,jmf(j),k,1) .eq. alphi_ext(i,jpf(j),k,1) ) then 
       dx=0d0
       dy=0d0
       dz=sdzf       
       endif

       IF (ALPHI_EXT(I,J,K,1) .GE. 0.)  THEN
       ALPHI_EXT(II,JJ,KK,1)=ALPHI_EXT(I,J,K,1)
     &                              +DSQRT(DX**2d0+DY**2d0+DZ**2d0) 
       ELSE
       ALPHI_EXT(II,JJ,KK,1)=ALPHI_EXT(I,J,K,1)
     &                              -DSQRT(DX**2d0+DY**2d0+DZ**2d0) 
       ENDIF

       ENDIF !IF( (MASK_S(II,JJ,KK).EQ. 0) .OR. (N .GE. 11) ) THEN

!C-----NEW_CELL ASSIGN
          NUMS=NUMS+1
          I_B(NUMS,N,LLVS)=II
          J_B(NUMS,N,LLVS)=JJ
          K_B(NUMS,N,LLVS)=KK
          MASK_S(II,JJ,KK)=1  

      ENDIF !IF ( MASK_S(II,JJ,KK) .NE. 1 )

 1007   CONTINUE
       ENDDO
 1006   CONTINUE
       ENDDO
 1005   CONTINUE
       ENDDO

 1004  CONTINUE
           NUMA(N,LLVS)=NUMS
!C==========================makes a I_B,J_B,K_B=========================C
 1003  CONTINUE

!C===========================BOUNDARY CONDITION=========================C
      !IF LEVEL-SET FUNCT. TOUCHES THE DOMAIN BOUNDARY, THIS ALGORITHM IS ACTIVATED.
      !THIS IS NOT IN TUBE, BUT USED WHEN CALCULATE TRANSPORT & RE-INITIALIZATION.
      IF (IX_TOUCH_N .NE. N_MAX+1 ) THEN
          WRITE(*,790) IX_TOUCH_N,IX_TOUCH1+IX_TOUCH2
 790  FORMAT ('X_BOUNDARY TOUCH. N_MIN=',I5,',  TOTAL TOUCH CELL=',I10)
       ENDIF
      IF (IY_TOUCH_N .NE. N_MAX+1 ) THEN
          WRITE(*,791) IY_TOUCH_N,IY_TOUCH1+IY_TOUCH2
 791  FORMAT ('Y_BOUNDARY TOUCH. N_MIN=',I5,',  TOTAL TOUCH CELL=',I10)
       ENDIF
      IF (IZ_TOUCH_N .NE. N_MAX+1 ) THEN
          WRITE(*,792) IZ_TOUCH_N,IZ_TOUCH1+IZ_TOUCH2
 792  FORMAT ('Z_BOUNDARY TOUCH. N_MIN=',I5,',  TOTAL TOUCH CELL=',I10)
       ENDIF

!C===========================BOUNDARY CONDITION=========================C
!      NUMA_NEW=0
!      DO N=1,N_MAX
!      NUMA_NEW=NUMA_NEW+NUMA(N,LLVS)
!      ENDDO
!
!       WRITE(*,50) NUMA_OLD,NUMA_NEW,NUMA_NEW-NUMA_OLD
! 50   FORMAT('PREVIOUS CELL: ',I8,' PRESENT CELL:',I8,' VARIATION: ',I8)

       NUMA_MAX=0
       DO N=1,N_MAX
         NUMA_MAX=max(NUMA_MAX,NUMA(N,LLVS))
       ENDDO
!       write(*,*) 'numa_max=',numa_max,mf_band       
!        WRITE(*,*) 'NUMA_MAX/MF_BAND=',FLOAT(NUMA_MAX)/FLOAT(MF_BAND)
       IF (NUMA_MAX .GT. MF_BAND) THEN
        WRITE(*,*) 'NUMA_MAX/MF_BAND=',FLOAT(NUMA_MAX)/FLOAT(MF_BAND)
        WRITE(*,*) 'MF_BAND IS NOT ENOUGH!!'
        WRITE(*,*) 'SIMULATION IS STOPED!!'
        STOP
       ENDIF
       
       ! if (llvs .eq. 1) the      ! if (llvs .eq. 1) then
! C=====SAVE
       ! OPEN(143,FILE='0alphi_AFTER_BANDG.DAT')
        ! WRITE(143,*) 'VARIABLE="N","NN","I_B","J_B","K_B","ALPHI"'
        ! DO N=1,14
        ! DO NN=1,NUMA(N,llvs)
          ! i=I_B(NN,N,llvs)
          ! j=J_B(NN,N,llvs)
          ! k=K_B(NN,N,llvs)
        ! WRITE(143,145) n,NN,i,j,k,alphi_ext(i,j,k,1)
        ! ENDDO
        ! ENDDO
         ! CLOSE(143)
 ! 145  FORMAT(5i8,1F15.8)

       
!C=====SAVE
!      write(*,*) '*************************************m=',m
!       OPEN(146,FILE='0alphi_contour_AFTER_BANDG.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","alphi"'
!      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!      DO K=1,N3FM
!      DO J=1,N2FM
!      DO I=1,N1FM
!        WRITE(146,147) XPF(I),YPF(J),ZPF(K),alphi_ext(i,j,k,1)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
! 147  FORMAT(4F15.8)
! 
       ! endif
       
      RETURN
      END

!C*******************************************************************
      SUBROUTINE REINITIALIZATION(NN1,NN2,ILIMIT,LLVS)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE

       REAL ALPHIO(-2:M1F+3,-2:M2F+3,-2:M3F+3)
       REAL ALPHI_RK3(-2:M1F+3,-2:M2F+3,-2:M3F+3)
       
       INTEGER NN1,NN2,ILIMIT,LLVS
       
       INTEGER I,J,K,N,NN,ITER_RE
       REAL SM_GRID,REDTL

       DO N=1,NN2
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         ALPHIO(I,J,K)=ALPHI_EXT(I,J,K,1)
        ENDDO
       ENDDO

      IF (N3FM .EQ. 1 ) THEN
        SM_GRID=MIN(SDXF,SDYF)
      ELSE
        SM_GRID=MIN(SDXF,SDYF,SDZF)
      ENDIF

        REDTL=0.5*SM_GRID   !CFL=<0.5 because of stability of Godunov scheme

!         ILIMIT=22                 !MAX_ITER=11/CFL (2008 herrmann)
        !ILIMIT=INT(1.5/SM_GRID)   !MAX_ITER=CFL*(3/DX)  (2006 herrmann)
!C==================MAIN REINITIALIZATION CAL._START=====================C
      DO 1000 ITER_RE=1,ILIMIT
!C-----3RD_ORDER TVD RK

       DO N=1,NN2
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         ALPHI_RK3(I,J,K)=ALPHI_EXT(I,J,K,1)
        ENDDO
       ENDDO

        CALL REINITIALIZATION_RK3(REDTL,ALPHIO,NN1,NN2,LLVS)
        CALL REINITIALIZATION_RK3(REDTL,ALPHIO,NN1,NN2,LLVS)

       DO N=1,NN2
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         ALPHI_EXT(I,J,K,1)=0.75*ALPHI_RK3(I,J,K)
     &                                     +0.25*ALPHI_EXT(I,J,K,1)
        ENDDO
      ENDDO

       CALL REINITIALIZATION_RK3(REDTL,ALPHIO,NN1,NN2,LLVS)

       DO N=1,NN2
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)                      
         ALPHI_EXT(I,J,K,1)=1./3.*ALPHI_RK3(I,J,K)
     &                                     +2./3.*ALPHI_EXT(I,J,K,1)          
        ENDDO
      ENDDO

!!       !STEADY-STATE EVALUATION
!      IF (MOD(ITER_RE,10) .EQ. 0 .OR. ITER_RE .EQ. ILIMIT  ) THEN
!
!        DALPHIMAX=0.
!        DALPHIMIN=1.
!
!       DO N=1,NN1                ! T-BAND ONLY (HERRMANN 2008)
!        DALPHIMAX_OMP=0.
!        DALPHIMIN_OMP=1.
!!$OMP PARALLEL DO private(I,J,K,DALPHIX,DALPHIY,DALPHIZ,DALPHI_SUM)
!!$OMP&reduction(MAX:DALPHIMAX_OMP)
!!$OMP&reduction(MIN:DALPHIMIN_OMP)
!        DO NN=1,NUMA(N,LLVS)
!         I=I_B(NN,N,LLVS)
!         J=J_B(NN,N,LLVS)
!         K=K_B(NN,N,LLVS)
!        DALPHIX=0.5*( ALPHI_EXT(IPF(I),J,K,LLVS)
!    &                               -ALPHI_EXT(IMF(I),J,K,LLVS) )*SSDXF
!        DALPHIY=0.5*( ALPHI_EXT(I,JPF(J),K,LLVS)
!    &                               -ALPHI_EXT(I,JMF(J),K,LLVS) )*SSDYF
!        DALPHIZ=0.5*( ALPHI_EXT(I,J,KPF(K),LLVS)
!    &                               -ALPHI_EXT(I,J,KMF(K),LLVS) )*SSDZF
!        DALPHI_SUM=DALPHIX**2+DALPHIY**2+DALPHIZ**2 
!
!        DALPHIMAX_OMP=MAX(DALPHI_SUM,DALPHIMAX_OMP)
!        DALPHIMIN_OMP=MIN(DALPHI_SUM,DALPHIMIN_OMP)
!       ENDDO
!         DALPHIMAX=MAX(DALPHIMAX,DALPHIMAX_OMP)
!         DALPHIMIN=MIN(DALPHIMIN,DALPHIMIN_OMP)
!       ENDDO
!         DALPHIMAX=SQRT(DALPHIMAX)
!         DALPHIMIN=SQRT(DALPHIMIN)
!
!      IF( DALPHIMAX-1. .LE. 0.1 .AND. 1.-DALPHIMIN .LE. 0.1) THEN
!        GOTO 1000
!      ENDIF
!      ENDIF
!
 1000 CONTINUE
! 
!         WRITE(*,900) DALPHIMAX,DALPHIMIN
! 900   FORMAT('DPHIMAX=',F10.5,'  DPHIMIN=',F10.5)
!C==================MAIN REINITIALIZATION CAL._END======================C  
!C=====SAVE
!       OPEN(146,FILE='0alphi_contour_AFTER_REINI.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","alphi"'
!      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!      DO K=1,N3FM
!      DO J=1,N2FM
!      DO I=1,N1FM
!        WRITE(146,147) XPF(I),YPF(J),ZPF(K),alphi(i,j,k)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
! 147  FORMAT(4F15.8) 

        RETURN
       END


!C*******************************************************************
      SUBROUTINE REINITIALIZATION_RK3(REDTL,ALPHIO,NN1,NN2,LLVS)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE

       REAL ALPHI_NEW(MF_BAND,M_MAX)
       REAL ALPHIO(-2:M1F+3,-2:M2F+3,-2:M3F+3)
       REAL DPHI_X(-2:M1F+3,-2:M2F+3,-2:M3F+3)
     &      ,DPHI_Y(-2:M1F+3,-2:M2F+3,-2:M3F+3)
     &      ,DPHI_Z(-2:M1F+3,-2:M2F+3,-2:M3F+3)

      INTEGER NN1,NN2,LLVS
      REAL REDTL
      
      INTEGER N,NN,I,J,K,IP1,IP2,IM1,IM2,IM3
      INTEGER JP1,JP2,JM1,JM2,JM3,KP1,KP2,KM1,KM2,KM3
      REAL EPS,DALPHIX,DALPHIY,DALPHIZ,DALPHI,ASIGN
      REAL V1,V2,V3,V4,V5,WENO_X,WENO_Y,WENO_Z,DRHS
      REAL ADLPHIXM,ADLPHIXP,ALPHI_X,ADLPHIYM,ADLPHIYP,ALPHI_Y
      REAL ADLPHIZM,ADLPHIZP,ALPHI_Z
      
        IF (N3F .EQ. 2 ) THEN
         EPS=MIN(SDXF,SDYF)**2
        ELSE
         EPS=MIN(SDXF,SDYF,SDZF)**2
        ENDIF

      IF (NN1 .NE. 0) THEN
!C======================TO AVOID REPEATED CAL._START====================C
      !DPHI_X,DPHI_Y,DPHI_Z IS DEFINED AT CELL-INTERFACE
       DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         DPHI_X(I,J,K)=(ALPHI_EXT(IPF(I),J,K,1)
     &                                    -ALPHI_EXT(I,J,K,1))*SSDXF
         DPHI_Y(I,J,K)=(ALPHI_EXT(I,JPF(J),K,1)
     &                                    -ALPHI_EXT(I,J,K,1))*SSDYF        
         DPHI_Z(I,J,K)=(ALPHI_EXT(I,J,KPF(K),1)
     &                                    -ALPHI_EXT(I,J,K,1))*SSDZF
         ENDDO
       ENDDO

      N=NN1               ! FOR THE LAST LAYER OF T-BAND
!$OMP PARALLEL DO private(I,J,K,IP1,IP2,IM1,IM2,IM3)
!$OMP&private(JP1,JP2,JM1,JM2,JM3,KP1,KP2,KM1,KM2,KM3)
       DO NN=1,NUMA(N,LLVS)  ! FOR ALL LLVS
        I=I_B(NN,N,LLVS)
        J=J_B(NN,N,LLVS)
        K=K_B(NN,N,LLVS)
        IM1=IMF(I)        !!IF INTERFACE TOUCHES THE BOUNDARY, 
        IM2=IMF(IM1)      !!BOUNDARY CONDITION IS NEEDED HERE!!
        IM3=IMF(IM2)
        IP1=IPF(I)
        IP2=IPF(IP1)
      DPHI_X(IM3,J,K)=(ALPHI_EXT(IM2,J,K,1)
     &                                   -ALPHI_EXT(IM3,J,K,1))*SSDXF            
      DPHI_X(IM2,J,K)=(ALPHI_EXT(IM1,J,K,1)
     &                                   -ALPHI_EXT(IM2,J,K,1))*SSDXF
      DPHI_X(IM1,J,K)=(ALPHI_EXT(I,J,K,1)
     &                                   -ALPHI_EXT(IM1,J,K,1))*SSDXF
      DPHI_X(IP1,J,K)=(ALPHI_EXT(IP2,J,K,1)
     &                                   -ALPHI_EXT(IP1,J,K,1))*SSDXF
      DPHI_X(IP2,J,K)=(ALPHI_EXT(IPF(IP2),J,K,1)
     &                                   -ALPHI_EXT(IP2,J,K,1))*SSDXF
        JM1=JMF(J)
        JM2=JMF(JM1)
        JM3=JMF(JM2)
        JP1=JPF(J)
        JP2=JPF(JP1)
      DPHI_Y(I,JM3,K)=(ALPHI_EXT(I,JM2,K,1)
     &                                   -ALPHI_EXT(I,JM3,K,1))*SSDYF
      DPHI_Y(I,JM2,K)=(ALPHI_EXT(I,JM1,K,1)
     &                                   -ALPHI_EXT(I,JM2,K,1))*SSDYF
      DPHI_Y(I,JM1,K)=(ALPHI_EXT(I,J,K,1)
     &                                   -ALPHI_EXT(I,JM1,K,1))*SSDYF
      DPHI_Y(I,JP1,K)=(ALPHI_EXT(I,JP2,K,1)
     &                                   -ALPHI_EXT(I,JP1,K,1))*SSDYF
      DPHI_Y(I,JP2,K)=(ALPHI_EXT(I,JPF(JP2),K,1)
     &                                   -ALPHI_EXT(I,JP2,K,1))*SSDYF
        KM1=KMF(K)
        KM2=KMF(KM1)
        KM3=KMF(KM2)
        KP1=KPF(K)
        KP2=KPF(KP1)
      DPHI_Z(I,J,KM3)=(ALPHI_EXT(I,J,KM2,1)
     &                                   -ALPHI_EXT(I,J,KM3,1))*SSDZF
      DPHI_Z(I,J,KM2)=(ALPHI_EXT(I,J,KM1,1)
     &                                   -ALPHI_EXT(I,J,KM2,1))*SSDZF
      DPHI_Z(I,J,KM1)=(ALPHI_EXT(I,J,K,1)
     &                                   -ALPHI_EXT(I,J,KM1,1))*SSDZF
      DPHI_Z(I,J,KP1)=(ALPHI_EXT(I,J,KP2,1)
     &                                   -ALPHI_EXT(I,J,KP1,1))*SSDZF
      DPHI_Z(I,J,KP2)=(ALPHI_EXT(I,J,KPF(KP2),1)
     &                                   -ALPHI_EXT(I,J,KP2,1))*SSDZF
       ENDDO
!C======================TO AVOID REPEATED CAL._END======================C
      ENDIF !IF ( NN1 .NE. 0 ) 

!C======================= T-TUBE :  5TH-ORDER WEO=======================C
       DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K,V1,V2,V3,V4,V5,IP1,IP2,IM1,IM2,IM3)
!$OMP&private(JP1,JP2,JM1,JM2,JM3,KP1,KP2,KM1,KM2,KM3)
!$OMP&private(DALPHIX,DALPHIY,DALPHIZ,DALPHI,ASIGN,DRHS)
!$OMP&private(WENO_X,WENO_Y,WENO_Z,ADLPHIXM,ADLPHIXP,ALPHI_X)
!$OMP&private(ADLPHIYM,ADLPHIYP,ALPHI_Y,ADLPHIZM,ADLPHIZP,ALPHI_Z)
        DO NN=1,NUMA(N,LLVS)

         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)

         IM1=IMF(I)
         IM2=IMF(IM1)
         IM3=IMF(IM2)
         IP1=IPF(I)
         IP2=IPF(IP1)

         JM1=JMF(J)
         JM2=JMF(JM1)
         JM3=JMF(JM2)
         JP1=JPF(J)
         JP2=JPF(JP1)

         KM1=KMF(K)
         KM2=KMF(KM1)
         KM3=KMF(KM2)
         KP1=KPF(K)
         KP2=KPF(KP1)

       DALPHIX=0.5*(ALPHI_EXT(IP1,J,K,1)
     &                                   -ALPHI_EXT(IM1,J,K,1))*SSDXF
       DALPHIY=0.5*(ALPHI_EXT(I,JP1,K,1)
     &                                   -ALPHI_EXT(I,JM1,K,1))*SSDYF
       DALPHIZ=0.5*(ALPHI_EXT(I,J,KP1,1)
     &                                   -ALPHI_EXT(I,J,KM1,1))*SSDZF
       DALPHI=DALPHIX**2+DALPHIY**2+DALPHIZ**2

        ASIGN=ALPHIO(I,J,K)/SQRT(ALPHIO(I,J,K)**2+DALPHI*EPS)
!C        !USE ALPHIO. BECAUSE ZERO-INTERFACE SHOULD BE FIXED.
!C        THEN SIGN AROUND INTERAFCE IS FIXED.
!C        !IF INTERFACE IS TOO FLAT, ABOVE eqn CAN HAVE PROBLEM!!
!C         PENG ET AL.(JCP,1999)<- probably fixed the above prob.
!     "Peng's choice of approximation to sign funct. solve the problem of the changing of sign of phi"

      ! X-DIR
        V1=DPHI_X(IM3,J,K)
        V2=DPHI_X(IM2,J,K)
        V3=DPHI_X(IM1,J,K)
        V4=DPHI_X(I,J,K)
        V5=DPHI_X(IP1,J,K)
         CALL HJWENO(V2-V1,V3-V2,V4-V3,V5-V4,WENO_X)
        ADLPHIXM=1./12.*(-V2+7.*V3+7.*V4-V5)-WENO_X
        V1=DPHI_X(IP2,J,K)
        V2=DPHI_X(IP1,J,K)
        V3=DPHI_X(I,J,K)
        V4=DPHI_X(IM1,J,K)
        V5=DPHI_X(IM2,J,K)
         CALL HJWENO(V1-V2,V2-V3,V3-V4,V4-V5,WENO_X)
        ADLPHIXP=1./12.*(-V5+7.*V4+7.*V3-V2)+WENO_X

      !Y-DIR
        V1=DPHI_Y(I,JM3,K)
        V2=DPHI_Y(I,JM2,K)
        V3=DPHI_Y(I,JM1,K)
        V4=DPHI_Y(I,J,K)
        V5=DPHI_Y(I,JP1,K)
         CALL HJWENO(V2-V1,V3-V2,V4-V3,V5-V4,WENO_Y)
        ADLPHIYM=1./12.*(-V2+7.*V3+7.*V4-V5)-WENO_Y 
        V1=DPHI_Y(I,JP2,K)
        V2=DPHI_Y(I,JP1,K)
        V3=DPHI_Y(I,J,K)  
        V4=DPHI_Y(I,JM1,K)
        V5=DPHI_Y(I,JM2,K)   
         CALL HJWENO(V1-V2,V2-V3,V3-V4,V4-V5,WENO_Y)
        ADLPHIYP=1./12.*(-V5+7.*V4+7.*V3-V2)+WENO_Y

      !Z-DIR
        V1=DPHI_Z(I,J,KM3)
        V2=DPHI_Z(I,J,KM2)
        V3=DPHI_Z(I,J,KM1)
        V4=DPHI_Z(I,J,K)
        V5=DPHI_Z(I,J,KP1)
         CALL HJWENO(V2-V1,V3-V2,V4-V3,V5-V4,WENO_Z)
        ADLPHIZM=1./12.*(-V2+7.*V3+7.*V4-V5)-WENO_Z       
        V1=DPHI_Z(I,J,KP2)
        V2=DPHI_Z(I,J,KP1)
        V3=DPHI_Z(I,J,K)
        V4=DPHI_Z(I,J,KM1)
        V5=DPHI_Z(I,J,KM2)
         CALL HJWENO(V1-V2,V2-V3,V3-V4,V4-V5,WENO_Z)
        ADLPHIZP=1./12.*(-V5+7.*V4+7.*V3-V2)+WENO_Z 

       CALL GODUNOV(ASIGN,ADLPHIXP,ADLPHIXM,ALPHI_X)
       CALL GODUNOV(ASIGN,ADLPHIYP,ADLPHIYM,ALPHI_Y)
       CALL GODUNOV(ASIGN,ADLPHIZP,ADLPHIZM,ALPHI_Z)

!C***** EULER METHOD 
       DRHS=ASIGN*(1.-SQRT(ALPHI_X**2+ALPHI_Y**2+ALPHI_Z**2))
       ALPHI_NEW(NN,N)=ALPHI_EXT(I,J,K,1)+REDTL*DRHS

        ENDDO
       ENDDO
!C======================= T-TUBE : 5TH-ORDER WENO========================C  

!C======================= N-TUBE :  1ST-ORDER UPWIND  ===================C
       DO N=NN1+1,NN2
!$OMP PARALLEL DO private(I,J,K,IM1,IP1,JM1,JP1,KM1,KP1)
!$OMP&private(DALPHIX,DALPHIY,DALPHIZ,DALPHI,ASIGN,DRHS)
!$OMP&private(ADLPHIXM,ADLPHIXP,ALPHI_X,ADLPHIYM,ADLPHIYP,ALPHI_Y)
!$OMP&private(ADLPHIZM,ADLPHIZP,ALPHI_Z)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)

         IM1=IMF(I)
         IP1=IPF(I)
         JM1=JMF(J)
         JP1=JPF(J)
         KM1=KMF(K)
         KP1=KPF(K)

       DALPHIX=0.5*(ALPHI_EXT(IP1,J,K,1)
     &                                   -ALPHI_EXT(IM1,J,K,1))*SSDXF
       DALPHIY=0.5*(ALPHI_EXT(I,JP1,K,1)
     &                                   -ALPHI_EXT(I,JM1,K,1))*SSDYF
       DALPHIZ=0.5*(ALPHI_EXT(I,J,KP1,1)
     &                                   -ALPHI_EXT(I,J,KM1,1))*SSDZF
       DALPHI=DALPHIX**2+DALPHIY**2+DALPHIZ**2

        ASIGN=ALPHIO(I,J,K)/SQRT(ALPHIO(I,J,K)**2+DALPHI*EPS)
        ADLPHIXM=( ALPHI_EXT(I,J,K,1)-ALPHI_EXT(IM1,J,K,1) )*SSDXF
        ADLPHIXP=( ALPHI_EXT(IP1,J,K,1)-ALPHI_EXT(I,J,K,1) )*SSDXF
        ADLPHIYM=( ALPHI_EXT(I,J,K,1)-ALPHI_EXT(I,JM1,K,1) )*SSDYF
        ADLPHIYP=( ALPHI_EXT(I,JP1,K,1)-ALPHI_EXT(I,J,K,1) )*SSDYF
        ADLPHIZM=( ALPHI_EXT(I,J,K,1)-ALPHI_EXT(I,J,KM1,1) )*SSDZF
        ADLPHIZP=( ALPHI_EXT(I,J,KP1,1)-ALPHI_EXT(I,J,K,1) )*SSDZF

!C***** FIND RHS
       CALL GODUNOV(ASIGN,ADLPHIXP,ADLPHIXM,ALPHI_X)
       CALL GODUNOV(ASIGN,ADLPHIYP,ADLPHIYM,ALPHI_Y)
       CALL GODUNOV(ASIGN,ADLPHIZP,ADLPHIZM,ALPHI_Z)

!C***** EULER METHOD
       DRHS=ASIGN*(1.-SQRT(ALPHI_X**2+ALPHI_Y**2+ALPHI_Z**2))
       ALPHI_NEW(NN,N)=ALPHI_EXT(I,J,K,1)+REDTL*DRHS

         ENDDO
        ENDDO
C======================= N-TUBE :  1ST-ORDER UPWIND  ===================C

       DO N=1,NN2
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
          I=I_B(NN,N,LLVS)
          J=J_B(NN,N,LLVS)
          K=K_B(NN,N,LLVS)
           ALPHI_EXT(I,J,K,1)=ALPHI_NEW(NN,N)   
        ENDDO
       ENDDO

      RETURN
      END


!C*******************************************************************
      SUBROUTINE GODUNOV(ASIGN,ADLPHIP,ADLPHIM,ADLPHI)
!C*******************************************************************
      IMPLICIT NONE

      REAL ASIGN,ADLPHIP,ADLPHIM,ADLPHI
      REAL APLUS,AMINUS

      !GODUNOV FLUX FUNC.
!c-----A PDE-based fast local level set method - peng et al.(1999)
      !mistake in paper...? i calculate sqrt(max(a,b)+max(c,d))
      APLUS=ASIGN*ADLPHIP
      AMINUS=ASIGN*ADLPHIM
      IF ( (APLUS .GE. 0.) .AND. (AMINUS .GE. 0.) ) THEN  !UPWIND
        ADLPHI=ADLPHIM
      ELSE IF( (APLUS .LE. 0.) .AND. (AMINUS .LE. 0.) ) THEN
        ADLPHI=ADLPHIP
      ELSE IF( (APLUS .GT. 0.) .AND. (AMINUS .LT. 0.) ) THEN  !EXPANSION
        ADLPHI=0.
      ELSE IF( (APLUS .LT. 0.) .AND. (AMINUS .GT. 0.) ) THEN  !SHOCK
        IF (APLUS+AMINUS .GT. 0. ) THEN
          ADLPHI=ADLPHIM
        ELSE
          ADLPHI=ADLPHIP
        ENDIF
      ENDIF

       END

!C*******************************************************************
      SUBROUTINE HJWENO(A,B,C,D,ADLPHI)
!C*******************************************************************
C-----WEIGHTED ENO SCHEMES FOR HAMILTON-JACOBI EQUATIONS-JIANG ET AL.(2000)
      IMPLICIT NONE

       REAL A,B,C,D,ADLPHI

       REAL EPS,S1,S2,S3,A1,A2,A3,W1,W3

      EPS=1.E-6

       S1=13.*(A-B)**2+3.*(A-3*B)**2
       S2=13.*(B-C)**2+3.*(B+C)**2
       S3=13.*(C-D)**2+3.*(3.*C-D)**2
        A1=1./(EPS+S1)**2
        A2=6./(EPS+S2)**2
        A3=3./(EPS+S3)**2
        W1=A1/(A1+A2+A3)
        W3=A3/(A1+A2+A3)

        ADLPHI=1./3.*W1*(A-2.*B+C)+1./6.*(W3-0.5)*(B-2.*C+D)

       RETURN
       END

!C*******************************************************************
      SUBROUTINE CAL_PSIF(INI_PSI,NN1,LLVS,U,V,W,PSI_CN,BUBCOL,bubcol2
     &,vol_tot_ori )
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
       use flow_var ! 2017-10-03
      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE

      REAL PSI_TMP(0:M1,0:M2,0:M3)
      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)
!      REAL psif1(M1L:M1U,M2L:M2U,M3L:M3U)               !   2017-08-19   추가
!      REAL psif2(M1L:M1U,M2L:M2U,M3L:M3U)               !   2017-08-19   추가
!      REAL psif3(M1L:M1U,M2L:M2U,M3L:M3U)               !   2017-08-19   추가
      INTEGER BUBCOL(MLVS)
      INTEGER BUBCOL2(MLVS)                              !   2017-08-20 추가
      REAL D_CRI                                         !   2017-08-20 추가
      real vol_tot_ori(mlvs)                             ! 2017-09-26ㅊㄱ
      INTEGER INI_PSI,NN1,LLVS
      
      INTEGER I,J,K,N,NN
      
      PSIF_BACKGROUND=1. !PSIF_BACKGROUND OPTION 1: water background
      ! IF ( PSIF_BACKGROUND .EQ. 0.) THEN                 !   2017-08-16 주석처리
      ! WRITE(*,*) 'PSIF_BACKGROUND=',PSIF_BACKGROUND      !   2017-08-16 주석처리
      ! ENDIF                                              !   2017-08-16 주석처리

      BUBCOL(LLVS)=0
      BUBCOL2(LLVS)=0

       IF (INI_PSI .EQ. 1) THEN
      !CAL PSIF_MTP
      CALL CAL_PSIF_MTP(1,NN1,LLVS,PSIF_MTP)
!C-----MAKE GLOBAL VOLUME FRACTION
      !below correction should be done in whole domain covering all the bubble volume.
      !local band usually cannot cover all the volume region.
      IF (LLVS .EQ. 1) THEN
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3FM
       DO J=1,N2FM
       DO I=1,N1FM
         PSIF(I,J,K)=0.        
       ENDDO
       ENDDO
       ENDDO
!$OMP PARALLEL DO private(I,J)       
       DO K=1,N3FM
       DO J=1,N2FM
       DO I=1,N1FM
         PSIF(I,J,K)=PSIF_MTP(I,J,K)        
       ENDDO
       ENDDO
       ENDDO

      ELSE
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3FM
       DO J=1,N2FM
       DO I=1,N1FM
          PSIF(I,J,K)=MIN(PSIF(I,J,K),PSIF_MTP(I,J,K)) ! MIN IF WATER-BACKGROUND, MAX IF AIR BACKGROUND
       ENDDO
       ENDDO
       ENDDO
      ENDIF !IF (LLVS .EQ. 1)

      ELSE !IF (INI_PSI .EQ. 1) THEN
      
       ! if (msub .eq. 3) CALL BREAKUP(LLVS,VOL_TOT_ORI) !2017-10-06
       !2018-08-02_CALL BREAKUP(LLVS,VOL_TOT_ORI) !2017-10-06
   
      !CAL PSIF_MTP
      CALL CAL_PSIF_MTP(0,NN1,LLVS,PSIF_MTP)   
      IF (LLVS .EQ. 1) THEN
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3FM
       DO J=1,N2FM
       DO I=1,N1FM
        MASK_GLOBAL(I,J,K)=0   
        MASK_COL(I,J,K)=0    !2017-08-20   
       ENDDO
       ENDDO
       ENDDO
      ENDIF
      
      
!!-----!2017-08-19-- checking how much psif------------------------------
!!       DO N=1,nn1 !      NN1=n_max-3
!!       DO NN=1,NUMA(N,LLVS)
!!       I=I_B(NN,N,LLVS)
!!       J=J_B(NN,N,LLVS)
!!       K=K_B(NN,N,LLVS)
!!       if (llvs .eq. 1) then       
!!       psif1(i,j,k)=psif_mtp(i,j,k)!2017-08-19
!!       else
!!       psif2(i,j,k)=psif_mtp(i,j,k)!2017-08-19
!!       ENdif
!!       ENDDO
!!       ENDDO
!!-----!2017-08-19------------------------------------------------------  

!C-----MAKE GLOBAL VOLUME FRACTION
      !below correction should be done in whole domain covering all the bubble volume.
      !local band usually cannot cover all the volume region.
      !IF (PSIF_BACKGROUND .EQ. 1.) THEN
       DO N=1,NN1 !      NN1=n_max-3
!$OMP PARALLEL DO private(I,J,K)
       DO NN=1,NUMA(N,LLVS)
       I=I_B(NN,N,LLVS)
       J=J_B(NN,N,LLVS)
       K=K_B(NN,N,LLVS)
!WRITE NEW VALUE AT FIRST TIME. FOR DISCARDING PREVIOUS VALUE NEAR THE INTERFACE
        IF (MASK_GLOBAL(I,J,K) .EQ. 0) THEN 
          PSIF(I,J,K)=PSIF_MTP(I,J,K)
          MASK_GLOBAL(I,J,K)=LLVS
!OVERLAP THE PSIF FOR MTP LVS.
        ELSE
         IF ( ABS(PSIF(I,J,K)-0.5) .LT. 0.5 
     &  .AND. ABS(PSIF_MTP(I,J,K)-0.5) .LT. 0.5) THEN !BUB_COLLISION_CRI      
            BUBCOL(MASK_GLOBAL(I,J,K))=1              ! this is when the two bubbles 
            BUBCOL(LLVS)=1                            ! cross each other
           
            BUBCOL2(MASK_COL(I,J,K))=1                
            BUBCOL2(LLVS)=1                           
            CONTACT(LLVS,MASK_COL(I,J,K))=1         
            CONTACT(MASK_COL(I,J,K),LLVS)=1  
         ENDIF
          PSIF(I,J,K)=MIN(PSIF(I,J,K),PSIF_MTP(I,J,K))!psif(i,j,k)+psif_mtp(i,j,k)!MAX(PSIF(I,J,K),PSIF_MTP(I,J,K)) ! MIN IF WATER-BACKGROUND, MAX IF AIR BACKGROUND
        ENDIF       
       ENDDO
       ENDDO
      
!!-----!2017-08-20------------------------------------------------------
       D_CRI=0.5*MIN(SDXF,SDYF,SDZF) 
       DO N=1,NN1 !   
!$OMP PARALLEL DO private(I,J,K)
       DO NN=1,NUMA(N,LLVS)
       I=I_B(NN,N,LLVS)
       J=J_B(NN,N,LLVS)
       K=K_B(NN,N,LLVS)
        IF (MASK_COL(I,J,K) .EQ. 0) THEN 
          ALPHI_COL(I,J,K)=ALPHI_EXT(I,J,K,1)
          MASK_COL(I,J,K)=LLVS
        ELSE
         IF ( (ABS(ALPHI_COL(I,J,K)) .LE. d_cri)       ! this is when the two bubbles
     &  .and. (ABS(ALPHI_COL(I,J,K)) .NE. 0.   )       ! touch each other
     &  .and. (ABS(ALPHI_ext(I,J,K,1)) .NE. 0.   )    
     &  .and. (abs(alphi_ext(i,j,k,1)) .le. d_cri)) THEN !BUB_COLLISION_CRI      
            !write(*,*)ABS(ALPHI_COL(I,J,K)),ABS(ALPHI_COL(I,J,K))
            BUBCOL2(MASK_COL(I,J,K))=1
            BUBCOL2(LLVS)=1
            
            CONTACT(LLVS,MASK_COL(I,J,K))=1         
            CONTACT(MASK_COL(I,J,K),LLVS)=1         
         ENDIF
         ALPHI_COL(I,J,K)=MIN(ALPHI_COL(I,J,K),ALPHI_EXT(I,J,K,1))
        ENDIF

       ENDDO
       ENDDO
       
!!-----!2017-08-20------------------------------------------------------         
        
      !ELSE 
      !  write(*,*) 'psif_background=0, revisit here'
!       DO N=1,NN1
!!$OMP PARALLEL DO private(I,J,K)
!       DO NN=1,NUMA(N,LLVS)
!       I=I_B(NN,N,LLVS)
!       J=J_B(NN,N,LLVS)
!       K=K_B(NN,N,LLVS)
!!WRITE NEW VALUE AT FIRST TIME. FOR DISCARDING PREVIOUS VALUE NEAR THE INTERFACE
!        IF (MASK_GLOBAL(I,J,K) .EQ. 0) THEN 
!          PSIF(I,J,K)=PSIF_MTP(I,J,K)
!          MASK_GLOBAL(I,J,K)=1
!!OVERLAP THE PSIF FOR MTP LVS.
!        ELSE
!        IF (PSIF_MTP(I,J,K) .LT. 0.999) THEN
!          PSIF(I,J,K)=PSIF(I,J,K)+(PSIF_MTP(I,J,K)-1.)
!        ELSE
!          PSIF(I,J,K)=1.
!        ENDIF
!        ENDIF
!       ENDDO
!       ENDDO
      !ENDIF !  IF (PSIF_BACKGROUND .EQ. 1.) THEN

      ENDIF ! IF (INI_PSI .EQ. 1) THEN
   
      !CAL PSI_TMP
      CALL CAL_PSI_TMP(INI_PSI,NN1,LLVS,PSI_TMP,PSIF_MTP,PSI_CN)

      !BUBBLE_POSITION_TRACKING
       ! CALL BUBBLE_TRACKING(LLVS,PSI_TMP,U,V,W)    

      !IF (LLVS .EQ. NLVS) THEN !AT LAST STEP.  !2017-08-30
        if ((llvs .ne. 1) .and. (bubcol2(llvs).eq.1)) then
        write(*,*) '----------LVS CONTACT---------------- '
        write(*,*) 'CONTACT of LVS=',LLVS !, 'bubcollision'
       !elseif ((llvs .ne. 1) .and. (bubcol2(llvs).eq.0)) then     
       !write(*,*) '----------no llvs contact------------'
        endif
      !2017-09-22_MOVED HERE FOR BUBBLE VELOCITY
        
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3FM
       DO J=1,N2FM
       DO I=1,N1FM
          IF (PSIF(I,J,K) .LT. 1d-5) THEN !FOR NUMERICAL ERROR
            PSIF(I,J,K)=0. 
          ELSE IF (PSIF(I,J,K) .GT. 0.99999d0) THEN
            PSIF(I,J,K)=1.
          ENDIF
       ENDDO
       ENDDO
       ENDDO      

      !BOUNDARY CONDITION
      IF ( IPX .EQ. 1 ) THEN
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F    !START FROM ZERO BECAUSE OF EDGE OF DOMAIN.(NEED FOR INTER_PSI)
      DO J=0,N2F
      DO I=0,N1F_BD
        PSIF(-I,J,K)=PSIF(N1FM-I,J,K)
        PSIF(N1F+I,J,K)=PSIF(1+I,J,K)
      ENDDO
      ENDDO
      ENDDO
      ELSE
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F    !START FROM ZERO BECAUSE OF EDGE OF DOMAIN.(NEED FOR INTER_PSI)
      DO J=0,N2F
      DO I=0,N1F_BD
        PSIF(-I,J,K)=PSIF(1,J,K)
        PSIF(N1F+I,J,K)=PSIF(N1FM,J,K)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      IF ( IPY .EQ. 1 ) THEN
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F
      DO J=0,N2F_BD
      DO I=0,N1F
        PSIF(I,-J,K)=PSIF(I,N2FM-J,K)
        PSIF(I,N2F+J,K)=PSIF(I,1+J,K)
      ENDDO
      ENDDO
      ENDDO
      ELSE
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F
      DO J=0,N2F_BD
      DO I=0,N1F
        PSIF(I,-J,K)=PSIF(I,1,K)
        PSIF(I,N2F+J,K)=PSIF(I,N2FM,K)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      IF ( IPZ .EQ. 1 ) THEN
      DO K=0,N3F_BD
!$OMP PARALLEL DO private(I)
      DO J=0,N2F
      DO I=0,N1F
        PSIF(I,J,-K)=PSIF(I,J,N3FM-K)
        PSIF(I,J,N3F+K)=PSIF(I,J,1+K)
      ENDDO
      ENDDO
      ENDDO
      ELSE
      DO K=0,N3F_BD
!$OMP PARALLEL DO private(I)
      DO J=0,N2F
      DO I=0,N1F
        PSIF(I,J,-K)=PSIF(I,J,1)
        PSIF(I,J,N3F+K)=PSIF(I,J,N3FM)
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      !ENDIF !IF (LLVS .EQ. NLVS) THEN !2017-08-30
   
      !CAL CURVATURE*DELTA_FUNCTIO-N USING PSIF_MTP
      CALL CAL_CUR(4,LLVS,PSI_TMP,PSIF_MTP) !1 FOR EVERY STEP REINITIALIZATION
!      CALL CAL_CUR_WHOLE(4,LLVS) !1 FOR EVERY STEP REINITIALIZATION

   

!C=====SAVE
!         if(ini_psi .eq. 0) then
!!         
!      DO N=1,nn1
!        DO NN=1,NUMA(N,LLVS)
!         I=I_B(NN,N,LLVS)
!         J=J_B(NN,N,LLVS)
!         K=K_B(NN,N,LLVS)
!        psif3(i,j,k)=psif(I,J,K)
!        ENDDO
!       ENDDO         
!       
!         if (llvs .eq. 1) then
!         if (llvs .eq. 2) then
!       OPEN(146,FILE='03PSIF.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","psif1","psif2","psif3"'!,"ALPHI_EXT"'
!      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!      DO K=1,N3FM
!      DO J=1,N2FM
!      DO I=1,N1FM
!        WRITE(146,147) XPF(I),YPF(J),ZPF(K),psif1(i,j,k)
!     &   ,psif2(i,j,k),psif3(i,j,k)
!      !&                               ,ALPHI_EXT(i,j,k,llvs)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
! ! 147  FORMAT(5F15.8)
! 147  FORMAT(6F15.8)
!        endif
!        endif
!

      RETURN
      END


!C*******************************************************************
      SUBROUTINE CAL_PSIF_MTP(IOUTER,NN1,LLVS,PSIF_MTP)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE

      INTEGER MASK(0:M1F,0:M2F,0:M3F)
      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)

      INTEGER, DIMENSION (:), ALLOCATABLE :: JC,KC
      
      INTEGER IOUTER,NN1,LLVS
      
      INTEGER I,J,K,N,NN,J_NUM,L,k_num
      REAL D1,D2,D3,E1,E2,E3,F1,F2,F3,F4,F5,EPS,ABSALPHI_EXT
      REAL PSIF_T,PSIF_TMP

!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F
      DO J=0,N2F
      DO I=0,N1F
      PSIF_MTP(I,J,K)=0.
      MASK(I,J,K)=0
      ENDDO
      ENDDO
      ENDDO

       EPS=1.E-8

      DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
!$OMP&private(D1,D2,D3,E1,E2,E3,F1,F2,F3,F4,F5,PSIF_T,ABSALPHI_EXT)
      DO NN=1,NUMA(N,LLVS)
       I=I_B(NN,N,LLVS)
       J=J_B(NN,N,LLVS)
       K=K_B(NN,N,LLVS)

       D1=0.5*ABS(ALPHI_EXT(IPF(I),J,K,1)-ALPHI_EXT(IMF(I),J,K,1))
       D2=0.5*ABS(ALPHI_EXT(I,JPF(J),K,1)-ALPHI_EXT(I,JMF(J),K,1))
       D3=0.5*ABS(ALPHI_EXT(I,J,KPF(K),1)-ALPHI_EXT(I,J,KMF(K),1))
        E1=MAX(D1,D2,D3)
        E3=MIN(D1,D2,D3)
        E2=D1+D2+D3-E1-E3

      IF ( (E1 .GT. EPS).AND.(E2 .GT. EPS).AND.(E3 .GT. EPS) ) THEN
        ABSALPHI_EXT=ABS(ALPHI_EXT(I,J,K,1))
        F1=MAX(0.5*(E1+E2+E3)-ABSALPHI_EXT,0.)
        F2=MAX(0.5*(E1+E2-E3)-ABSALPHI_EXT,0.)        
        F3=MAX(0.5*(E1-E2+E3)-ABSALPHI_EXT,0.)
        F4=MAX(0.5*(-E1+E2+E3)-ABSALPHI_EXT,0.)        
        F5=MAX(0.5*(E1-E2-E3)-ABSALPHI_EXT,0.)
        PSIF_T=(F1**3-F2**3-F3**3-F4**3+F5**3)/(6*E1*E2*E3)
      ELSE IF( (E1 .GT. EPS).AND.(E2 .GT. EPS).AND.(E3 .LE. EPS) ) THEN
        ABSALPHI_EXT=ABS(ALPHI_EXT(I,J,K,1))
        F1=MAX(0.5*(E1+E2+E3)-ABSALPHI_EXT,0.)
        F3=MAX(0.5*(E1-E2+E3)-ABSALPHI_EXT,0.)         
        PSIF_T=(F1**2-F3**2)/(2*E1*E2)
      ELSE IF( (E1 .GT. EPS).AND.(E2 .LE. EPS).AND.(E3 .LE. EPS) ) THEN  
        F1=MAX(0.5*(E1+E2+E3)-ABS(ALPHI_EXT(I,J,K,1)),0.)        
        PSIF_T=F1/E1
      ELSE IF( (E1 .LE. EPS).AND.(E2 .LE. EPS).AND.(E3 .LE. EPS)
     &          .AND.(ALPHI_EXT(I,J,K,1) .NE. 0) ) THEN  
         PSIF_T=0.
      ELSE IF( (E1 .LE. EPS).AND.(E2 .LE. EPS).AND.(E3 .LE. EPS)
     &          .AND.(ALPHI_EXT(I,J,K,1) .EQ. 0) ) THEN  
         PSIF_T=0.5
      ELSE
        WRITE(*,*) 'DANGER - CAL_PSI HAS BUG!!'
        write(*,*) '(i,j,k,llvs)',i,j,k,llvs
        stop
      ENDIF
      IF ( ALPHI_EXT(I,J,K,1) .LE. 0 ) THEN
       PSIF_MTP(I,J,K)=PSIF_T
      ELSE
       PSIF_MTP(I,J,K)=1.-PSIF_T  
      ENDIF

         MASK(I,J,K)=1  

      ENDDO
      ENDDO

       DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
      DO NN=1,NUMA(N,LLVS)
       I=I_B(NN,N,LLVS)
       J=J_B(NN,N,LLVS)
       K=K_B(NN,N,LLVS)
          IF (PSIF_MTP(I,J,K) .LT. 1.d-5) THEN !FOR NUMERICAL ERROR
            PSIF_MTP(I,J,K)=0. 
          ELSE IF (PSIF_MTP(I,J,K) .GT. 0.99999d0) THEN
            PSIF_MTP(I,J,K)=1.
          ENDIF
      ENDDO
      ENDDO

!C======================FAR FROM THE INTERFACE_START====================C
      IF (IOUTER .EQ. 1) THEN
!!
!!!$OMP PARALLEL private(I,J,J_NUM,JC,L,PSIF_TMP)
!!      ALLOCATE(JC(N2F))
!!!$OMP DO
!!      DO K=1,N3FM
!!      DO I=1,N1FM
!!
!!        !I SELECT R-DIR BECAUSE WE ALSO WANT TO SIMULATE N1F=2 AND N3F=2.
!!        !R-DIRECTION EXSIT FOR BOTH TWO CASES.
!!        JC(:)=0
!!        J_NUM=0
!!        !FIND INTERFACE
!!       DO J=1,N2FM
!!
!!        IF(ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,JPF(J),K,1) .LT. 0.)THEN
!!          J_NUM=J_NUM+1
!!          JC(J_NUM)=J
!!        ENDIF       
!!       ENDDO
!!              
!!      IF ( J_NUM .EQ. 0 ) THEN
!!        PSIF_TMP=PSIF_BACKGROUND 
!!        DO J=1,N2FM
!!            IF(MASK(I,J,K).EQ.0) PSIF_MTP(I,J,K)=PSIF_TMP
!!        ENDDO
!!      ELSE
!!        !ASSIGN PSIF-INSIDE
!!       DO L=1,J_NUM
!!          IF( ALPHI_EXT(I,JC(L),K,1) .LT. 0. ) THEN
!!            PSIF_TMP=0.
!!          ELSE
!!            PSIF_TMP=1.
!!          ENDIF
!!
!!         IF( L .EQ. 1 ) THEN  !WE CAN HANDLE WHEN JC(L)=1
!!          IF ( ABS(PSIF_MTP(I,1,K)-0.5) .LT. 0.5 ) THEN
!!           DO J=2,JC(L)             !IF INTERFACE EXIST BETWEEN ALPHI_EXT(I,1,K) AND ALPHI_EXT(I,0,K)
!!            IF(MASK(I,J,K).EQ.0) PSIF_MTP(I,J,K)=PSIF_TMP   !PSIF(I,1,K)=ORIGINAL PSIF(I,1,K) AND CALCULATE FROM J=2
!!           ENDDO
!!          ELSE
!!           DO J=1,JC(L)
!!            IF(MASK(I,J,K).EQ.0) PSIF_MTP(I,J,K)=PSIF_TMP
!!           ENDDO
!!          ENDIF
!!         ELSE
!!          DO J=JC(L-1)+1,JC(L)
!!            IF(MASK(I,J,K).EQ.0) PSIF_MTP(I,J,K)=PSIF_TMP
!!          ENDDO
!!         ENDIF
!!       ENDDO
!!
!!        !ASSIGN PSIF-OUTER
!!       IF (JC(L) .NE. N2FM) THEN
!!          IF( ALPHI_EXT(I,JC(J_NUM)+1,K,1) .LT. 0. ) THEN
!!            PSIF_TMP=0.
!!          ELSE
!!            PSIF_TMP=1.
!!          ENDIF
!!          DO J=JC(J_NUM)+1,N2FM
!!            IF(MASK(I,J,K).EQ.0) PSIF_MTP(I,J,K)=PSIF_TMP
!!          ENDDO        
!!       ENDIF
!!
!!      ENDIF
!!      ENDDO
!!      ENDDO
!!!$OMP END DO
!!      DEALLOCATE(JC)
!!!$OMP END PARALLEL

!$OMP PARALLEL private(I,K,K_NUM,KC,L,PSIF_TMP)
      ALLOCATE(KC(N3F))
!$OMP DO
      DO J=1,N2FM
      DO I=1,N1FM

        !!this is in k direction: I SELECT R-DIR BECAUSE WE ALSO WANT TO SIMULATE N1F=2 AND N3F=2.
        KC(:)=0
        K_NUM=0
        !FIND INTERFACE
       DO K=N3FM,1,-1 !1,N3FM ! the order is changed due to 1-d configuration

        IF (ALPHI_EXT(I,J,K,LLVS)*ALPHI_EXT(I,J,KPF(K),LLVS).LT.0.)THEN
          K_NUM=K_NUM+1
          KC(K_NUM)=K
        ENDIF
       ENDDO

      IF ( K_NUM .EQ. 0 ) THEN
        PSIF_TMP=PSIF_BACKGROUND
        DO K=1,N3FM
            IF(MASK(I,J,K).EQ.0) PSIF_MTP(I,J,K)=PSIF_TMP
        ENDDO
      ELSE
        !ASSIGN PSIF-INSIDE
       DO L=1,K_NUM
          IF( ALPHI_EXT(I,J,KC(L),LLVS) .LT. 0. ) THEN
            PSIF_TMP=0.
          ELSE
            PSIF_TMP=1.
          ENDIF

         IF( L .EQ. 1 ) THEN  !WE CAN HANDEL WHEN JC(L)=1
          IF ( ABS(PSIF_MTP(I,J,1)-0.5) .LT. 0.5 ) THEN
           DO K=2,KC(L)             !IF INTERFACE EXIST BETWEEN ALPHI_EXT(I,1,K) AND ALPHI_EXT(I,0,K)
            IF(MASK(I,J,K).EQ.0) PSIF_MTP(I,J,K)=PSIF_TMP   !PSIF(I,1,K)=ORIGINAL PSIF(I,1,K) AND CALCULATE FROM J=2
           ENDDO
          ELSE
           DO K=1,KC(L)
            IF(MASK(I,J,K).EQ.0) PSIF_MTP(I,J,K)=PSIF_TMP
           ENDDO
          ENDIF
         ELSE
          DO K=KC(L-1),KC(L)
            IF(MASK(I,J,K).EQ.0) PSIF_MTP(I,J,K)=PSIF_TMP
          ENDDO
         ENDIF
       ENDDO

        !ASSIGN PSIF-OUTER
       IF (KC(L) .NE. N3FM) THEN
          IF( ALPHI_EXT(I,J,KC(K_NUM)+1,LLVS) .LT. 0. ) THEN
            PSIF_TMP=0.
          ELSE
            PSIF_TMP=1.
          ENDIF
          DO K=KC(K_NUM)+1,N3FM
            IF(MASK(I,J,K).EQ.0) PSIF_MTP(I,J,K)=PSIF_TMP
          ENDDO        
       ENDIF

      ENDIF
      ENDDO
      ENDDO
!$OMP END DO
      DEALLOCATE(KC)
!$OMP END PARALLEL

      ENDIF !IF (IOUTER .EQ. 1)

!C=======================FAR FROM THE INTERFACE_END=====================C
          ! IF (MSUB .EQ. 1) THEN !(2017-08-14)
!          if (llvs .eq. 2) then
! !C=====SAVE
!        OPEN(146,FILE='0PSIF_MTP1.DAT')
!        WRITE(146,*) 'VARIABLES="X","Y","Z","psif"'
!       WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!       DO K=1,N3FM
!       DO J=1,N2FM
!       DO I=1,N1FM
!         WRITE(146,147) XPF(I),YPF(J),ZPF(K),PSIF_MTP(i,j,k)
!        ! WRITE(146,147) XPF(I),YPF(J),ZPF(K),2.*MASK(i,j,k)
!       ENDDO
!       ENDDO
!       ENDDO
!        CLOSE(146)
!  147  FORMAT(5F15.8)
!        endif

       ! ENDIF     !IF (MSUB .EQ. 1) THEN !(2017-08-14)    

      !BOUNDARY CONDITION
      IF ( IPX .EQ. 1 ) THEN
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F    !START FROM ZERO BECAUSE OF EDGE OF DOMAIN.(NEED FOR INTER_PSI)
      DO J=0,N2F
      DO I=0,N1F_BD
        PSIF_MTP(-I,J,K)=PSIF_MTP(N1FM-I,J,K)
        PSIF_MTP(N1F+I,J,K)=PSIF_MTP(1+I,J,K)
      ENDDO
      ENDDO
      ENDDO
      ELSE
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F    !START FROM ZERO BECAUSE OF EDGE OF DOMAIN.(NEED FOR INTER_PSI)
      DO J=0,N2F
      DO I=0,N1F_BD
        PSIF_MTP(-I,J,K)=PSIF_MTP(1,J,K)
        PSIF_MTP(N1F+I,J,K)=PSIF_MTP(N1FM,J,K)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      IF ( IPY .EQ. 1 ) THEN
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F
      DO J=0,N2F_BD
      DO I=0,N1F
        PSIF_MTP(I,-J,K)=PSIF_MTP(I,N2FM-J,K)
        PSIF_MTP(I,N2F+J,K)=PSIF_MTP(I,1+J,K)
      ENDDO
      ENDDO
      ENDDO
      ELSE
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F
      DO J=0,N2F_BD
      DO I=0,N1F
        PSIF_MTP(I,-J,K)=PSIF_MTP(I,1,K)
        PSIF_MTP(I,N2F+J,K)=PSIF_MTP(I,N2FM,K)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      IF ( IPZ .EQ. 1 ) THEN
      DO K=0,N3F_BD
!$OMP PARALLEL DO private(I)
      DO J=0,N2F
      DO I=0,N1F
        PSIF_MTP(I,J,-K)=PSIF_MTP(I,J,N3FM-K)
        PSIF_MTP(I,J,N3F+K)=PSIF_MTP(I,J,1+K)
      ENDDO
      ENDDO
      ENDDO
      ELSE
      DO K=0,N3F_BD
!$OMP PARALLEL DO private(I)
      DO J=0,N2F
      DO I=0,N1F
        PSIF_MTP(I,J,-K)=PSIF_MTP(I,J,1)
        PSIF_MTP(I,J,N3F+K)=PSIF_MTP(I,J,N3FM)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      
      RETURN
      END

!C*******************************************************************
      SUBROUTINE CAL_PSI_TMP(INI_PSI,NN1,LLVS,PSI_TMP,PSIF_MTP,PSI_CN)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING
      
      IMPLICIT NONE

      REAL PSI_TMP(0:M1,0:M2,0:M3)
      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)

      INTEGER MASK(M1M,M2M,M3M)

      REAL PSI_CN(0:M1,0:M2,0:M3)
      
      INTEGER INI_PSI,NN1,LLVS
      
      INTEGER I,J,K,N,NN,M_BD
      REAL PSI_SUM
      INTEGER I1,I2,J1,J2,K1,K2
      REAL X1,X2,Y1,Y2,Z1,Z2

      IF (INI_PSI.EQ.1) THEN
!$OMP PARALLEL DO private(I,J)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       MASK(I,J,K)=1 !INITIALLY ALL THE FIELD SHOULD BE CALCULATED.
       MASK_BUB(I,J,K)=0 !OUT OF BUBBLE=0
      ENDDO
      ENDDO
      ENDDO

      ELSE !ONLY INTERFACE CHANGES.

!$OMP PARALLEL DO private(I,J)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       MASK(I,J,K)=0
      ENDDO
      ENDDO
      ENDDO

      IF (N3M .EQ. 1) THEN
       M_BD=max(M1F_BD,M2F_BD)
        if (n2m .eq. 1) then 
        m_bd=m1f_bd
        endif
      ELSE
       M_BD=max(M1F_BD,M2F_BD,M3F_BD)
      ENDIF

      DO N=1,NN1-M_BD
!$OMP PARALLEL DO private(I,J,K)
      DO NN=1,NUMA(N,LLVS)
       I=I_B(NN,N,LLVS)
       J=J_B(NN,N,LLVS)
       K=K_B(NN,N,LLVS)
       MASK(ICOUMP_VEL(I),JCOUMP_VEL(J),KCOUMP_VEL(K))=1
      ENDDO
      ENDDO

      ENDIF

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        IF (MASK_BUB(I,J,K) .EQ. LLVS) THEN
          PSI_TMP(I,J,K)=PSI_CN(I,J,K)
        ELSE
         PSI_TMP(I,J,K)=1.
        ENDIF
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO 
!$OMP&private(I,J,X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)
       DO 4000 K=1,N3M
       DO 4000 J=1,N2M
       DO 4000 I=1,N1M

       IF (MASK(I,J,K) .EQ. 1) THEN !CALCULATE ONLY NEAR THE BAND
         X1=X(I)
         X2=X(I+1)
         I1=ICOU2(I)
         I2=ICOU1(I+1)
         Y1=Y(J)
         Y2=Y(J+1)
         J1=JCOU2(J)
         J2=JCOU1(J+1)         
         Z1=Z(K)
         Z2=Z(K+1)
         K1=KCOU2(K)
         K2=KCOU1(K+1)

       CALL INTER_PSI_MTP(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM
     &,PSIF_MTP)

        IF (PSI_SUM .LE. 1.E-3) THEN
           PSI_TMP(I,J,K)=0.
           MASK_BUB(I,J,K)=LLVS
        ELSE IF (PSI_SUM .GE. 0.999) THEN
           PSI_TMP(I,J,K)=1.
           MASK_BUB(I,J,K)=0
        ELSE
           PSI_TMP(I,J,K)=PSI_SUM
           MASK_BUB(I,J,K)=LLVS
        ENDIF

       ENDIF

 4000   CONTINUE

!!$OMP PARALLEL DO private(I,J)
!       DO K=1,N3M
!       DO J=1,N2M
!       DO I=1,N1M
!
!       IF (MASK(I,J,K) .EQ. 1) THEN
!
!       ENDIF
!
!       ENDDO
!       ENDDO
!       ENDDO
 
!        if (llvs .eq. 1) then
!       OPEN(146,FILE='0psi_tmp.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","psi_tmp"'
!      WRITE(146,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!        WRITE(146,143) XP(I),YP(J),ZP(K),psi_tmp(i,j,k)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
! 143  format(4f15.7)
!       endif

      RETURN
      END
!C*******************************************************************
      SUBROUTINE CAL_CUR(N_CUR,LLVS,PSI_TMP,PSIF_MTP)
!C*******************************************************************      
      USE TWO_PHASE_PROPERTY

      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE

      REAL CURF_TMP(M1L:M1U,M2L:M2U,M3L:M3U)
      REAL CURF(M1L:M1U,M2L:M2U,M3L:M3U)

      REAL PSI_TMP(0:M1,0:M2,0:M3)
      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)
      
      INTEGER N_CUR,LLVS
      
      INTEGER I,J,K,IP1,IM1,JP1,JM1,KP1,KM1,N,NN
      REAL GX,GY,GZ,GXX,GYY,GZZ,GXY,GXZ,GYZ
      REAL DPX,DPY,DPZ,AMDPI,DIS,X_BASE,Y_BASE,Z_BASE,XD,YD,ZD
      INTEGER I1,I2,J1,J2,K1,K2
      REAL C00,C10,C01,C11,C0,C1
      REAL CUR

      IF (LLVS .EQ. 1) THEN

      ALLOCATE(SUR_X(M1M,M2M,M3M),SUR_Y(M1M,M2M,M3M),SUR_Z(M1M,M2M,M3M))
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        SUR_X(I,J,K)=0.
        SUR_Y(I,J,K)=0.
        SUR_Z(I,J,K)=0.
       ENDDO
       ENDDO
       ENDDO
    
      ENDIF

      DO N=1,N_CUR+1
!$OMP PARALLEL DO private(I,J,K,IP1,IM1,JP1,JM1,KP1,KM1)
!$OMP&private(GX,GY,GZ,GXX,GYY,GZZ,GXY,GXZ,GYZ)
      DO NN=1,NUMA(N,LLVS)
       I=I_B(NN,N,LLVS)
       J=J_B(NN,N,LLVS)
       K=K_B(NN,N,LLVS)
        IP1=IPF(I)
        IM1=IMF(I)
        JP1=JPF(J)
        JM1=JMF(J)
        KP1=KPF(K)
        KM1=KMF(K)

      GX=0.5*( ALPHI_EXT(IP1,J,K,1)-ALPHI_EXT(IM1,J,K,1) )*SSDXF
      GY=0.5*( ALPHI_EXT(I,JP1,K,1)-ALPHI_EXT(I,JM1,K,1) )*SSDYF
      GZ=0.5*( ALPHI_EXT(I,J,KP1,1)-ALPHI_EXT(I,J,KM1,1) )*SSDZF
      GXX=( ALPHI_EXT(IP1,J,K,1)-2.*ALPHI_EXT(I,J,K,1)
     &                           +ALPHI_EXT(IM1,J,K,1) )*SSDXF*SSDXF
      GYY=( ALPHI_EXT(I,JP1,K,1)-2.*ALPHI_EXT(I,J,K,1)
     &                           +ALPHI_EXT(I,JM1,K,1) )*SSDYF*SSDYF
      GZZ=( ALPHI_EXT(I,J,KP1,1)-2.*ALPHI_EXT(I,J,K,1)
     &                           +ALPHI_EXT(I,J,KM1,1) )*SSDZF*SSDZF
      GXY=0.25*( (ALPHI_EXT(IP1,JP1,K,1)-ALPHI_EXT(IM1,JP1,K,1))
     &          -(ALPHI_EXT(IP1,JM1,K,1)-ALPHI_EXT(IM1,JM1,K,1))
     &                                                  )*SSDXF*SSDYF
      GXZ=0.25*( (ALPHI_EXT(IP1,J,KP1,1)-ALPHI_EXT(IM1,J,KP1,1))
     &          -(ALPHI_EXT(IP1,J,KM1,1)-ALPHI_EXT(IM1,J,KM1,1))
     &                                                  )*SSDXF*SSDZF
      GYZ=0.25*( (ALPHI_EXT(I,JP1,KP1,1)-ALPHI_EXT(I,JP1,KM1,1))
     &          -(ALPHI_EXT(I,JM1,KP1,1)-ALPHI_EXT(I,JM1,KM1,1))
     &                                                  )*SSDYF*SSDZF
      CURF_TMP(I,J,K)=-( 
     &    (GXX*(GY**2+GZ**2)+GYY*(GX**2+GZ**2)+GZZ*(GX**2+GY**2))
     &    -2.*(GXY*GX*GY+GXZ*GX*GZ+GYZ*GY*GZ) )/(GX**2+GY**2+GZ**2)**1.5

!       CURX=0.5*( ANOR(IP1,J,K,1)-ANOR(IM1,J,K,1) )*SSDXF
!       CURY=0.5*( ANOR(I,JP1,K,2)-ANOR(I,JM1,K,2) )*SSDYF  
!       CURZ=0.5*( ANOR(I,J,KP1,3)-ANOR(I,J,KM1,3) )*SSDZF
!       CURF_TMP(I,J,K)=-(CURX+CURY+CURZ)

      ENDDO
      ENDDO

      !!Trilinear interpolation
      DO N=1,N_CUR
!$OMP PARALLEL DO private(I,J,K,XD,YD,ZD,I1,I2,J1,J2,K1,K2)
!$OMP&private(IP1,IM1,JP1,JM1,KP1,KM1)
!$OMP&private(DPX,DPY,DPZ,AMDPI,DIS,X_BASE,Y_BASE,Z_BASE)
!$OMP&private(C00,C10,C01,C11,C0,C1)
      DO NN=1,NUMA(N,LLVS)
       I=I_B(NN,N,LLVS)
       J=J_B(NN,N,LLVS)
       K=K_B(NN,N,LLVS)

         IP1=IPF(I)
         IM1=IMF(I)
         JP1=JPF(J)
         JM1=JMF(J)
         KP1=KPF(K)
         KM1=KMF(K)

         DPX=0.5*(ALPHI_EXT(IP1,J,K,1)-ALPHI_EXT(IM1,J,K,1))*SSDXF
         DPY=0.5*(ALPHI_EXT(I,JP1,K,1)-ALPHI_EXT(I,JM1,K,1))*SSDYF
         DPZ=0.5*(ALPHI_EXT(I,J,KP1,1)-ALPHI_EXT(I,J,KM1,1))*SSDZF

         AMDPI=1./SQRT(DPX**2+DPY**2+DPZ**2)
         DIS=ALPHI_EXT(I,J,K,1)*AMDPI

         X_BASE=XPF(I)-DIS*(DPX*AMDPI)
         Y_BASE=YPF(J)-DIS*(DPY*AMDPI)
         Z_BASE=ZPF(K)-DIS*(DPZ*AMDPI)

       IF ( X_BASE .GE. XPF(I) ) THEN
          XD=(X_BASE-XPF(I))*SSDXF
          I1=I
          I2=IP1
       ELSE
          XD=( X_BASE-(XPF(I)-SDXF) )*SSDXF
!          XD=( X_BASE(I,J,K)-XPF(IMF(I)) )*SSDXF
          I1=IM1
          I2=I
       ENDIF

       IF ( Y_BASE .GE. YPF(J) ) THEN
           YD=(Y_BASE-YPF(J))*SSDYF
           J1=J
           J2=JP1
       ELSE
          YD=( Y_BASE-(YPF(J)-SDYF) )*SSDYF
!          YD=( Y_BASE(I,J,K)-YPF(JMF(J)) )*SSDYF
          J1=JM1
          J2=J
       ENDIF

       IF ( Z_BASE .GE. ZPF(K) ) THEN
           ZD=(Z_BASE-ZPF(K))*SSDZF
           K1=K
           K2=KP1
       ELSE
          ZD=( Z_BASE-(ZPF(K)-SDZF) )*SSDZF
!          ZD=( Z_BASE(I,J,K)-ZPF(KMF(K)) )*SSDZF    !PROBLEM WHEN PERIODIC CONDITION
          K1=KM1
          K2=K
       ENDIF

        C00=CURF_TMP(I1,J1,K1)*(1-XD)+CURF_TMP(I2,J1,K1)*XD
        C10=CURF_TMP(I1,J2,K1)*(1-XD)+CURF_TMP(I2,J2,K1)*XD
        C01=CURF_TMP(I1,J1,K2)*(1-XD)+CURF_TMP(I2,J1,K2)*XD
        C11=CURF_TMP(I1,J2,K2)*(1-XD)+CURF_TMP(I2,J2,K2)*XD  

        C0=C00*(1-YD)+C10*YD
        C1=C01*(1-YD)+C11*YD

        CURF(I,J,K)=(C0*(1.-ZD)+C1*ZD)

      ENDDO
      ENDDO

      !BOUNDARY CONDITION
      IF ( IPX .EQ. 1 ) THEN
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F    !START FROM ZERO BECAUSE OF EDGE OF DOMAIN.(NEED FOR INTER_PSI)
      DO J=0,N2F
      DO I=0,N1F_BD
        CURF(-I,J,K)=CURF(N1FM-I,J,K)
        CURF(N1F+I,J,K)=CURF(1+I,J,K)
      ENDDO
      ENDDO
      ENDDO
      ELSE
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F    !START FROM ZERO BECAUSE OF EDGE OF DOMAIN.(NEED FOR INTER_PSI)
      DO J=0,N2F
      DO I=0,N1F_BD
        CURF(-I,J,K)=CURF(1,J,K)
        CURF(N1F+I,J,K)=CURF(N1FM,J,K)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      IF ( IPY .EQ. 1 ) THEN
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F
      DO J=0,N2F_BD
      DO I=0,N1F
        CURF(I,-J,K)=CURF(I,N2FM-J,K)
        CURF(I,N2F+J,K)=CURF(I,1+J,K)
      ENDDO
      ENDDO
      ENDDO
      ELSE
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F
      DO J=0,N2F_BD
      DO I=0,N1F
        CURF(I,-J,K)=CURF(I,1,K)
        CURF(I,N2F+J,K)=CURF(I,N2FM,K)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      IF ( IPZ .EQ. 1 ) THEN
      DO K=0,N3F_BD
!$OMP PARALLEL DO private(I)
      DO J=0,N2F
      DO I=0,N1F
        CURF(I,J,-K)=CURF(I,J,N3FM-K)
        CURF(I,J,N3F+K)=CURF(I,J,1+K)
      ENDDO
      ENDDO
      ENDDO
      ELSE
      DO K=0,N3F_BD
!$OMP PARALLEL DO private(I)
      DO J=0,N2F
      DO I=0,N1F
        CURF(I,J,-K)=CURF(I,J,1)
        CURF(I,J,N3F+K)=CURF(I,J,N3FM)
      ENDDO
      ENDDO
      ENDDO
      ENDIF

C-----FOR MULTIPLE LVS
!$OMP PARALLEL DO private(I,J,IM1,CUR)
       DO K=1,N3M
       DO J=1,N2M
       DO I=IBG,N1M
         IM1=IMV(I)
      IF (ABS(PSI_TMP(I,J,K)-PSI_TMP(IM1,J,K)) .NE. 0. ) THEN
        CALL CAL_CUR_INTEGRAL(I,J,K,1,CUR,PSIF_MTP,CURF)
        SUR_X(I,J,K)=SUR_X(I,J,K)
     &    +SURF*CUR*(PSI_TMP(I,J,K)-PSI_TMP(IM1,J,K))*VVDX(I)
      ENDIF
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J,JM1,CUR)
       DO K=1,N3M
       DO J=JBG,N2M
         JM1=JMV(J)
       DO I=1,N1M
      IF (ABS(PSI_TMP(I,J,K)-PSI_TMP(I,JM1,K)) .NE. 0. ) THEN
        CALL CAL_CUR_INTEGRAL(I,J,K,2,CUR,PSIF_MTP,CURF)
        SUR_Y(I,J,K)=SUR_Y(I,J,K)
     &    +SURF*CUR*(PSI_TMP(I,J,K)-PSI_TMP(I,JM1,K))*VVDY(J)
      ENDIF
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J,KM1,CUR)
       DO K=KBG,N3M
         KM1=KMV(K)
       DO J=1,N2M
       DO I=1,N1M
      IF (ABS(PSI_TMP(I,J,K)-PSI_TMP(I,J,KM1)) .NE. 0. ) THEN
        CALL CAL_CUR_INTEGRAL(I,J,K,3,CUR,PSIF_MTP,CURF)  
        SUR_Z(I,J,K)=SUR_Z(I,J,K)
     &     +SURF*CUR*(PSI_TMP(I,J,K)-PSI_TMP(I,J,KM1))*VVDZ(K)
      ENDIF
       ENDDO
       ENDDO
       ENDDO

!C=====SAVE(BASE)
!      OPEN(60,FILE='0BASE.DAT')
!       WRITE(60,*) 'VARIABLES="I","J","K","X_BASE","y_BASE","z_BASE"
!     &,"XPF","YPF","ZPF","X_BASE-XPF","y_BASE-YPF","z_BASE-ZPF"'
!       DO NN=1,NUMA(1)
!              I=I_B(NN,1)
!             J=J_B(NN,1)
!             K=K_B(NN,1)
!      WRITE(60,61) I,J,K,X_BASE(I,J,K),Y_BASE(I,J,K),Z_BASE(I,J,K)
!     &                ,XPF(I),YPF(J),ZPF(K),X_BASE(I,J,K)-XPF(I)
!     &                ,y_BASE(I,J,K)-YPF(J),Z_BASE(I,J,K)-ZPF(K)
!       ENDDO
!      CLOSE(60)
!  61   FORMAT(3I8,9es15.5)
!C=====SAVE
!       OPEN(146,FILE='0curf_tmp.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","curf_tmp"'
!      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!      DO K=1,N3FM
!      DO J=1,N2FM
!      DO I=1,N1FM
!        WRITE(146,147) XPF(I),YPF(J),ZPF(K),curf_tmp(i,j,k)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
! 147  FORMAT(4F15.8)
!C=====SAVE
!       OPEN(146,FILE='0curf_contour.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","curf"'
!      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!      DO K=1,N3FM
!      DO J=1,N2FM
!      DO I=1,N1FM
!        WRITE(146,143) XPF(I),YPF(J),ZPF(K),curf(i,j,k)
!!        write(146,148) i,j,k,curf(i,j,k)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
! 143  format(4f15.7)
! 148  format(3i5,f15.5)

!      if (llvs .eq. nlvs) then
!C=====SAVE
!       OPEN(146,FILE='0sur_contour.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","sur1","sur2","sur3"'
!      WRITE(146,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!        WRITE(146,148) XPF(I),YPF(J),ZPF(K)
!     &            ,sur_x(i,j,k),sur_y(i,j,k),sur_z(i,j,k)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
! 148  format(6f15.5)
!          endif

      RETURN
      END


! =====================================================================
      SUBROUTINE CAL_CUR_INTEGRAL(I,J,K,IDIR,CUR,PSIF_MTP,CURF)
! =====================================================================
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_COUPLING
      
      IMPLICIT NONE

      REAL CURF(M1L:M1U,M2L:M2U,M3L:M3U)

      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)

      INTEGER I,J,K,IDIR
      REAL CUR

      INTEGER I1,I2,J1,J2,K1,K2
      REAL X1,X2,Y1,Y2,Z1,Z2
      REAL VOLF
      
C-----X-DIR      
      IF (IDIR .EQ. 1) THEN
        X1=XP(I-1)
        X2=XP(I)
        Y1=Y(J)
        Y2=Y(J+1)
        Z1=Z(K)
        Z2=Z(K+1)

        I1=ICOUMP2(I-1)
        I2=ICOUMP1(I)
        J1=JCOU2(J)
        J2=JCOU1(J+1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

         CALL INTER_CUR(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,CUR,VOLF
     &,PSIF_MTP,CURF)
       IF ( VOLF .EQ. 0. ) THEN
!        WRITE(*,*) 'CURX DON`T HAVE INTERFACE ALGORITM HAS ACTIVATED!!'
        X1=X(I-1)
        X2=X(I+1)
        Y1=Y(J)
        Y2=Y(J+1)
        Z1=Z(K)
        Z2=Z(K+1)

        I1=ICOU2(I-1)
        I2=ICOU1(I+1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

         CALL INTER_CUR(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,CUR,VOLF
     &,PSIF_MTP,CURF)
       ENDIF

C-----Y-DIR              
        ELSE IF ( IDIR .EQ. 2) THEN

        X1=X(I)
        X2=X(I+1)
        Y1=YP(J-1)
        Y2=YP(J)
        Z1=Z(K)
        Z2=Z(K+1)

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        J1=JCOUMP2(J-1)
        J2=JCOUMP1(J)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

         CALL INTER_CUR(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,CUR,VOLF
     &,PSIF_MTP,CURF)

       IF ( VOLF .EQ. 0. ) THEN
!        WRITE(*,*) 'CURY DON`T HAVE INTERFACE ALGORITM HAS ACTIVATED!!'
        X1=X(I)
        X2=X(I+1)
        Y1=Y(J-1)
        Y2=Y(J+1)
        Z1=Z(K)
        Z2=Z(K+1)

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        J1=JCOU2(J-1)
        J2=JCOU1(J+1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

         CALL INTER_CUR(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,CUR,VOLF
     &,PSIF_MTP,CURF)
       ENDIF

C-----Z-DIR      
      ELSE IF ( IDIR .EQ. 3) THEN

        X1=X(I)
        X2=X(I+1)
        Y1=Y(J)
        Y2=Y(J+1)
        Z1=ZP(K-1)
        Z2=ZP(K)

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)
        K1=KCOUMP2(K-1)
        K2=KCOUMP1(K)

         CALL INTER_CUR(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,CUR,VOLF
     &,PSIF_MTP,CURF)

       IF ( VOLF .EQ. 0. ) THEN
!        WRITE(*,*) 'CURZ DON`T HAVE INTERFACE ALGORITM HAS ACTIVATED!!'
        X1=X(I)
        X2=X(I+1)
        Y1=Y(J)
        Y2=Y(J+1)
        Z1=Z(K-1)
        Z2=Z(K+1)         

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)
        K1=KCOU2(K-1)
        K2=KCOU1(K+1)

         CALL INTER_CUR(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,CUR,VOLF
     &,PSIF_MTP,CURF)
       ENDIF

      ENDIF
   
       RETURN
       END

!C*******************************************************************
      SUBROUTINE INTER_CUR(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,CUR,VOLF
     &,PSIF_MTP,CURF)
!C*******************************************************************
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE

      REAL CURF(M1L:M1U,M2L:M2U,M3L:M3U)

      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)

      INTEGER I1,I2,J1,J2,K1,K2
      REAL X1,X2,Y1,Y2,Z1,Z2
      REAL CUR,VOLF
      
      INTEGER II,JJ,KK
      REAL C1,C2,CUR_SUM
      REAL AI_CUT,AJ_CUT

!       C1=0. 
!       C2=1.-C1
!
!       CUR=0.
!       VOLF=0.
!
!       DO KK=K1,K2
!       DO JJ=J1,J2
!       DO II=I1,I2
!      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
!     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
!        CUR=CUR+CURF(II,JJ,KK)
!        VOLF=VOLF+1.
!      ENDIF
!       ENDDO
!       ENDDO
!       ENDDO
!  
!       IF ( VOLF .NE. 0.) THEN
!         CUR=CUR/VOLF
!        ELSE
!          CUR=0.
!        ENDIF

       C1=0. 
       C2=1.-C1

       CUR_SUM=0.
       VOLF=0.

      DO II=I1,I2
       !IIIIIIII0
        IF ( (XF(II) .LE. X1).AND. (XF(II+1) .GE. X2 ) ) THEN
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)
        VOLF=VOLF+1.
      ENDIF
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)
        VOLF=VOLF+1.
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ1
        ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
      ENDIF
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)
        VOLF=VOLF+1.
      ENDIF
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)
        VOLF=VOLF+1.
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
      ENDIF
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ4
        ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 40    CONTINUE
       ENDDO

       !IIIIIIII1
        ELSE IF ( (XF(II) .LT. X1).AND. (XF(II+1) .GT. X1 ) ) THEN  
            AI_CUT=(XF(II+1)-X1)*SSDXF
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT
        VOLF=VOLF+AI_CUT
      ENDIF
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT
        VOLF=VOLF+AI_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ1
        ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*AI_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT*AI_CUT
      ENDIF
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*AI_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT*AI_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*AI_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT
        VOLF=VOLF+AI_CUT
      ENDIF
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT
        VOLF=VOLF+AI_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*AI_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT*AI_CUT
      ENDIF
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*AI_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT*AI_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*AI_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ4
        ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 50    CONTINUE
       ENDDO

       !IIIIIIII2
       ELSE IF( (XF(II) .GE. X1).AND.(XF(II+1) .LE. X2) )THEN
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)
        VOLF=VOLF+1.
      ENDIF
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)
        VOLF=VOLF+1.
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF
         ENDDO
        !JJJJJJJJJ1
        ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
      ENDIF
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN 
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)
        VOLF=VOLF+1.
      ENDIF
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)
        VOLF=VOLF+1.
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
      ENDIF
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF
         ENDDO
        !JJJJJJJJJ4
        ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 60    CONTINUE
       ENDDO

       !IIIIIIII3
       ELSE IF( (XF(II) .LT. X2) .AND.((XF(II+1) .GT. X2)) )THEN
            AI_CUT=(X2-XF(II))*SSDXF
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT
        VOLF=VOLF+AI_CUT
      ENDIF
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT
        VOLF=VOLF+AI_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ1
        ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*AI_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT*AI_CUT
      ENDIF
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*AI_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT*AI_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*AI_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT
        VOLF=VOLF+AI_CUT
      ENDIF
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT
        VOLF=VOLF+AI_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*AI_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT*AI_CUT
      ENDIF
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
      ENDIF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+AJ_CUT*AI_CUT*CURF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT*AI_CUT
      ENDIF
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
      IF ( PSIF_MTP(II,JJ,KK) .GT. C1 .AND. PSIF_MTP(II,JJ,KK) .LT. C2 
     &                              .AND. CURF(II,JJ,KK) .NE. 0. ) THEN
        CUR_SUM=CUR_SUM+CURF(II,JJ,KK)*AJ_CUT*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*AI_CUT*(Z2-ZF(KK))*SSDZF
      ENDIF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ4
        ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 70    CONTINUE
       ENDDO

       !IIIIIIII4
       ELSE IF( XF(II) .GE. X2 ) THEN
         GOTO 101
       ENDIF
 100   ENDDO

  101   IF ( VOLF .NE. 0.) THEN
         CUR=CUR_SUM/VOLF
        ELSE
         CUR=0.
        ENDIF

       RETURN
       END
   

!C*******************************************************************
      SUBROUTINE INTER_PSI_MTP(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2
     &,PSI_SUM,PSIF_MTP)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE

      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)

      INTEGER I1,I2,J1,J2,K1,K2
      REAL X1,X2,Y1,Y2,Z1,Z2
      REAL PSI_SUM
      
      INTEGER II,JJ,KK
      REAL VOLF,AI_CUT,AJ_CUT
      
!         !ASSUMING TWO GRIDS ARE ALWAYS ALIGN.
!          PSI_SUM=0.
!          VOLF=0.
!          DO KK=K1,K2
!          DO JJ=J1,J2
!          DO II=I1,I2
!              PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)
!              VOLF=VOLF+1.
!          ENDDO
!          ENDDO
!          ENDDO
!           PSI_SUM=PSI_SUM/VOLF

          PSI_SUM=0.
          VOLF=0.

      DO II=I1,I2
       !IIIIIIII0
        IF ( (XF(II) .LE. X1).AND. (XF(II+1) .GE. X2 ) ) THEN  
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ1
        ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ4
        ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 40    CONTINUE
       ENDDO

       !IIIIIIII1
        ELSE IF ( (XF(II) .LT. X1).AND. (XF(II+1) .GT. X1 ) ) THEN  
            AI_CUT=(XF(II+1)-X1)*SSDXF
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ1
        ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 50 
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM
     &+PSIF_MTP(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM
     &+PSIF_MTP(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM
     &+PSIF_MTP(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM
     &+PSIF_MTP(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ4
        ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 50    CONTINUE
       ENDDO

       !IIIIIIII2
       ELSE IF( (XF(II) .GE. X1).AND.(XF(II+1) .LE. X2) )THEN
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF
         ENDDO
        !JJJJJJJJJ1
       ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF     
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF  
       ENDDO

        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF
         ENDDO
        !JJJJJJJJJ4
           ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 60    CONTINUE
       ENDDO

       !IIIIIIII3
       ELSE IF( (XF(II) .LT. X2) .AND.((XF(II+1) .GT. X2)) )THEN
            AI_CUT=(X2-XF(II))*SSDXF
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ1
       ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM
     &+PSIF_MTP(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM
     &+PSIF_MTP(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO        
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF_MTP(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF_MTP(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM
     &+PSIF_MTP(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF_MTP(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM
     &+PSIF_MTP(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ4
           ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 70    CONTINUE
       ENDDO

       !IIIIIIII4
       ELSE IF( XF(II) .GE. X2 ) THEN
         GOTO 101
       ENDIF
 100   ENDDO

 101   CONTINUE

         PSI_SUM=PSI_SUM/VOLF

       RETURN
       END

!C*******************************************************************
      SUBROUTINE GRID_COUPLING(INI_PSI,PSI_XN,PSI_YN,PSI_ZN,PSI_CN)
!C*******************************************************************  
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING
      
      IMPLICIT NONE

      REAL PSI_XN(0:M1,0:M2,0:M3)
      REAL PSI_YN(0:M1,0:M2,0:M3)
      REAL PSI_ZN(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)

      INTEGER MASK(M1M,M2M,M3M)
      
      INTEGER INI_PSI
      
      INTEGER I,J,K
      
      INTEGER I1,I2,J1,J2,K1,K2
      REAL X1,X2,Y1,Y2,Z1,Z2
      REAL PSI_SUM

      IF (INI_PSI.EQ.1) THEN
!$OMP PARALLEL DO private(I,J)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       MASK(I,J,K)=1. !INITIALLY ALL THE FIELD SHOULD BE CALCULATED.
      ENDDO
      ENDDO
      ENDDO

      ELSE !ONLY INTERFACE CHANGES.
!$OMP PARALLEL DO private(I,J)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       MASK(I,J,K)=0
      ENDDO
      ENDDO
      ENDDO

!$OMP PARALLEL DO private(I,J)
      DO K=1,N3FM
      DO J=1,N2FM
      DO I=1,N1FM
       IF ( MASK_GLOBAL(I,J,K) .NE. 0) THEN
        MASK(ICOUMP_VEL(I),JCOUMP_VEL(J),KCOUMP_VEL(K))=1
       ENDIF
      ENDDO
      ENDDO
      ENDDO

      ENDIF

!!!C-----X-DIR
!!!$OMP PARALLEL DO 
!!!$OMP&private(I,J,X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)
!!       DO 1000 K=1,N3M
!!         Z1=Z(K)
!!         Z2=Z(K+1)
!!         K1=KCOU2(K)
!!         K2=KCOU1(K+1)
!!       DO 1000 J=1,N2M
!!         Y1=Y(J)
!!         Y2=Y(J+1)
!!         J1=JCOU2(J)
!!         J2=JCOU1(J+1)
!!       DO 1000 I=IBG,N1M
!!         X1=XP(I-1)
!!         X2=XP(I)
!!         I1=ICOUMP2(I-1)
!!         I2=ICOUMP1(I)
!!       IF (MASK(I,J,K).EQ.1) THEN
!!       CALL INTER_PSI(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)
!!        PSI_XN(I,J,K)=PSI_SUM
!!       ENDIF
!! 1000  CONTINUE
!!
!!!$OMP PARALLEL DO private(I,J)
!!       DO 1001 K=1,N3M
!!       DO 1001 J=1,N2M
!!       DO 1001 I=IBG,N1M
!!        IF (PSI_XN(I,J,K) .LT. 1.E-3) THEN  !FOR NUMERICAL ERROR
!!           PSI_XN(I,J,K)=0.
!!        ELSE IF (PSI_XN(I,J,K) .GT. 0.999) THEN
!!           PSI_XN(I,J,K)=1.
!!        ENDIF
!! 1001  CONTINUE
!!
!!!C-----Y-DIR
!!!$OMP PARALLEL DO 
!!!$OMP&private(I,J,X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)
!!       DO 2000 K=1,N3M
!!         Z1=Z(K)
!!         Z2=Z(K+1)
!!         K1=KCOU2(K)
!!         K2=KCOU1(K+1)
!!       DO 2000 J=JBG,N2M
!!         Y1=YP(J-1)
!!         Y2=YP(J)
!!         J1=JCOUMP2(J-1)
!!         J2=JCOUMP1(J)
!!       DO 2000 I=1,N1M
!!         X1=X(I)
!!         X2=X(I+1)
!!         I1=ICOU2(I)
!!         I2=ICOU1(I+1)
!!       IF (MASK(I,J,K).EQ.1) THEN
!!       CALL INTER_PSI(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)
!!         PSI_YN(I,J,K)=PSI_SUM
!!       ENDIF
!! 2000   CONTINUE
!!
!!!$OMP PARALLEL DO private(I,J)
!!       DO 2001 K=1,N3M
!!       DO 2001 J=JBG,N2M
!!       DO 2001 I=1,N1M
!!        IF (PSI_YN(I,J,K) .LT. 1.E-3) THEN  !FOR NUMERICAL ERROR
!!           PSI_YN(I,J,K)=0.
!!        ELSE IF (PSI_YN(I,J,K) .GT. 0.999) THEN
!!           PSI_YN(I,J,K)=1.
!!        ENDIF
!! 2001  CONTINUE
!!
!C-----Z-DIR
!$OMP PARALLEL DO 
!$OMP&private(I,J,X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)
       DO 3000 K=KBG,N3M
         Z1=ZP(K-1)
         Z2=ZP(K)
         K1=KCOUMP2(K-1)
         K2=KCOUMP1(K)
       DO 3000 J=1,N2M
         Y1=Y(J)
         Y2=Y(J+1)
         J1=JCOU2(J)
         J2=JCOU1(J+1)
       DO 3000 I=1,N1M
         X1=X(I)
         X2=X(I+1)
         I1=ICOU2(I)
         I2=ICOU1(I+1)
       IF (MASK(I,J,K).EQ.1) THEN
       CALL INTER_PSI(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)
         PSI_ZN(I,J,K)=PSI_SUM
       ENDIF
 3000  CONTINUE

!$OMP PARALLEL DO private(I,J)
       DO 3001 K=1,N3M
       DO 3001 J=JBG,N2M
       DO 3001 I=1,N1M
        IF (PSI_ZN(I,J,K) .LT. 1.E-3) THEN  !FOR NUMERICAL ERROR
           PSI_ZN(I,J,K)=0.
        ELSE IF (PSI_ZN(I,J,K) .GT. 0.999) THEN
           PSI_ZN(I,J,K)=1.
        ENDIF
 3001  CONTINUE
 
!C-----CAL PSI_CN AT CELL CENTER
!$OMP PARALLEL DO 
!$OMP&private(I,J,X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)
       DO 4000 K=1,N3M
         Z1=Z(K)
         Z2=Z(K+1)
         K1=KCOU2(K)
         K2=KCOU1(K+1)
       DO 4000 J=1,N2M
         Y1=Y(J)
         Y2=Y(J+1)
         J1=JCOU2(J)
         J2=JCOU1(J+1)
       DO 4000 I=1,N1M
         X1=X(I)
         X2=X(I+1)
         I1=ICOU2(I)
         I2=ICOU1(I+1)
       IF (MASK(I,J,K).EQ.1) THEN
       CALL INTER_PSI(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)
         PSI_CN(I,J,K)=PSI_SUM
       ENDIF
 4000   CONTINUE

!$OMP PARALLEL DO private(I,J)
       DO 4001 K=1,N3M
       DO 4001 J=1,N2M
       DO 4001 I=1,N1M
        IF (PSI_CN(I,J,K) .LT. 1.E-3) THEN  !FOR NUMERICAL ERROR
           PSI_CN(I,J,K)=0.
        ELSE IF (PSI_CN(I,J,K) .GT. 0.999) THEN
           PSI_CN(I,J,K)=1.
        ENDIF
 4001  CONTINUE
 
!C-----BOUNDARY CONDITION FOR DEVIDE REGION
      IF ( IPX .NE. 1 ) THEN
!$OMP PARALLEL DO private(J)
      DO K=1,N3M           !ASSUMING NO-INTERFACE DOMAIN-BOUNDARY CELL
      DO J=1,N2M
       PSI_CN(0,J,K)=PSI_CN(1,J,K)
       PSI_CN(N1,J,K)=PSI_CN(N1M,J,K)
      ENDDO
      ENDDO
      ENDIF
      IF ( IPY .NE. 1 ) THEN
!$OMP PARALLEL DO private(I)
      DO K=1,N3M
      DO I=1,N1M
       PSI_CN(I,0,K)=PSI_CN(I,1,K)
       PSI_CN(I,N2,K)=PSI_CN(I,N2M,K)
      ENDDO
      ENDDO
      ENDIF
      IF ( IPZ .NE. 1 ) THEN
!$OMP PARALLEL DO private(I)
      DO J=1,N2M
      DO I=1,N1M                    
       PSI_CN(I,J,0)=PSI_CN(I,J,1)
       PSI_CN(I,J,N3)=PSI_CN(I,J,N3M)
      ENDDO
      ENDDO
      ENDIF

!C=====SAVE
!       if(msub .eq. 1) then
!       OPEN(142,FILE='0PSI.DAT')
!       WRITE(142,*) 'VARIABLES="X","Y","Z","PSIX","PSIY","PSIZ","PSIC"'
!      WRITE(142,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
!       DO K=1,N3M
!       DO J=1,N2M
!       DO I=1,N1M
!         psi1=0.5*(PSI_XN(I,J,K)+PSI_XN(IPV(I),J,K))
!         psi2=0.5*(PSI_YN(I,J,K)+PSI_YN(I,JPV(J),K))
!         psi3=0.5*(PSI_ZN(I,J,K)+PSI_ZN(I,J,Kpv(k)))
!         psi4=PSI_CN(I,J,K)         
!!         IF(i.ne.1.and.I.NE.N1M) psi1=PSI_XN(I,J,K)
!!         IF(j.ne.1.and.J.NE.N2M) psi2=PSI_YN(I,J,K)
!!         psi3=PSI_ZN(I,J,K)
!!         psi4=PSI_CN(I,J,K)         
!       WRITE(142,148) XP(I),YP(J),ZP(K),PSI1,PSI2,PSI3,PSI4  
!        ENDDO
!        ENDDO
!        ENDDO
!       CLOSE(142)
! 148  FORMAT(7F15.8)
!       endif

        RETURN
        END

        

!C*******************************************************************
      SUBROUTINE INTER_PSI(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      implicit none
      INTEGER I1,I2,J1,J2,K1,K2
      REAL X1,X2,Y1,Y2,Z1,Z2
      REAL PSI_SUM
      
      INTEGER II,JJ,KK
      REAL VOLF,AI_CUT,AJ_CUT

!          PSI_SUM=0.
!          VOLF=0.
!          DO KK=K1,K2
!          DO JJ=J1,J2
!          DO II=I1,I2
!              PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
!              VOLF=VOLF+1.
!          ENDDO
!          ENDDO
!          ENDDO
!           PSI_SUM=PSI_SUM/VOLF

          PSI_SUM=0.
          VOLF=0.

      DO II=I1,I2
       !IIIIIIII0
        IF ( (XF(II) .LE. X1).AND. (XF(II+1) .GE. X2 ) ) THEN  
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ1
        ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 40
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 40
           ENDIF
         ENDDO
        !JJJJJJJJJ4
        ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 40    CONTINUE
       ENDDO

       !IIIIIIII1
        ELSE IF ( (XF(II) .LT. X1).AND. (XF(II+1) .GT. X1 ) ) THEN  
            AI_CUT=(XF(II+1)-X1)*SSDXF
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ1
        ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 50 
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 50
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 50
           ENDIF
         ENDDO
        !JJJJJJJJJ4
        ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 50    CONTINUE
       ENDDO

       !IIIIIIII2
       ELSE IF( (XF(II) .GE. X1).AND.(XF(II+1) .LE. X2) )THEN
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF
         ENDDO
        !JJJJJJJJJ1
       ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF     
         ENDDO
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
              VOLF=VOLF+1.
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)
        VOLF=VOLF+1.
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF  
       ENDDO

        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AJ_CUT
              GOTO 60
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 60
           ENDIF
         ENDDO
        !JJJJJJJJJ4
           ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 60    CONTINUE
       ENDDO

       !IIIIIIII3
       ELSE IF( (XF(II) .LT. X2) .AND.((XF(II+1) .GT. X2)) )THEN
            AI_CUT=(X2-XF(II))*SSDXF
       DO JJ=J1,J2
        !JJJJJJJJJ0
        IF ( YF(JJ) .LE. Y1 .AND. YF(JJ+1) .GE. Y2 ) THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ1
       ELSE IF ( (YF(JJ) .LT. Y1).AND. (YF(JJ+1) .GT. Y1 ) ) THEN
            AJ_CUT=(YF(JJ+1)-Y1)*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO        
        !JJJJJJJJJ2
       ELSE IF( (YF(JJ) .GE. Y1).AND.(YF(JJ+1) .LE. Y2) )THEN
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ3
       ELSE IF( (YF(JJ) .LT. Y2) .AND.((YF(JJ+1) .GT. Y2)) )THEN
            AJ_CUT=(Y2-YF(JJ))*SSDYF
         DO KK=K1,K2
           IF ( ZF(KK) .LE. Z1 .AND. ZF(KK+1) .GE. Z2 ) THEN
              PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
              VOLF=VOLF+AI_CUT*AJ_CUT
              GOTO 70
           ELSE IF ( (ZF(KK) .LT. Z1).AND. (ZF(KK+1) .GT. Z1 ) ) THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(ZF(KK+1)-Z1)*SSDZF
           ELSE IF( (ZF(KK) .GE. Z1).AND.(ZF(KK+1) .LE. Z2) )THEN
        PSI_SUM=PSI_SUM+AI_CUT*AJ_CUT*PSIF(II,JJ,KK)
        VOLF=VOLF+AI_CUT*AJ_CUT
           ELSE IF( (ZF(KK) .LT. Z2) .AND.((ZF(KK+1) .GT. Z2)) )THEN
        PSI_SUM=PSI_SUM+PSIF(II,JJ,KK)*AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
        VOLF=VOLF+AI_CUT*AJ_CUT*(Z2-ZF(KK))*SSDZF
           ELSE IF( ZF(KK) .GE. Z2 ) THEN
         GOTO 70
           ENDIF
         ENDDO
        !JJJJJJJJJ4
           ELSE IF( YF(JJ) .GE. Y2 ) THEN
         GOTO 100
       ENDIF
 70    CONTINUE
       ENDDO

       !IIIIIIII4
       ELSE IF( XF(II) .GE. X2 ) THEN
         GOTO 101
       ENDIF
 100   ENDDO

 101   CONTINUE

         PSI_SUM=PSI_SUM/VOLF

       RETURN
       END
       
!C*******************************************************************
      SUBROUTINE VEL_INTP(IVEL_UPDATE,NN2,LLVS,UF,VF,WF,U,V,W,UT,VT,WT)
!C*******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING
      
      USE PHASE_BOILING
      USE TWO_PHASE_PROPERTY     
      IMPLICIT NONE

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 UT(-1:M1+2,-1:M2+2,-1:M3+2),VT(-1:M1+2,-1:M2+2,-1:M3+2),
     &                                    WT(-1:M1+2,-1:M2+2,-1:M3+2)
      REAL*8 UF(M1F,M2F,M3F),VF(M1F,M2F,M3F),WF(M1F,M2F,M3F) 
      
      INTEGER*8 IVEL_UPDATE,NN2,LLVS,IM1,IP1,JM1,JP1,KM1,KP1
      
      INTEGER*8  I,J,K,N,NN
      REAL*8     SM_GRID
      REAL*8     U_INTP,V_INTP,W_INTP
      REAL*8     GX,GY,GZ,AGI
      REAL*8     NORM_X,NORM_Y,NORM_Z
      integer*8   nn1
      
      nn1=n_max-3
      
        IF ( M .EQ. 1 .AND. MSUB .EQ. 1 .AND. LLVS .EQ. 1) THEN
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
         UO(I,J,K)=U(I,J,K)
         VO(I,J,K)=V(I,J,K)
         WO(I,J,K)=W(I,J,K)
       ENDDO
       ENDDO
       ENDDO
        ENDIF


!C========================VEL_TRANSPORT_START===========================C
      IF ( IVEL_UPDATE .EQ. 1 .OR. IVEL_UPDATE .EQ. 2 ) THEN
      !UO,VO,WO : N_STEP VEL.
      !U,V,W : (N+1)_STEP VEL.
      !UT,VT,WT : LVS TRANSPORT VEL.

!$OMP PARALLEL DO private(I,J)
      DO K=0,N3
      DO J=0,N2
      DO I=0,N1
      UT(I,J,K)=0.5*(U(I,J,K)+UO(I,J,K))
      VT(I,J,K)=0.5*(V(I,J,K)+VO(I,J,K))
      WT(I,J,K)=0.5*(W(I,J,K)+WO(I,J,K))
      ENDDO
      ENDDO
      ENDDO

      IF ( IVEL_UPDATE .EQ. 2 ) THEN  !VEL_SAVE

!$OMP PARALLEL DO private(I,J)
      DO K=0,N3
      DO J=0,N2
      DO I=0,N1
        UO(I,J,K)=U(I,J,K)
        VO(I,J,K)=V(I,J,K)
        WO(I,J,K)=W(I,J,K)
      ENDDO
      ENDDO
      ENDDO

      ENDIF !IF ( IVEL_UPDATE .EQ. 2 ) THEN

      ENDIF !IF ( IVEL_UPDATE .EQ. 1 .OR. IVEL_UPDATE .EQ. 2 ) THEN
!C==========================VEL_TRANSPORT_END===========================C

!C=========================VELOCITY_INTP_START==========================C
      !I,J,K : ACTIVATED REFINED GRID
      !II,JJ,KK : FLOW_SOLVER GIRD MATCHED TO ACTIVATED REFINED GRID
      !UF,VF,WF : INTERPOLATED REFINED VEL.

        IF ( N3M .EQ. 1 ) THEN
         SM_GRID=MIN(SDXF,SDYF)
          if (n2m .eq. 1) sm_grid=sdxf
       ELSE
         SM_GRID=MIN(SDXF,SDYF,SDZF)
       ENDIF

       DO N=1,NN2 !NN2 EQUAL N_MAX-1
!$OMP PARALLEL DO private(I,J,K,U_INTP,V_INTP,W_INTP)
!$OMP&private(IP1,IM1,JP1,JM1,KP1,KM1,GX,GY,GZ,AGI,NORM_X,NORM_Y,NORM_Z)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)

         CALL VELASSIGN(I,J,K,U_INTP,V_INTP,W_INTP,UT,VT,WT)

        IF (N .LT. NN2-1) THEN
         UF(I,J,K)=U_INTP
         VF(I,J,K)=V_INTP
         WF(I,J,K)=W_INTP
        ELSE
         UF(I,J,K)=U_INTP*0.5
         VF(I,J,K)=V_INTP*0.5
         WF(I,J,K)=W_INTP*0.5
        ENDIF

        ENDDO
       ENDDO

       N=N_MAX
!$OMP PARALLEL DO private(I,J,K)
       DO NN=1,NUMA(N,LLVS)
        I=I_B(NN,N,LLVS)
        J=J_B(NN,N,LLVS)
        K=K_B(NN,N,LLVS)
         UF(I,J,K)=0.
         VF(I,J,K)=0.
         WF(I,J,K)=0.
       ENDDO

!C===========================VELOCITY_INTP_END==========================C

!C=====SAVE
!       OPEN(146,FILE='0intp_vel.DAT')
!       WRITE(146,*) 'VARIABLES="X","Y","Z","UF","VF","WF"'
!      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!      DO K=1,N3FM
!      DO J=1,N2FM
!      DO I=1,N1FM
!       WRITE(146,147) XPF(I),YPF(J),ZPF(K),UF(i,j,k),VF(i,j,k),WF(i,j,k)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(146)
! 147  FORMAT(6F15.8)

       RETURN
       END

!C*******************************************************************
      SUBROUTINE VELASSIGN(I,J,K,U_INTP,V_INTP,W_INTP,UT,VT,WT)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_GEOM_VAR
      USE LVS_COUPLING
      
      IMPLICIT NONE

      REAL UT(-1:M1+2,-1:M2+2,-1:M3+2),VT(-1:M1+2,-1:M2+2,-1:M3+2),
     &                                    WT(-1:M1+2,-1:M2+2,-1:M3+2)
     
      INTEGER I,J,K
      REAL U_INTP,V_INTP,W_INTP

      INTEGER II,JJ,KK
      INTEGER IM1,IP1,IP2,JM1,JP1,JP2,KM1,KP1,KP2
      REAL XX,X1,X2,X3,X4,YY,Y1,Y2,Y3,Y4,ZZ,Z1,Z2,Z3,Z4
      REAL CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,CZ1,CZ2,CZ3,CZ4
      REAL AZ1,AZ2,AZ3,AZ4

      !THIS INTERPOLATION GUARANTEES C^2 CONTINUITY
      !INTERPOLATED VEL. IS DEFINED AT CELL CENTER?

!C***** UF
         II=ICOU_VEL(I)
         JJ=JCOUMP_VEL(J)
         KK=KCOUMP_VEL(K)
      !IB WALL TREATMENT
!       IF( FUNCBODY(X(II),YP(JJ),ZP(KK)) .LT. 0.) THEN
!            U_INTP=0.
!      ELSE

         XX=XF(I)
         X2=XP(II)
       IF (IPX .EQ. 1) THEN
        IF ( II .EQ. 1 ) THEN
         X1=X2-VDX(N1M)
        ELSE
         X1=X2-VDX(II)
        ENDIF
         X3=X2+VDX(IPV(II))
         X4=X3+VDX(IPV(IPV(II)))

       ELSE
         X1=X2-VDX(II)
       IF (II.EQ.N1M-1) THEN   !MIRROR
         X3=X2+VDX(IPV(II))
         X4=X3+2.*VDX(N1)  
       ELSE IF (II.EQ.N1M) THEN
         X3=X2+2.*VDX(N1)  
         X4=X3+VDX(N1M)
       ELSE
         X3=X2+VDX(IPV(II))
         X4=X3+VDX(IPV(IPV(II)))
       ENDIF
       ENDIF
        CX1=((XX-X2)*(XX-X3)*(XX-X4))/((X1-X2)*(X1-X3)*(X1-X4))
        CX2=((XX-X1)*(XX-X3)*(XX-X4))/((X2-X1)*(X2-X3)*(X2-X4))
        CX3=((XX-X1)*(XX-X2)*(XX-X4))/((X3-X1)*(X3-X2)*(X3-X4))
        CX4=((XX-X1)*(XX-X2)*(XX-X3))/((X4-X1)*(X4-X2)*(X4-X3))

         YY=YPF(J)
         Y2=YP(JJ)
       IF (IPY .EQ. 1) THEN
        IF ( JJ .EQ. 1 ) THEN
         Y1=Y2-VDY(N2M)
        ELSE
         Y1=Y2-VDY(JJ)
        ENDIF
         Y3=Y2+VDY(JPV(JJ))
         Y4=Y3+VDY(JPV(JPV(JJ)))

        ELSE           
         Y1=Y2-VDY(JJ)
       IF (JJ.EQ.N2M-1) THEN   !MIRROR
         Y3=Y2+VDY(JPV(JJ))
         Y4=Y3+2.*VDY(N2)  
       ELSE IF (JJ.EQ.N2M) THEN
         Y3=Y2+2.*VDY(N2)  
         Y4=Y3+VDY(N2M)
       ELSE
         Y3=Y2+VDY(JPV(JJ))
         Y4=Y3+VDY(JPV(JPV(JJ)))
       ENDIF
       ENDIF
        CY1=((YY-Y2)*(YY-Y3)*(YY-Y4))/((Y1-Y2)*(Y1-Y3)*(Y1-Y4))
        CY2=((YY-Y1)*(YY-Y3)*(YY-Y4))/((Y2-Y1)*(Y2-Y3)*(Y2-Y4))
        CY3=((YY-Y1)*(YY-Y2)*(YY-Y4))/((Y3-Y1)*(Y3-Y2)*(Y3-Y4))
        CY4=((YY-Y1)*(YY-Y2)*(YY-Y3))/((Y4-Y1)*(Y4-Y2)*(Y4-Y3))

      IF (N3M .EQ. 1) THEN
       CALL VEL_R1(UT,U_INTP,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KK)
      ELSE
         KM1=KMV(KK)
         KP1=KPV(KK)
         KP2=KPV(KP1)

         ZZ=ZPF(K)
         Z2=ZP(KK)
       IF (IPZ .EQ. 1) THEN
        IF ( KK .EQ. 1 ) THEN
         Z1=Z2-VDZ(N3M)
        ELSE
         Z1=Z2-VDZ(KK)
        ENDIF
         Z3=Z2+VDZ(KP1)
         Z4=Z3+VDZ(KP2)

       ELSE
         Z1=Z2-VDZ(KK)
       IF (KK.EQ.N3M-1) THEN   !MIRROR
         Z3=Z2+VDZ(KP1)
         Z4=Z3+2.*VDZ(N3)  
       ELSE IF (KK.EQ.N3M) THEN
         Z3=Z2+2.*VDZ(N3)  
         Z4=Z3+VDZ(N3M)
       ELSE
         Z3=Z2+VDZ(KP1)
         Z4=Z3+VDZ(KP2)
       ENDIF
       ENDIF
        CZ1=((ZZ-Z2)*(ZZ-Z3)*(ZZ-Z4))/((Z1-Z2)*(Z1-Z3)*(Z1-Z4))
        CZ2=((ZZ-Z1)*(ZZ-Z3)*(ZZ-Z4))/((Z2-Z1)*(Z2-Z3)*(Z2-Z4))
        CZ3=((ZZ-Z1)*(ZZ-Z2)*(ZZ-Z4))/((Z3-Z1)*(Z3-Z2)*(Z3-Z4))
        CZ4=((ZZ-Z1)*(ZZ-Z2)*(ZZ-Z3))/((Z4-Z1)*(Z4-Z2)*(Z4-Z3))

       CALL VEL_R1(UT,AZ1,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KM1)
       CALL VEL_R1(UT,AZ2,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KK)
       CALL VEL_R1(UT,AZ3,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KP1)
       CALL VEL_R1(UT,AZ4,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KP2)
        U_INTP=CZ1*AZ1+CZ2*AZ2+CZ3*AZ3+CZ4*AZ4
      ENDIF !      IF (N3M .EQ. 1) THEN
!      ENDIF

C***** VF
!      IF (FUNCBODY(XP(II),Y(JJ),ZP(KK)) .LT. 0. ) THEN
!            V_INTP=0.
!      ELSE
         II=ICOUMP_VEL(I)
         JJ=JCOU_VEL(J)
         KK=KCOUMP_VEL(K)

         XX=XPF(I)
         X2=XP(II)
       IF (IPX .EQ. 1) THEN
        IF ( II .EQ. 1 ) THEN
         X1=X2-VDX(N1M)
        ELSE
         X1=X2-VDX(II)
        ENDIF
         X3=X2+VDX(IPV(II))
         X4=X3+VDX(IPV(IPV(II)))

       ELSE
         X1=X2-VDX(II)
       IF (II.EQ.N1M-1) THEN   !MIRROR
         X3=X2+VDX(IPV(II))
         X4=X3+2.*VDX(N1)  
       ELSE IF (II.EQ.N1M) THEN
         X3=X2+2.*VDX(N1)  
         X4=X3+VDX(N1M)
       ELSE
         X3=X2+VDX(IPV(II))
         X4=X3+VDX(IPV(IPV(II)))
       ENDIF
       ENDIF
        CX1=((XX-X2)*(XX-X3)*(XX-X4))/((X1-X2)*(X1-X3)*(X1-X4))
        CX2=((XX-X1)*(XX-X3)*(XX-X4))/((X2-X1)*(X2-X3)*(X2-X4))
        CX3=((XX-X1)*(XX-X2)*(XX-X4))/((X3-X1)*(X3-X2)*(X3-X4))
        CX4=((XX-X1)*(XX-X2)*(XX-X3))/((X4-X1)*(X4-X2)*(X4-X3))

         YY=YF(J)
         Y2=YP(JJ)
       IF (IPY .EQ. 1) THEN
        IF ( JJ .EQ. 1 ) THEN
         Y1=Y2-VDY(N2M)
        ELSE
         Y1=Y2-VDY(JJ)
        ENDIF
         Y3=Y2+VDY(JPV(JJ))
         Y4=Y3+VDY(JPV(JPV(JJ)))

        ELSE           
         Y1=Y2-VDY(JJ)
       IF (JJ.EQ.N2M-1) THEN   !MIRROR
         Y3=Y2+VDY(JPV(JJ))
         Y4=Y3+2.*VDY(N2)  
       ELSE IF (JJ.EQ.N2M) THEN
         Y3=Y2+2.*VDY(N2)  
         Y4=Y3+VDY(N2M)
       ELSE
         Y3=Y2+VDY(JPV(JJ))
         Y4=Y3+VDY(JPV(JPV(JJ)))
       ENDIF
       ENDIF
        CY1=((YY-Y2)*(YY-Y3)*(YY-Y4))/((Y1-Y2)*(Y1-Y3)*(Y1-Y4))
        CY2=((YY-Y1)*(YY-Y3)*(YY-Y4))/((Y2-Y1)*(Y2-Y3)*(Y2-Y4))
        CY3=((YY-Y1)*(YY-Y2)*(YY-Y4))/((Y3-Y1)*(Y3-Y2)*(Y3-Y4))
        CY4=((YY-Y1)*(YY-Y2)*(YY-Y3))/((Y4-Y1)*(Y4-Y2)*(Y4-Y3))

      IF (N3M .EQ. 1) THEN
       CALL VEL_R2(VT,V_INTP,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KK)
      ELSE
         KM1=KMV(KK)
         KP1=KPV(KK)
         KP2=KPV(KP1)
         
         ZZ=ZPF(K)
         Z2=ZP(KK)
       IF (IPZ .EQ. 1) THEN
        IF ( KK .EQ. 1 ) THEN
         Z1=Z2-VDZ(N3M)
        ELSE
         Z1=Z2-VDZ(KK)
        ENDIF
         Z3=Z2+VDZ(KP1)
         Z4=Z3+VDZ(KP2)

       ELSE
         Z1=Z2-VDZ(KK)
       IF (KK.EQ.N3M-1) THEN   !MIRROR
         Z3=Z2+VDZ(KP1)
         Z4=Z3+2.*VDZ(N3)  
       ELSE IF (KK.EQ.N3M) THEN
         Z3=Z2+2.*VDZ(N3)  
         Z4=Z3+VDZ(N3M)
       ELSE
         Z3=Z2+VDZ(KP1)
         Z4=Z3+VDZ(KP2)
       ENDIF
       ENDIF
        CZ1=((ZZ-Z2)*(ZZ-Z3)*(ZZ-Z4))/((Z1-Z2)*(Z1-Z3)*(Z1-Z4))
        CZ2=((ZZ-Z1)*(ZZ-Z3)*(ZZ-Z4))/((Z2-Z1)*(Z2-Z3)*(Z2-Z4))
        CZ3=((ZZ-Z1)*(ZZ-Z2)*(ZZ-Z4))/((Z3-Z1)*(Z3-Z2)*(Z3-Z4))
        CZ4=((ZZ-Z1)*(ZZ-Z2)*(ZZ-Z3))/((Z4-Z1)*(Z4-Z2)*(Z4-Z3))

       CALL VEL_R2(VT,AZ1,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KM1)
       CALL VEL_R2(VT,AZ2,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KK)
       CALL VEL_R2(VT,AZ3,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KP1)
       CALL VEL_R2(VT,AZ4,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KP2)
          V_INTP=CZ1*AZ1+CZ2*AZ2+CZ3*AZ3+CZ4*AZ4
      ENDIF !      IF (N3M .EQ. 1) THEN
!      ENDIF

C***** WF 
      IF (N3M .EQ. 1) THEN
        W_INTP=0.
      ELSE
!      IF (FUNCBODY(XP(II),YP(JJ),Z(KK)) .LT. 0. ) THEN
!            W_INTP=0.
!      ELSE
         II=ICOUMP_VEL(I)
         JJ=JCOUMP_VEL(J)
         KK=KCOU_VEL(K)

         XX=XPF(I)
         X2=XP(II)
       IF (IPX .EQ. 1) THEN
        IF ( II .EQ. 1 ) THEN
         X1=X2-VDX(N1M)
        ELSE
         X1=X2-VDX(II)
        ENDIF
         X3=X2+VDX(IPV(II))
         X4=X3+VDX(IPV(IPV(II)))

       ELSE
         X1=X2-VDX(II)
       IF (II.EQ.N1M-1) THEN   !MIRROR
         X3=X2+VDX(IPV(II))
         X4=X3+2.*VDX(N1)  
       ELSE IF (II.EQ.N1M) THEN
         X3=X2+2.*VDX(N1)  
         X4=X3+VDX(N1M)
       ELSE
         X3=X2+VDX(IPV(II))
         X4=X3+VDX(IPV(IPV(II)))
       ENDIF
       ENDIF
        CX1=((XX-X2)*(XX-X3)*(XX-X4))/((X1-X2)*(X1-X3)*(X1-X4))
        CX2=((XX-X1)*(XX-X3)*(XX-X4))/((X2-X1)*(X2-X3)*(X2-X4))
        CX3=((XX-X1)*(XX-X2)*(XX-X4))/((X3-X1)*(X3-X2)*(X3-X4))
        CX4=((XX-X1)*(XX-X2)*(XX-X3))/((X4-X1)*(X4-X2)*(X4-X3))

         YY=YPF(J)
         Y2=YP(JJ)
       IF (IPY .EQ. 1) THEN
        IF ( JJ .EQ. 1 ) THEN
         Y1=Y2-VDY(N2M)
        ELSE
         Y1=Y2-VDY(JJ)
        ENDIF
         Y3=Y2+VDY(JPV(JJ))
         Y4=Y3+VDY(JPV(JPV(JJ)))

        ELSE           
         Y1=Y2-VDY(JJ)
       IF (JJ.EQ.N2M-1) THEN   !MIRROR
         Y3=Y2+VDY(JPV(JJ))
         Y4=Y3+2.*VDY(N2)  
       ELSE IF (JJ.EQ.N2M) THEN
         Y3=Y2+2.*VDY(N2)  
         Y4=Y3+VDY(N2M)
       ELSE
         Y3=Y2+VDY(JPV(JJ))
         Y4=Y3+VDY(JPV(JPV(JJ)))
       ENDIF
       ENDIF
        CY1=((YY-Y2)*(YY-Y3)*(YY-Y4))/((Y1-Y2)*(Y1-Y3)*(Y1-Y4))
        CY2=((YY-Y1)*(YY-Y3)*(YY-Y4))/((Y2-Y1)*(Y2-Y3)*(Y2-Y4))
        CY3=((YY-Y1)*(YY-Y2)*(YY-Y4))/((Y3-Y1)*(Y3-Y2)*(Y3-Y4))
        CY4=((YY-Y1)*(YY-Y2)*(YY-Y3))/((Y4-Y1)*(Y4-Y2)*(Y4-Y3))

         KM1=KMV(KK)
         KP1=KPV(KK)
         KP2=KPV(KP1)

         ZZ=ZF(K)
         Z2=ZP(KK)
       IF (IPZ .EQ. 1) THEN
        IF ( KK .EQ. 1 ) THEN
         Z1=Z2-VDZ(N3M)
        ELSE
         Z1=Z2-VDZ(KK)
        ENDIF
         Z3=Z2+VDZ(KP1)
         Z4=Z3+VDZ(KP2)

       ELSE
         Z1=Z2-VDZ(KK)
       IF (KK.EQ.N3M-1) THEN   !MIRROR
         Z3=Z2+VDZ(KP1)
         Z4=Z3+2.*VDZ(N3)  
       ELSE IF (KK.EQ.N3M) THEN
         Z3=Z2+2.*VDZ(N3)  
         Z4=Z3+VDZ(N3M)
       ELSE
         Z3=Z2+VDZ(KP1)
         Z4=Z3+VDZ(KP2)
       ENDIF
       ENDIF

        CZ1=((ZZ-Z2)*(ZZ-Z3)*(ZZ-Z4))/((Z1-Z2)*(Z1-Z3)*(Z1-Z4))
        CZ2=((ZZ-Z1)*(ZZ-Z3)*(ZZ-Z4))/((Z2-Z1)*(Z2-Z3)*(Z2-Z4))
        CZ3=((ZZ-Z1)*(ZZ-Z2)*(ZZ-Z4))/((Z3-Z1)*(Z3-Z2)*(Z3-Z4))
        CZ4=((ZZ-Z1)*(ZZ-Z2)*(ZZ-Z3))/((Z4-Z1)*(Z4-Z2)*(Z4-Z3))

       CALL VEL_R3(WT,AZ1,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KM1)
       CALL VEL_R3(WT,AZ2,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KK)
       CALL VEL_R3(WT,AZ3,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KP1)
       CALL VEL_R3(WT,AZ4,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KP2)
          W_INTP=CZ1*AZ1+CZ2*AZ2+CZ3*AZ3+CZ4*AZ4
!      ENDIF
      ENDIF!IF (N3M .EQ. 1)

        RETURN
        END




!C*******************************************************************      
      SUBROUTINE VEL_R1(UT,AZZ,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KK)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL UT(-1:M1+2,-1:M2+2,-1:M3+2)
      
      INTEGER II,JJ,KK
      REAL AZZ
      
      INTEGER IM1,IP1,IP2,IP3,JM1,JP1,JP2
      REAL CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4
      REAL AY1,AY2,AY3,AY4

        IM1=IMV(II)
        IP1=IPV(II)
        IP2=IPV(IP1)
        IP3=IPV(IP2)

        JM1=JMV(JJ)
        JP1=JPV(JJ)
        JP2=JPV(JP1)

        AY1=CX1*0.5*( UT(IM1,JM1,KK)+UT(II,JM1,KK) )
     &     +CX2*0.5*( UT(II,JM1,KK)+UT(IP1,JM1,KK) )
     &     +CX3*0.5*( UT(IP1,JM1,KK)+UT(IP2,JM1,KK) )
     &     +CX4*0.5*( UT(IP2,JM1,KK)+UT(IP3,JM1,KK) )
        AY2=CX1*0.5*( UT(IM1,JJ,KK)+UT(II,JJ,KK) )
     &     +CX2*0.5*( UT(II,JJ,KK)+UT(IP1,JJ,KK) )
     &     +CX3*0.5*( UT(IP1,JJ,KK)+UT(IP2,JJ,KK) )
     &     +CX4*0.5*( UT(IP2,JJ,KK)+UT(IP3,JJ,KK) )
        AY3=CX1*0.5*( UT(IM1,JP1,KK)+UT(II,JP1,KK) )
     &     +CX2*0.5*( UT(II,JP1,KK)+UT(IP1,JP1,KK) )
     &     +CX3*0.5*( UT(IP1,JP1,KK)+UT(IP2,JP1,KK) )
     &     +CX4*0.5*( UT(IP2,JP1,KK)+UT(IP3,JP1,KK) )
        AY4=CX1*0.5*( UT(IM1,JP2,KK)+UT(II,JP2,KK) )
     &     +CX2*0.5*( UT(II,JP2,KK)+UT(IP1,JP2,KK) )
     &     +CX3*0.5*( UT(IP1,JP2,KK)+UT(IP2,JP2,KK) )
     &     +CX4*0.5*( UT(IP2,JP2,KK)+UT(IP3,JP2,KK) )

        AZZ=CY1*AY1+CY2*AY2+CY3*AY3+CY4*AY4

      RETURN
      END 

!C*******************************************************************      
      SUBROUTINE VEL_R2(VT,AZZ,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KK)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE
      
      REAL VT(-1:M1+2,-1:M2+2,-1:M3+2)
      
      INTEGER II,JJ,KK
      REAL AZZ

      INTEGER JM1,JP1,JP2,JP3,IM1,IP1,IP2
      REAL CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4
      REAL AY1,AY2,AY3,AY4
      
        JM1=JMV(JJ)
        JP1=JPV(JJ)
        JP2=JPV(JP1)
        JP3=JPV(JP2)

        IM1=IMV(II)
        IP1=IPV(II)
        IP2=IPV(IP1)

        AY1=CX1*0.5*( VT(IM1,JM1,KK)+VT(IM1,JJ,KK) )
     &     +CX2*0.5*( VT(II,JM1,KK)+VT(II,JJ,KK) )
     &     +CX3*0.5*( VT(IP1,JM1,KK)+VT(IP1,JJ,KK) )
     &     +CX4*0.5*( VT(IP2,JM1,KK)+VT(IP2,JJ,KK) )
        AY2=CX1*0.5*( VT(IM1,JJ,KK)+VT(IM1,JP1,KK) )
     &     +CX2*0.5*( VT(II,JJ,KK)+VT(II,JP1,KK) )
     &     +CX3*0.5*( VT(IP1,JJ,KK)+VT(IP1,JP1,KK) )
     &     +CX4*0.5*( VT(IP2,JJ,KK)+VT(IP2,JP1,KK) )
        AY3=CX1*0.5*( VT(IM1,JP1,KK)+VT(IM1,JP2,KK) )
     &     +CX2*0.5*( VT(II,JP1,KK)+VT(II,JP2,KK) )
     &     +CX3*0.5*( VT(IP1,JP1,KK)+VT(IP1,JP2,KK) )
     &     +CX4*0.5*( VT(IP2,JP1,KK)+VT(IP2,JP2,KK) )
        AY4=CX1*0.5*( VT(IM1,JP2,KK)+VT(IM1,JP3,KK) )
     &     +CX2*0.5*( VT(II,JP2,KK)+VT(II,JP3,KK) )
     &     +CX3*0.5*( VT(IP1,JP2,KK)+VT(IP1,JP3,KK) )
     &     +CX4*0.5*( VT(IP2,JP2,KK)+VT(IP2,JP3,KK) )

        AZZ=CY1*AY1+CY2*AY2+CY3*AY3+CY4*AY4

        RETURN
        END

!C*******************************************************************
      SUBROUTINE VEL_R3(WT,AZZ,CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,II,JJ,KK)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE      
      
      REAL WT(-1:M1+2,-1:M2+2,-1:M3+2)

      INTEGER II,JJ,KK
      REAL AZZ

      INTEGER JM1,JP1,JP2,KP1,IM1,IP1,IP2
      REAL CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4
      REAL AY1,AY2,AY3,AY4
      
        JM1=JMV(JJ)
        JP1=JPV(JJ)
        JP2=JPV(JP1)

        KP1=KPV(KK)

        IM1=IMV(II)
        IP1=IPV(II)
        IP2=IPV(IP1)

       AY1=CX1*0.5*( WT(IM1,JM1,KK)+WT(IM1,JM1,KP1) )
     &    +CX2*0.5*( WT(II,JM1,KK)+WT(II,JM1,KP1) )
     &    +CX3*0.5*( WT(IP1,JM1,KK)+WT(IP1,JM1,KP1) )
     &    +CX4*0.5*( WT(IP2,JM1,KK)+WT(IP2,JM1,KP1) )
       AY2=CX1*0.5*( WT(IM1,JJ,KK)+WT(IM1,JJ,KP1) )
     &    +CX2*0.5*( WT(II,JJ,KK)+WT(II,JJ,KP1) )
     &    +CX3*0.5*( WT(IP1,JJ,KK)+WT(IP1,JJ,KP1) )
     &    +CX4*0.5*( WT(IP2,JJ,KK)+WT(IP2,JJ,KP1) )
       AY3=CX1*0.5*( WT(IM1,JP1,KK)+WT(IM1,JP1,KP1) )
     &    +CX2*0.5*( WT(II,JP1,KK)+WT(II,JP1,KP1) )
     &    +CX3*0.5*( WT(IP1,JP1,KK)+WT(IP1,JP1,KP1) )
     &    +CX4*0.5*( WT(IP2,JP1,KK)+WT(IP2,JP1,KP1) )
       AY4=CX1*0.5*( WT(IM1,JP2,KK)+WT(IM1,JP2,KP1) )
     &    +CX2*0.5*( WT(II,JP2,KK)+WT(II,JP2,KP1) )
     &    +CX3*0.5*( WT(IP1,JP2,KK)+WT(IP1,JP2,KP1) )
     &    +CX4*0.5*( WT(IP2,JP2,KK)+WT(IP2,JP2,KP1) )

        AZZ=CY1*AY1+CY2*AY2+CY3*AY3+CY4*AY4

        RETURN
        END


 
!!C*******************************************************************
      SUBROUTINE GLOBAL_MASS_CORRECTION(NN1,LLVS,VOL_TOT_ORI,dtconst)
!!C*******************************************************************
      !!ZHANG ET AL.(2010)
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      use two_phase_property
      USE PHASE_BOILING
      
      IMPLICIT NONE

      REAL*8 PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)
      REAL*8 VOL_TOT_ORI(MLVS)
      REAL*8 GPSI
      INTEGER NN1,LLVS

      INTEGER*8 I,J,K,N,NN,L
      REAL*8 RECOR,DVOL,VOL_TOT,VOL_BAND,VOL_IN,COR
      real*8 dtconst

      IF (N3FM .EQ. 1 ) THEN !CFL=1.
        RECOR=MIN(SDXF,SDYF)
        DVOL=SDXF*SDYF
      ELSE
        RECOR=MIN(SDXF,SDYF,SDZF)
        DVOL=SDXF*SDYF*SDZF
      ENDIF
        RECOR=RECOR*1.5 !ACCELERATOR
       
       
      ! ADJUST THE VOL_TOT_ORI    
      ! INCREASED AMOUNT OF VOLUME RHO_G/RHO_F*MFDOT*DT
      !CHECK

      IF (IPHS .EQ. 1) THEN 
      WRITE(*,*) 'OLD VOL=',VOL_TOT_ORI(LLVS)
      VOL_TOT_ORI(LLVS)=VOL_TOT_IPHS(LLVS)      
      !CHECK
      WRITE(*,*) 'NEW VOL=',VOL_TOT_ORI(LLVS)
      ENDIF  

      
       DO L=1,15

      !CALCULATE_ORIGINAL_TOTAL_MASS OF BUBBLE
       IF (L.EQ.1) THEN
        CALL CAL_PSIF_MTP(1,NN1,LLVS,PSIF_MTP)
        VOL_TOT=0D0
!$OMP PARALLEL DO private(I,J)
!$OMP&reduction(+:VOL_TOT)
        DO K=1,N3FM
        DO J=1,N2FM
        DO I=1,N1FM
         VOL_TOT=VOL_TOT+(1d0-PSIF_MTP(I,J,K))*DVOL
        ENDDO
        ENDDO
        ENDDO      
      write(*,*) 'vol_tot=',vol_tot
      
                       !stop

        VOL_BAND=0.
        DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
!$OMP&reduction(+:VOL_BAND)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
          VOL_BAND=VOL_BAND+(1.-PSIF_MTP(I,J,K))*DVOL
        ENDDO
        ENDDO

        VOL_IN=VOL_TOT-VOL_BAND
       ELSE
        CALL CAL_PSIF_MTP(0,NN1,LLVS,PSIF_MTP)
        VOL_BAND=0.
        DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
!$OMP&reduction(+:VOL_BAND)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
          VOL_BAND=VOL_BAND+(1.-PSIF_MTP(I,J,K))*DVOL
        ENDDO
        ENDDO

        VOL_TOT=VOL_IN+VOL_BAND
       ENDIF

       !VOLUME ERROR EXIST
       IF (ABS(VOL_TOT-VOL_TOT_ORI(LLVS)) .GE. 1.D-5) THEN
         !PSI=0 IS REFERENCE AREA. SO MINUS IS NEEDED. REFER TO ZHANG ET AL.(2010)
         COR=RECOR*(-(VOL_TOT_ORI(LLVS)-VOL_TOT)/VOL_TOT_ORI(LLVS) )

       DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         ALPHI_EXT(I,J,K,1)=ALPHI_EXT(I,J,K,1)+COR
        ENDDO
       ENDDO

       ELSE
        GOTO 100
       ENDIF

       ENDDO !DO L=1,20

 100   CONTINUE
      IF ( ABS(VOL_TOT-VOL_TOT_ORI(LLVS)) .GE. 1.D-5 ) THEN
      WRITE(*,*) 'GLOBAL_MASS_CONSERVATION IS NOT CONSERVED!!!',llvs
     &,ABS(VOL_TOT-VOL_TOT_ORI(LLVS))
      ENDIF

        RETURN
        END

!!C*******************************************************************
      SUBROUTINE BUBBLE_TRACKING(LLVS,PSI_TMP,U,V,W)
!!C*******************************************************************

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR
      USE LVS_VAR
      
      IMPLICIT NONE

      REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL PSI_TMP(0:M1,0:M2,0:M3)

      INTEGER LLVS
      
      INTEGER ICACY,IZ_PERIODIC
      INTEGER I,J,K
      REAL PI,VOL,VOLF,VOLI
      REAL XX,YY,ZZ,RR,TT,TT_TMP
      REAL UU,VV,WW,UR,UT
      REAL X_CEN,Y_CEN,Z_CEN
      REAL AR_BUB,ANG_BUB1,ANG_BUB2

       ICACY=1 !1:U,V, 2:UR,UT

      IF (MSUB .EQ. 3) THEN
       PI=ACOS(-1.)

        IZ_PERIODIC=0  !THIS IS FOR THE PERIODIC B.C
!FIND THE POSITION OF CM(CENTER OF MASS)
!$OMP PARALLEL DO private(I)
        DO J=1,N2M
        DO I=1,N1M
         IF (PSI_TMP(I,J,N3M) .NE. 1.) THEN
           IZ_PERIODIC=1
!           GOTO 10 !cannot used in PARALLEL COMPUTING.
         ENDIF
        ENDDO
        ENDDO
! 10    CONTINUE

      IF (ICACY .EQ. 1) THEN
       XX=0.
       YY=0.
       ZZ=0.
       UU=0.
       VV=0.
       WW=0.

       VOL=0.

       IF (IPZ .EQ. 1 .AND. IZ_PERIODIC .EQ. 1) THEN
!FIND THE POSITION OF CM(CENTER OF MASS)
!$OMP PARALLEL DO private(I,J,VOLF)
!$OMP&reduction(+:XX,YY,ZZ,UU,VV,WW,VOL)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         IF (PSI_TMP(I,J,K) .NE. 1.) THEN
           VOLF=(1.-PSI_TMP(I,J,K))*SDX(I)*SDY(J)*SDZ(K)

          XX=XX+XP(I)*VOLF
          YY=YY+YP(J)*VOLF
          IF (ZP(K) .LE. 0.5*ZL) THEN
           ZZ=ZZ+(ZP(K)+ZL)*VOLF
          ELSE
           ZZ=ZZ+ZP(K)*VOLF
          ENDIF

          UU=UU+0.5*(U(I,J,K)+U(IPV(I),J,K))*VOLF
          VV=VV+0.5*(V(I,J,K)+V(I,JPV(J),K))*VOLF
          WW=WW+0.5*(W(I,J,K)+W(I,J,KPV(K)))*VOLF
          VOL=VOL+VOLF
         ENDIF
        ENDDO
        ENDDO
        ENDDO

        ELSE      
!FIND THE POSITION OF CM(CENTER OF MASS)
!$OMP PARALLEL DO private(I,J,VOLF)
!$OMP&reduction(+:XX,YY,ZZ,UU,VV,WW,VOL)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         IF (PSI_TMP(I,J,K) .NE. 1.) THEN
           VOLF=(1.-PSI_TMP(I,J,K))*SDX(I)*SDY(J)*SDZ(K)

          XX=XX+XP(I)*VOLF
          YY=YY+YP(J)*VOLF
          ZZ=ZZ+ZP(K)*VOLF

          UU=UU+0.5*(U(I,J,K)+U(IPV(I),J,K))*VOLF
          VV=VV+0.5*(V(I,J,K)+V(I,JPV(J),K))*VOLF
          WW=WW+0.5*(W(I,J,K)+W(I,J,KPV(K)))*VOLF
          VOL=VOL+VOLF
         ENDIF
        ENDDO
        ENDDO
        ENDDO
       ENDIF
         IF (VOL .NE. 0.) THEN
           VOLI=1./VOL
           XX=XX*VOLI
           YY=YY*VOLI
           ZZ=ZZ*VOLI
           UU=UU*VOLI
           VV=VV*VOLI
           WW=WW*VOLI
         ENDIF

       X_CEN=XX
       Y_CEN=YY
       Z_CEN=ZZ

       CALL BUBBLE_EIGENV(X_CEN,Y_CEN,Z_CEN,IZ_PERIODIC,VOL
     &,PSI_TMP,AR_BUB,ANG_BUB1,ANG_BUB2)  

      OPEN(100,FILE='0bub_tracking.dat',POSITION='APPEND')
      WRITE(100,101)TIME,LLVS,XX,YY,ZZ,UU,VV,WW,AR_BUB,ANG_BUB1,ANG_BUB2
      CLOSE(100)
 101   FORMAT(1ES15.7,1I5,9ES15.7)

      ELSE
       XX=0.
       YY=0.
       ZZ=0.
       UU=0.
       VV=0.
       WW=0.

       VOL=0.

       IF (IPZ .EQ. 1 .AND. IZ_PERIODIC .EQ. 1) THEN

!FIND THE POSITION OF CM(CENTER OF MASS)
!$OMP PARALLEL DO private(I,J,VOLF,TT_TMP)
!$OMP&reduction(+:XX,YY,ZZ,UU,VV,WW,VOL)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         IF (PSI_TMP(I,J,K) .NE. 1.) THEN
          VOLF=(1.-PSI_TMP(I,J,K))*SDX(I)*SDY(J)*SDZ(K)

          XX=XX+XP(I)*VOLF
          YY=YY+YP(J)*VOLF
          IF (ZP(K) .LE. 0.5*ZL) THEN
           ZZ=ZZ+(ZP(K)+ZL)*VOLF
          ELSE
           ZZ=ZZ+ZP(K)*VOLF
          ENDIF

          UU=UU+0.5*(U(I,J,K)+U(IPV(I),J,K))*VOLF
          VV=VV+0.5*(V(I,J,K)+V(I,JPV(J),K))*VOLF
          WW=WW+0.5*(W(I,J,K)+W(I,J,KPV(K)))*VOLF
          VOL=VOL+VOLF
         ENDIF
        ENDDO
        ENDDO
        ENDDO

        ELSE  
!FIND THE POSITION OF CM(CENTER OF MASS)
!$OMP PARALLEL DO private(I,J,VOLF,TT_TMP)
!$OMP&reduction(+:XX,YY,ZZ,UU,VV,WW,VOL)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         IF (PSI_TMP(I,J,K) .NE. 1.) THEN
         VOLF=(1.-PSI_TMP(I,J,K))*SDX(I)*SDY(J)*SDZ(K)

          XX=XX+XP(I)*VOLF
          YY=YY+YP(J)*VOLF
          ZZ=ZZ+ZP(K)*VOLF

          UU=UU+0.5*(U(I,J,K)+U(IPV(I),J,K))*VOLF
          VV=VV+0.5*(V(I,J,K)+V(I,JPV(J),K))*VOLF
          WW=WW+0.5*(W(I,J,K)+W(I,J,KPV(K)))*VOLF
          VOL=VOL+VOLF
         ENDIF
        ENDDO
        ENDDO
        ENDDO
       ENDIF

         IF (VOL .NE. 0.) THEN
           VOLI=1./VOL
           XX=XX*VOLI
           YY=YY*VOLI
           ZZ=ZZ*VOLI
           UU=UU*VOLI
           VV=VV*VOLI
           WW=WW*VOLI
         ENDIF

         RR=SQRT(XX**2+YY**2)
         TT_TMP=ATAN(YY/XX)
         IF (XX .LE. 0. .AND. YY .GE. 0.) THEN !2ND QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .LE. 0. .AND. YY .LE. 0.) THEN !3RD QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .GE. 0. .AND. YY .LE. 0.) THEN !4TH QUADRANT
           TT_TMP=TT_TMP+2.*PI
         ENDIF
         TT=TT_TMP

       UR=UU*COS(TT)+VV*SIN(TT)
       UT=-UU*SIN(TT)+VV*COS(TT)

       X_CEN=RR*COS(TT)
       Y_CEN=RR*SIN(TT)
       Z_CEN=ZZ

       
       CALL BUBBLE_EIGENV(X_CEN,Y_CEN,Z_CEN,IZ_PERIODIC,VOL
     &,PSI_TMP,AR_BUB,ANG_BUB1,ANG_BUB2)  

      OPEN(200,FILE='0bub_tracking.dat',POSITION='APPEND')
      WRITE(200,201)TIME,LLVS,RR,TT,ZZ,UR,UT,WW,AR_BUB,ANG_BUB1,ANG_BUB2
      CLOSE(200)
 201   FORMAT(1ES15.7,1I5,9ES15.7)

      ENDIF
      
              !2017-09-22_T_CON SAVE
        VELVS(LLVS,1)=UU
        VELVS(LLVS,2)=VV
        VELVS(LLVS,3)=WW

      ENDIF !IF (MSUB .EQ. 3) THEN
          
        
        RETURN
        END

!!C*******************************************************************
      SUBROUTINE BUBBLE_EIGENV(X_CEN,Y_CEN,Z_CEN,IZ_PERIODIC,VOL
     &,PSI_TMP,AR_BUB,ANG_BUB1,ANG_BUB2)    
!!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL PSI_TMP(0:M1,0:M2,0:M3)
      REAL AMI(3,3),EIGENV(3,3) !2ND MOMENT OF INERTIA
      
      REAL X_CEN,Y_CEN,Z_CEN,VOL,AR_BUB,ANG_BUB1,ANG_BUB2
      INTEGER IZ_PERIODIC
      
      REAL AMI_XX,AMI_YY,AMI_ZZ,AMI_XY,AMI_YZ,AMI_ZX
      INTEGER I,J,K,ISMALL
      REAL VOLF,DELX,DELY,DELZ,VOLI,DENOM

!c-----calculate 2nd moment of inertia.
       AMI_XX=0.
       AMI_YY=0.
       AMI_ZZ=0.
       AMI_XY=0.
       AMI_YZ=0.
       AMI_ZX=0.

       IF (IPZ .EQ. 1 .AND. IZ_PERIODIC .EQ. 1) THEN
!FIND THE POSITION OF CM(CENTER OF MASS)
!$OMP PARALLEL DO private(I,J,VOLF,DELX,DELY,DELZ)
!$OMP&reduction(+:AMI_XX,AMI_YY,AMI_ZZ,AMI_XY,AMI_YZ,AMI_ZX)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         IF (PSI_TMP(I,J,K) .NE. 1.) THEN
          VOLF=(1.-PSI_TMP(I,J,K))*SDX(I)*SDY(J)*SDZ(K)

         DELX=XP(I)-X_CEN
         DELY=YP(J)-Y_CEN

          IF (ZP(K) .LE. 0.5*ZL) THEN
         DELZ=(ZP(K)+ZL)-Z_CEN
          ELSE
         DELZ=ZP(K)-Z_CEN
          ENDIF

         AMI_XX=AMI_XX+DELX**2*VOLF
         AMI_YY=AMI_YY+DELY**2*VOLF
         AMI_ZZ=AMI_ZZ+DELZ**2*VOLF
         AMI_XY=AMI_XY+DELX*DELY*VOLF
         AMI_YZ=AMI_YZ+DELY*DELZ*VOLF
         AMI_ZX=AMI_ZX+DELZ*DELX*VOLF

        ENDIF !IF (PSI_TMP(I,J,K) .NE. 1.) THEN

        ENDDO
        ENDDO
        ENDDO

        ELSE  
!FIND THE POSITION OF CM(CENTER OF MASS)
!$OMP PARALLEL DO private(I,J,VOLF,DELX,DELY,DELZ)
!$OMP&reduction(+:AMI_XX,AMI_YY,AMI_ZZ,AMI_XY,AMI_YZ,AMI_ZX)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         IF (PSI_TMP(I,J,K) .NE. 1.) THEN
         VOLF=(1.-PSI_TMP(I,J,K))*SDX(I)*SDY(J)*SDZ(K)

         DELX=XP(I)-X_CEN
         DELY=YP(J)-Y_CEN
         DELZ=ZP(K)-Z_CEN

         AMI_XX=AMI_XX+DELX**2*VOLF
         AMI_YY=AMI_YY+DELY**2*VOLF
         AMI_ZZ=AMI_ZZ+DELZ**2*VOLF
         AMI_XY=AMI_XY+DELX*DELY*VOLF
         AMI_YZ=AMI_YZ+DELY*DELZ*VOLF
         AMI_ZX=AMI_ZX+DELZ*DELX*VOLF

         ENDIF !IF (PSI_TMP(I,J,K) .NE. 1.) THEN
        ENDDO
        ENDDO
        ENDDO
       ENDIF !IF (IPZ .EQ. 1 .AND. IZ_PERIODIC .EQ. 1) THEN
         IF (VOL .NE. 0.) THEN
           VOLI=1./VOL
           AMI(1,1)=AMI_XX*VOLI
           AMI(2,2)=AMI_YY*VOLI
           AMI(3,3)=AMI_ZZ*VOLI
           AMI(1,2)=AMI_XY*VOLI
           AMI(2,3)=AMI_YZ*VOLI
           AMI(3,1)=AMI_ZX*VOLI
         ENDIF
          AMI(2,1)=AMI(1,2)
          AMI(3,2)=AMI(2,3)
          AMI(1,3)=AMI(3,1)

        CALL jacobi_bub(AMI,EIGENV,1.E-12,3) !CALCULATE EIGENVALUE & EIGENVECTOR
      ! write(*,*) ami(1,1),ami(2,2),ami(3,3)
        !THE RATIO OF LARGEST AND SMALLEST EIGENVALUE
        IF (AMI(1,1) .GE. AMI(2,2)) THEN
          IF (AMI(1,1) .GE. AMI(3,3)) THEN
           IF (AMI(2,2) .GE. AMI(3,3)) THEN
             AR_BUB=SQRT(AMI(1,1)/AMI(3,3))
             ISMALL=3
           ELSE
             AR_BUB=SQRT(AMI(1,1)/AMI(2,2))
             ISMALL=2
           ENDIF
          ELSE
           AR_BUB=SQRT(AMI(3,3)/AMI(2,2))
             ISMALL=2
          ENDIF
        ELSE
          IF (AMI(2,2) .GE. AMI(3,3)) THEN
           IF (AMI(1,1) .GE. AMI(3,3)) THEN
            AR_BUB=SQRT(AMI(2,2)/AMI(3,3)  )
             ISMALL=3
           ELSE
             AR_BUB=SQRT(AMI(2,2)/AMI(1,1))
             ISMALL=1
           ENDIF  
          ELSE
           AR_BUB=SQRT(AMI(3,3)/AMI(1,1))
             ISMALL=1
          ENDIF
        ENDIF

        !ANGLE BETWEEN EIGENVECTOR OF SMALLEST EIGENVALUE AND A VERTICAL DIRECTION.
        DENOM=1./SQRT(EIGENV(1,ISMALL)**2
     &                     +EIGENV(2,ISMALL)**2+EIGENV(3,ISMALL)**2)
        ANG_BUB1=ACOS(EIGENV(1,ISMALL)*DENOM)*57.295779513082321  !180/PI
        ANG_BUB2=ACOS(EIGENV(3,ISMALL)*DENOM)*57.295779513082321  !180/PI

        IF (ANG_BUB1 .GE. 90.) ANG_BUB1=180.-ANG_BUB1
        IF (ANG_BUB2 .GE. 90.) ANG_BUB2=180.-ANG_BUB2

!       WRITE(*,*) AR_BUB,ANG_BUB1,ANG_BUB2
!       WRITE(*,*) EIGENV(1,ISMALL),EIGENV(2,ISMALL),EIGENV(3,ISMALL)

!        WRITE(*,*) AMI(1,1),AMI(2,2),AMI(3,3)
!        WRITE(*,*) EIGENV(1,1),EIGENV(1,2),EIGENV(1,3)
!        WRITE(*,*) EIGENV(2,1),EIGENV(2,2),EIGENV(2,3)
!        WRITE(*,*) EIGENV(3,1),EIGENV(3,2),EIGENV(3,3)

        RETURN
        END       

  
!*******************************************************************
      SUBROUTINE CAL_COALESCENCE(LVS2,LVS1,PSI_CN,VOL_TOT_ORI)
!*******************************************************************   
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE
      
      REAL ALPHI_C(M1FM,M2FM,M3FM)      
      REAL D_CRI,TEMP                                         
      REAL DX,DY,DZ,DVOL 
      REAL XX,YY,ZZ,ALPHI_EXT2(-2:M1F+3,-2:M2F+3,-2:M3F+3)     
      REAL VOL_TOT_ORI(MLVS)
      REAL PSI_TMP(0:M1,0:M2,0:M3)
      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)
      REAL PSI_CN(0:M1,0:M2,0:M3)
      REAL psif3(M1L:M1U,M2L:M2U,M3L:M3U)
      REAL ALPHI2
      REAL VOL_TEMP
      
      INTEGER MASK_LVS(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      INTEGER MASK_S(-2:M1F+3,-2:M2F+3,-2:M3F+3)   
      INTEGER NN1,NN2,LVS1,LVS2,NM,LLVS,L     
      INTEGER I,J,K,N,NN,II,JJ,KK,I2,J2,K2
      INTEGER NUMS,NUMS_INITIAL,NUMS_CRI
      
      nn1=n_max-3
      nn2=n_max-1
! INITIALIZE THE MASKS, AND WORKING ALPHI_EXT
!$OMP PARALLEL DO            
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       ALPHI_EXT(I,J,K,1)=0.
       MASK_LVS(I,J,K)=0
       PSIF3(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO

! ALPHI_EXT=ALPHI(LVS1)      
      DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)      
      DO NN=1,NUMA(N,LVS1)
      I=I_B(NN,N,LVS1)
      J=J_B(NN,N,LVS1)
      K=K_B(NN,N,LVS1)
      ALPHI_EXT(I,J,K,1)=ALPHI(NN,N,LVS1)
      MASK_LVS(I,J,K)=1
      ENDDO
      ENDDO            
! ALPHI_EXT=ALPHI(LVS2)      
      DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K,TEMP)      
      DO NN=1,NUMA(N,LVS2)
      I=I_B(NN,N,LVS2)
      J=J_B(NN,N,LVS2)
      K=K_B(NN,N,LVS2)
      
      ALPHI2=ALPHI(NN,N,LVS2)
      
      IF (MASK_LVS(I,J,K) .EQ. 0) THEN ! WHEN THERE IS NOTHING IN THE WAY
      ALPHI_EXT(I,J,K,1)=ALPHI2
      MASK_LVS(I,J,K)=1
      ELSE ! (MASK_LVS .EQ. 1)         ! WHEN THERE IS ALPHI(LVS1) IN THE CELL
       IF     ((ALPHI_EXT(I,J,K,1).LT.0.) .AND. (ALPHI2.LT.0.)) THEN      
       ALPHI_EXT(I,J,K,1)=(-1)*MIN(ABS(ALPHI_EXT(I,J,K,1)),ABS(ALPHI2))
       ELSEIF ((ALPHI_EXT(I,J,K,1).GT.0.) .AND. (ALPHI2.LT.0.)) THEN
       ALPHI_EXT(I,J,K,1)=ALPHI2!ALPHI_EXT(I,J,K,1)
       ELSEIF ((ALPHI_EXT(I,J,K,1).LT.0.) .AND. (ALPHI2.GT.0.)) THEN
       ALPHI_EXT(I,J,K,1)=ALPHI_EXT(I,J,K,1)!ALPHI2
       ELSEIF ((ALPHI_EXT(I,J,K,1).GT.0.) .AND. (ALPHI2.GT.0.)) THEN
       ALPHI_EXT(I,J,K,1)=MIN(ABS(ALPHI_EXT(I,J,K,1)),ABS(ALPHI2))!MAX(ALPHI_EXT(I,J,K,1),ALPHI2)
       ELSE
       ALPHI_EXT(I,J,K,1)=ALPHI2
       ENDIF
      ENDIF
      
      ENDDO
      ENDDO 

      
      ! SAVE THE WHOLE PRETTY FIELD IN ALPHI_EXT2
!$OMP PARALLEL DO private(I,J)        
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       ALPHI_EXT2(I,J,K)=ALPHI_EXT(I,J,K,1)
      ENDDO
      ENDDO
      ENDDO     

      
      !EMPTY EXISTING LVS1
       DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LVS1)
         I=I_B(NN,N,LVS1)
         J=J_B(NN,N,LVS1)
         K=K_B(NN,N,LVS1)
          ALPHI(NN,N,LVS1)=0
        ENDDO
       ENDDO  
      !EMPTY EXISTING LVS1 BAND, ACTIVATED CELLS       
       DO N=1,N_MAX
        NUMA(N,LVS1)=0
       ENDDO        
       
       
       ! COUNT CELLS WITHIN D_CRI DISTANCE
       D_CRI=1.2*MAX(SDXF,SDYF,SDZF)
       NUMS_INITIAL=0    
       !2017-08-30 일단 어떻게 될지 모르니 많은 CELL들을 대상으로 시작하겠다
!$OMP PARALLEL DO private(I,J)
!$OMP&reduction(+:NUMS_INITIAL)       
       DO K=1,N3FM
       DO J=1,N2FM
       DO I=1,N1FM
       IF (MASK_LVS(I,J,K) .EQ. 1) NUMS_INITIAL=NUMS_INITIAL+1      
       ENDDO
       ENDDO
       ENDDO
      
      ! 임시적인 BAND 를 만들겠습니다. 
      ! DISTRIBUTION IS ARBITRARY AT THE MOMENT
      NUMS_CRI=NUMS_INITIAL/(N_MAX-1) 
      NUMS_INITIAL=0
      N=1
               
      DO K=1,N3FM
      DO J=1,N2FM
      DO I=1,N1FM
       IF (MASK_LVS(I,J,K).EQ.1) THEN
       ! IF (ABS(ALPHI_EXT(I,J,K,1)) .LE. D_CRI) THEN  !2017-08-30 
       NUMS_INITIAL=NUMS_INITIAL+1           
       I_B(NUMS_INITIAL,N,LVS1)=I
       J_B(NUMS_INITIAL,N,LVS1)=J
       K_B(NUMS_INITIAL,N,LVS1)=K
       IF(NUMS_INITIAL .EQ. NUMS_CRI) THEN
       NUMA(N,LVS1)=NUMS_CRI  ! INITIALLY COUNTED ACTIVATED CELLS PER BAND
       NUMS_INITIAL=0
       N=N+1
       ENDIF
       !ENDIF                                          !2017-08-30 
       ENDIF       
      ENDDO
      ENDDO
      ENDDO
      NUMA(N,LVS1)=NUMS_INITIAL  ! THE REMAININGS
      
      ! CALL BAND_GENERATION(LVS1,NN2)

!$OMP PARALLEL DO private(I,J)
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
      MASK_S(I,J,K)=0
      ENDDO
      ENDDO
      ENDDO
      
!C===============================MAKE S-TUBE============================C
      LLVS=LVS1
      NUMS=0
       DO N=1,NN2
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
       IF(    (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(IPF(I),J,K,1) .LT. 0.)    
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(IMF(I),J,K,1) .LT. 0.)
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,JPF(J),K,1) .LT. 0.)
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,JMF(J),K,1) .LT. 0.)     
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,J,KPF(K),1) .LT. 0.)     
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,J,KMF(K),1) .LT. 0.)
     &                                                            ) THEN
           NUMS=NUMS+1
           I_B(NUMS,1,LLVS)=I
           J_B(NUMS,1,LLVS)=J
           K_B(NUMS,1,LLVS)=K
           MASK_S(I,J,K)=1
       ELSE
           MASK_S(I,J,K)=2
       ENDIF

        ENDDO
       ENDDO
       NUMA(1,LLVS)=NUMS           
!C===============================MAKE S-TUBE============================C
        DO 1003 N=2,N_MAX
          NUMS=0
          NM=N-1
         DO 1004 NN=1,NUMA(NM,LLVS)  !AT S-TUBE
          I=I_B(NN,NM,LLVS)
          J=J_B(NN,NM,LLVS)
          K=K_B(NN,NM,LLVS)  
       IF (MASK_LVS(I,J,K) .EQ. 1) THEN
         DO K2=K-1,K+1
            KK=K2
            IF (IPZ .EQ. 1) THEN
             IF (KK .EQ. 0) KK=N3FM
             IF (KK .EQ. N3F) KK=1  
            ELSE
             IF(K2 .EQ. 1) THEN        !IF INTERFACE TOUCHES THE BOUNDARY
              goto 1005
             ELSE IF(K2 .EQ. N3F) THEN
              goto 1005
             ENDIF
            ENDIF
         DO J2=J-1,J+1
            JJ=J2
            IF (IPY .EQ. 1) THEN
             IF (JJ .EQ. 0) JJ=N2FM
             IF (JJ .EQ. N2F) JJ=1  
            ELSE
             IF(J2 .EQ. 1) THEN 
              goto 1006
             ELSE IF(J2 .EQ. N2F) THEN
              goto 1006
             ENDIF
            ENDIF
          DO I2=I-1,I+1
            II=I2
            IF (IPX .EQ. 1) THEN
             IF (II .EQ. 0) II=N1FM
             IF (II .EQ. N1F) II=1             
            ELSE
             IF(I2 .EQ. 1) THEN 
              goto 1007
             ELSE IF (I2 .EQ. N1FM) THEN
              goto 1007
             ENDIF
            ENDIF
            
       IF ( MASK_S(II,JJ,KK) .NE. 1 ) THEN

        !IF ( MASK_S(II,JJ,KK) .EQ. 2 )  -> ALREADY IN I_B,J_B,K_B & HAVE ALPHI VALUE 

!C-----NEW_CELL ASSIGN
      IF( (MASK_S(II,JJ,KK).EQ. 0) .OR. (N .GE. NN2) ) THEN

      IF ( K2 .EQ. 0) THEN
        IF (J2 .EQ. 0) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ENDIF
        ELSE IF (J2 .EQ. N2F) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ENDIF
        ELSE
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
          ENDIF
        ENDIF
       ELSE IF (K2 .EQ. N3F) THEN
         IF (J2 .EQ. 0) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ENDIF
        ELSE IF (J2 .EQ. N2F) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ENDIF
        ELSE
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-(ZPF(1)-SDZF)
          ENDIF
        ENDIF
        ELSE
        IF (J2 .EQ. 0) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ENDIF
        ELSE IF (J2 .EQ. N2F) THEN
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-(YPF(1)-SDYF)
            DZ=ZPF(KK)-ZPF(K)
          ENDIF
        ELSE
          IF (I2 .EQ. 0) THEN
            DX=XPF(II)-(XPF(N1FM)+SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-ZPF(K)
          ELSE IF (I2 .EQ. N1F) THEN
            DX=XPF(II)-(XPF(1)-SDXF)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-ZPF(K)
          ELSE
            DX=XPF(II)-XPF(I)
            DY=YPF(JJ)-YPF(J)
            DZ=ZPF(KK)-ZPF(K)
          ENDIF
        ENDIF
      ENDIF

       IF (N3FM .EQ. 1) DZ=0. !2D CASE
       IF (ALPHI_EXT(I,J,K,1) .GE. 0.)  THEN
       ALPHI_EXT(II,JJ,KK,1)=ALPHI_EXT(I,J,K,1)+SQRT(DX**2+DY**2+DZ**2) 
       ELSE
       ALPHI_EXT(II,JJ,KK,1)=ALPHI_EXT(I,J,K,1)-SQRT(DX**2+DY**2+DZ**2) 
       ENDIF

       ENDIF !IF( (MASK_S(II,JJ,KK).EQ. 0) .OR. (N .GE. 11) ) THEN

!C-----NEW_CELL ASSIGN
          NUMS=NUMS+1
          I_B(NUMS,N,LLVS)=II
          J_B(NUMS,N,LLVS)=JJ
          K_B(NUMS,N,LLVS)=KK
          MASK_S(II,JJ,KK)=1  

      ENDIF !IF ( MASK_S(II,JJ,KK) .NE. 1 )

 1007   CONTINUE
       ENDDO
 1006   CONTINUE
       ENDDO
 1005   CONTINUE
       ENDDO
       ENDIF
 1004  CONTINUE
           NUMA(N,LLVS)=NUMS
!C==========================makes a I_B,J_B,K_B=========================C
 1003  CONTINUE

       NUMA_MAX=0
       DO N=1,N_MAX
         NUMA_MAX=max(NUMA_MAX,NUMA(N,LLVS))
       ENDDO
!       write(*,*) 'numa_max=',numa_max,mf_band       
!        WRITE(*,*) 'NUMA_MAX/MF_BAND=',FLOAT(NUMA_MAX)/FLOAT(MF_BAND)
       IF (NUMA_MAX .GT. MF_BAND) THEN
        WRITE(*,*) 'NUMA_MAX/MF_BAND=',FLOAT(NUMA_MAX)/FLOAT(MF_BAND)
        WRITE(*,*) 'MF_BAND IS NOT ENOUGH!!'
        WRITE(*,*) 'SIMULATION IS STOPED!!'
        STOP
       ENDIF
       
      
        DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LVS1)
         I=I_B(NN,N,LVS1)
         J=J_B(NN,N,LVS1)
         K=K_B(NN,N,LVS1)
          ALPHI(NN,N,LVS1)=ALPHI_EXT(I,J,K,1)
        ENDDO
       ENDDO
!----------------------------------------------------------------------------!
!-----------------EMPTY LVS2 ALPHI AND NUMA(LVS2)----------------------------!       
       DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LVS2)
         I=I_B(NN,N,LVS2)
         J=J_B(NN,N,LVS2)
         K=K_B(NN,N,LVS2)
          ALPHI(NN,N,LVS2)=0.
        ENDDO
       ENDDO     
       DO N=1,N_MAX
        NUMA(N,LVS2)=0
       ENDDO  
      
       IF (LVS2 .NE. 1) THEN 
       DO L=1,LVS2-1 
      IF (  CONTACT_OLD(LVS1,L) .EQ. 0)    
     &      CONTACT_OLD(LVS1,L)=  CONTACT_OLD(LVS2,L) 
      IF (      CONTACT(LVS1,L) .EQ. 0)        
     &          CONTACT(LVS1,L)=      CONTACT(LVS2,L)
      IF ( CONTACT_TIME(LVS1,L) .EQ. 0.)   
     &     CONTACT_TIME(LVS1,L)= CONTACT_TIME(LVS2,L)
      IF (CONTACT_START(LVS1,L) .EQ. 0.)  
     &    CONTACT_START(LVS1,L)=CONTACT_START(LVS2,L)
      IF (        VELVS(LVS1,L) .EQ. 0.)          
     &            VELVS(LVS1,L)=        VELVS(LVS2,L)
      IF (        T_CON(LVS1,L) .EQ. 0.)          
     &            T_CON(LVS1,L)=        T_CON(LVS2,L)
       ENDDO
       ENDIF
       
       IF (LVS2 .NE. NLVS) THEN
       DO L=LVS2+1,NLVS 
      IF (  CONTACT_OLD(L,LVS1) .EQ. 0)    
     &      CONTACT_OLD(L,LVS1)=  CONTACT_OLD(L,LVS2) 
      IF (      CONTACT(L,LVS1) .EQ. 0)        
     &          CONTACT(L,LVS1)=      CONTACT(L,LVS2)
      IF ( CONTACT_TIME(L,LVS1) .EQ. 0.)   
     &     CONTACT_TIME(L,LVS1)= CONTACT_TIME(L,LVS2)
      IF (CONTACT_START(L,LVS1) .EQ. 0.)  
     &    CONTACT_START(L,LVS1)=CONTACT_START(L,LVS2)
      IF (        VELVS(L,LVS1) .EQ. 0.)          
     &            VELVS(L,LVS1)=        VELVS(L,LVS2)
      IF (        T_CON(L,LVS1) .EQ. 0.)          
     &            T_CON(L,LVS1)=        T_CON(L,LVS2)      
       ENDDO
       ENDIF

       
!----------------------------------------------------------------------------!
!----------------------------------------------------------------------------!
       LVSON(LVS1)=1 
       LVSON(LVS2)=0 
       write(*,*) 'at contact start'
       call cal_psif_mtp(1,nn1,lvs1,psif_mtp)
      write(*,*) 'at contact end'      
       IF (N3FM .EQ. 1) THEN
        DVOL=SDXF*SDYF
       ELSE
        DVOL=SDXF*SDYF*SDZF
       ENDIF

        VOL_TOT_ORI(LVS1)=0.
        VOL_TOT_ORI(LVS2)=0.
        VOL_TEMP=0.
!$OMP PARALLEL DO private(I,J)
!$OMP&reduction(+:VOL_TEMP)         
        DO K=1,N3FM
        DO J=1,N2FM
        DO I=1,N1FM
         VOL_TEMP=VOL_TEMP+(1.-PSIF_MTP(I,J,K))*DVOL
        ENDDO
        ENDDO
        ENDDO
        VOL_TOT_ORI(LVS1)=VOL_TEMP
      WRITE(*,*) 'TOTAL BUBBLE VOL = ',LVS1,VOL_TOT_ORI(LVS1)
      
      ! call cal_psif_mtp(0,nn1,lvs1,psif_mtp)
      
      
      
       CALL CAL_PSI_TMP(0,NN1,LVS1,PSI_TMP,PSIF_MTP,PSI_CN)
      DEALLOCATE(SUR_X,SUR_Y,SUR_Z)
      if (lvs1 .ne. 1) then                             
      ALLOCATE(SUR_X(M1M,M2M,M3M),SUR_Y(M1M,M2M,M3M),SUR_Z(M1M,M2M,M3M))                                  
      endif             
       CALL CAL_CUR(4,LVS1,PSI_TMP,PSIF_MTP)       
      
      
      
      ! ONLY CHANGE THE PSIF OF COALESCED BANDS
      DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LVS1)
         I=I_B(NN,N,LVS1)
         J=J_B(NN,N,LVS1)
         K=K_B(NN,N,LVS1)
         PSIF(I,J,K)=PSIF_MTP(I,J,K)
        ENDDO
       ENDDO  
      
!$OMP PARALLEL DO private(I,J)      
       DO K=1,N3FM
       DO J=1,N2FM
       DO I=1,N1FM
          IF (PSIF(I,J,K) .LT. 1.E-3) THEN !FOR NUMERICAL ERROR
            PSIF(I,J,K)=0. 
          ELSE IF (PSIF(I,J,K) .GT. 0.999) THEN
            PSIF(I,J,K)=1.
          ENDIF
       ENDDO
       ENDDO
       ENDDO  

 
!!2017-08-19_ to see how alphi changes according to band#...
      ! DO N=1,11!(change)
        ! DO NN=1,NUMA(N,LVS1)
         ! I=I_B(NN,N,LVS1)
         ! J=J_B(NN,N,LVS1)
         ! K=K_B(NN,N,LVS1)
        ! alphi_ext(i,j,k,1)=ALPHI(NN,N,LVS1)
        ! ENDDO
       ! ENDDO
       
       ! OPEN(146,FILE='LVS1alphi_contour_AFTER_BANDG.DAT')
       ! WRITE(146,*) 'VARIABLES="X","Y","Z","alphi"'
      ! WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
      ! DO K=1,N3FM
      ! DO J=1,N2FM
      ! DO I=1,N1FM
        ! WRITE(146,147) XPF(I),YPF(J),ZPF(K),ALPHI_EXT(i,j,k,1)
      ! ENDDO
      ! ENDDO
      ! ENDDO          
       ! CLOSE(146)
 ! 147  FORMAT(4F15.8)
      ! stop
       ! DO N=1,nn1
         ! DO NN=1,NUMA(N,LVS1)
          ! I=I_B(NN,N,LVS1)
          ! J=J_B(NN,N,LVS1)
          ! K=K_B(NN,N,LVS1)
         ! psif3(i,j,k)=PSIF(I,J,K)
         ! ENDDO
        ! ENDDO         
        
        ! OPEN(148,FILE='03PSIF.DAT')
        ! WRITE(148,*) 'VARIABLES="X","Y","Z","psif"'
       ! WRITE(148,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
       ! DO K=1,N3FM
       ! DO J=1,N2FM
       ! DO I=1,N1FM
         ! WRITE(148,149) XPF(I),YPF(J),ZPF(K),psif(i,j,k)
       ! ENDDO
       ! ENDDO
       ! ENDDO
        ! CLOSE(148)
  ! 149  FORMAT(4F15.8)       
      
      RETURN
      END           

!***********************************************************************
      SUBROUTINE CAL_TCON(L,LL,vol_tot_ori)
!***********************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE TWO_PHASE_PROPERTY
      
      IMPLICIT NONE
      REAL VOL_TOT_ORI(MLVS)     
      REAL UREL
      REAL UU1,UU2,VV1,VV2,WW1,WW2
      INTEGER L,LL
      real length
      real dchar,vols
      ! length=0.4  !
      
      UU1=VELVS(L,1)
      VV1=VELVS(L,2)
      WW1=VELVS(L,3)
             
      UU2=VELVS(LL,1)
      VV2=VELVS(LL,2)
      WW2=VELVS(LL,3)
      
      vols=min(vol_tot_ori(ll),vol_tot_ori(l))
      dchar=0.5*((6*vols/3.14159)**(1./3.))
      UREL=((UU1-UU2)**2+(VV1-VV2)**2+(WW1-WW2)**2)**0.5
      !denr=5.8
      !T_CON(L,LL)=DENR*SURF_J*length*urel
        T_CON(L,LL)=SURF_J*(urel**2)*dchar                                
      
       write(*,*) uu1,uu2
       write(*,*) vv1,vv2
       write(*,*) ww1,ww2
       write(*,*) urel,surf_j,t_con(l,ll)
       !t_con(l,ll)=0.0001
      
      RETURN
      END        
!!*******************************************************************
!      SUBROUTINE CAL_COALESCENCE(LVS2,LVS1,PSI_CN,VOL_TOT_ORI)
!!*******************************************************************   
!      USE PARAM_VAR
!      USE FLOW_GEOM_VAR
!      USE FLOW_VAR
!
!      USE LVS_VAR
!      USE LVS_GEOM_VAR
!      
!      IMPLICIT NONE
!      
!      REAL ALPHI_C(M1FM,M2FM,M3FM)      
!      REAL D_CRI,TEMP                                         
!      REAL DX,DY,DZ,DVOL 
!      REAL XX,YY,ZZ,ALPHI_EXT2(-2:M1F+3,-2:M2F+3,-2:M3F+3)     
!      REAL VOL_TOT_ORI(MLVS)
!      REAL PSI_TMP(0:M1,0:M2,0:M3)
!      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)
!      REAL PSI_CN(0:M1,0:M2,0:M3)
!      REAL psif3(M1L:M1U,M2L:M2U,M3L:M3U)
!      REAL ALPHI2
!      
!      INTEGER MASK_LVS(-2:M1F+3,-2:M2F+3,-2:M3F+3)
!      INTEGER MASK_S(-2:M1F+3,-2:M2F+3,-2:M3F+3)   
!      INTEGER NN1,NN2,LVS1,LVS2,NM,LLVS      
!      INTEGER I,J,K,N,NN,II,JJ,KK,I2,J2,K2
!      INTEGER NUMS,NUMS_INITIAL,NUMS_CRI
!      
!      nn1=n_max-3
!      nn2=n_max-1
!! INITIALIZE THE MASKS, AND WORKING ALPHI_EXT
!!$OMP PARALLEL DO            
!      DO K=-2,N3F+3
!      DO J=-2,N2F+3
!      DO I=-2,N1F+3
!       ALPHI_EXT(I,J,K,1)=0.
!       MASK_LVS(I,J,K)=0
!       PSIF3(I,J,K)=0.
!      ENDDO
!      ENDDO
!      ENDDO
!
!! ALPHI_EXT=ALPHI(LVS1)      
!      DO N=1,N_MAX
!!$OMP PARALLEL DO private(I,J,K)      
!      DO NN=1,NUMA(N,LVS1)
!      I=I_B(NN,N,LVS1)
!      J=J_B(NN,N,LVS1)
!      K=K_B(NN,N,LVS1)
!      ALPHI_EXT(I,J,K,1)=ALPHI(NN,N,LVS1)
!      MASK_LVS(I,J,K)=1
!      ENDDO
!      ENDDO            
!! ALPHI_EXT=ALPHI(LVS2)      
!      DO N=1,N_MAX
!!$OMP PARALLEL DO private(I,J,K,TEMP)      
!      DO NN=1,NUMA(N,LVS2)
!      I=I_B(NN,N,LVS2)
!      J=J_B(NN,N,LVS2)
!      K=K_B(NN,N,LVS2)
!      
!      ALPHI2=ALPHI(NN,N,LVS2)
!      
!      IF (MASK_LVS(I,J,K) .EQ. 0) THEN ! WHEN THERE IS NOTHING IN THE WAY
!      ALPHI_EXT(I,J,K,1)=ALPHI2
!      MASK_LVS(I,J,K)=1
!      ELSE ! (MASK_LVS .EQ. 1)         ! WHEN THERE IS ALPHI(LVS1) IN THE CELL
!       IF     ((ALPHI_EXT(I,J,K,1).LT.0.) .AND. (ALPHI2.LT.0.)) THEN      
!       ALPHI_EXT(I,J,K,1)=(-1)*MIN(ABS(ALPHI_EXT(I,J,K,1)),ABS(ALPHI2))
!       ELSEIF ((ALPHI_EXT(I,J,K,1).GT.0.) .AND. (ALPHI2.LT.0.)) THEN
!       ALPHI_EXT(I,J,K,1)=ALPHI2!ALPHI_EXT(I,J,K,1)
!       ELSEIF ((ALPHI_EXT(I,J,K,1).LT.0.) .AND. (ALPHI2.GT.0.)) THEN
!       ALPHI_EXT(I,J,K,1)=ALPHI_EXT(I,J,K,1)!ALPHI2
!       ELSEIF ((ALPHI_EXT(I,J,K,1).GT.0.) .AND. (ALPHI2.GT.0.)) THEN
!       ALPHI_EXT(I,J,K,1)=MIN(ABS(ALPHI_EXT(I,J,K,1)),ABS(ALPHI2))!MAX(ALPHI_EXT(I,J,K,1),ALPHI2)
!       ELSE
!       ALPHI_EXT(I,J,K,1)=ALPHI2
!       ENDIF
!      ENDIF
!      
!      ENDDO
!      ENDDO 
!
!      
!      ! SAVE THE WHOLE PRETTY FIELD IN ALPHI_EXT2
!!$OMP PARALLEL DO private(I,J)        
!      DO K=-2,N3F+3
!      DO J=-2,N2F+3
!      DO I=-2,N1F+3
!       ALPHI_EXT2(I,J,K)=ALPHI_EXT(I,J,K,1)
!      ENDDO
!      ENDDO
!      ENDDO     
!
!      
!      !EMPTY EXISTING LVS1
!       DO N=1,N_MAX
!!$OMP PARALLEL DO private(I,J,K)
!        DO NN=1,NUMA(N,LVS1)
!         I=I_B(NN,N,LVS1)
!         J=J_B(NN,N,LVS1)
!         K=K_B(NN,N,LVS1)
!          ALPHI(NN,N,LVS1)=0
!        ENDDO
!       ENDDO  
!      !EMPTY EXISTING LVS1 BAND, ACTIVATED CELLS       
!       DO N=1,N_MAX
!        NUMA(N,LVS1)=0
!       ENDDO        
!       
!       
!       ! COUNT CELLS WITHIN D_CRI DISTANCE
!       D_CRI=1.2*MAX(SDXF,SDYF,SDZF)
!       NUMS_INITIAL=0
!!!$OMP PARALLEL DO 
!!!$OMP&REDUCTION(+:NUMS_INITIAL)            
!       ! DO K=1,N3FM
!       ! DO J=1,N2FM
!       ! DO I=1,N1FM
!       ! IF (MASK_LVS(I,J,K) .EQ. 1) THEN
!       ! IF (ABS(ALPHI_EXT(I,J,K,1)).LE.D_CRI) NUMS_INITIAL=NUMS_INITIAL+1
!       ! ENDIF       
!       ! ENDDO
!       ! ENDDO
!       ! ENDDO     
!       !2017-08-30 일단 어떻게 될지 모르니 많은 CELL들을 대상으로 시작하겠다
!!$OMP PARALLEL DO private(I,J)
!!$OMP&reduction(+:NUMS_INITIAL)       
!       DO K=1,N3FM
!       DO J=1,N2FM
!       DO I=1,N1FM
!       IF (MASK_LVS(I,J,K) .EQ. 1) NUMS_INITIAL=NUMS_INITIAL+1      
!       ENDDO
!       ENDDO
!       ENDDO
!      
!      ! 임시적인 BAND 를 만들겠습니다. 
!      ! DISTRIBUTION IS ARBITRARY AT THE MOMENT
!      NUMS_CRI=NUMS_INITIAL/(N_MAX-1) 
!      NUMS_INITIAL=0
!      N=1
!               
!      DO K=1,N3FM
!      DO J=1,N2FM
!      DO I=1,N1FM
!       IF (MASK_LVS(I,J,K).EQ.1) THEN
!       ! IF (ABS(ALPHI_EXT(I,J,K,1)) .LE. D_CRI) THEN  !2017-08-30 
!       NUMS_INITIAL=NUMS_INITIAL+1           
!       I_B(NUMS_INITIAL,N,LVS1)=I
!       J_B(NUMS_INITIAL,N,LVS1)=J
!       K_B(NUMS_INITIAL,N,LVS1)=K
!       IF(NUMS_INITIAL .EQ. NUMS_CRI) THEN
!       NUMA(N,LVS1)=NUMS_CRI  ! INITIALLY COUNTED ACTIVATED CELLS PER BAND
!       NUMS_INITIAL=0
!       N=N+1
!       ENDIF
!       !ENDIF                                          !2017-08-30 
!       ENDIF       
!      ENDDO
!      ENDDO
!      ENDDO
!      NUMA(N,LVS1)=NUMS_INITIAL  ! THE REMAININGS
!      
!      ! CALL BAND_GENERATION(LVS1,NN2)
!
!!$OMP PARALLEL DO private(I,J)
!      DO K=-2,N3F+3
!      DO J=-2,N2F+3
!      DO I=-2,N1F+3
!      MASK_S(I,J,K)=0
!      ENDDO
!      ENDDO
!      ENDDO
!      
!!C===============================MAKE S-TUBE============================C
!      LLVS=LVS1
!      NUMS=0
!       DO N=1,NN2
!        DO NN=1,NUMA(N,LLVS)
!         I=I_B(NN,N,LLVS)
!         J=J_B(NN,N,LLVS)
!         K=K_B(NN,N,LLVS)
!       IF(    (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(IPF(I),J,K,1) .LT. 0.)    
!     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(IMF(I),J,K,1) .LT. 0.)
!     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,JPF(J),K,1) .LT. 0.)
!     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,JMF(J),K,1) .LT. 0.)     
!     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,J,KPF(K),1) .LT. 0.)     
!     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,J,KMF(K),1) .LT. 0.)
!     &                                                            ) THEN
!           NUMS=NUMS+1
!           I_B(NUMS,1,LLVS)=I
!           J_B(NUMS,1,LLVS)=J
!           K_B(NUMS,1,LLVS)=K
!           MASK_S(I,J,K)=1
!       ELSE
!           MASK_S(I,J,K)=2
!       ENDIF
!
!        ENDDO
!       ENDDO
!       NUMA(1,LLVS)=NUMS           
!!C===============================MAKE S-TUBE============================C
!        DO 1003 N=2,N_MAX
!          NUMS=0
!          NM=N-1
!         DO 1004 NN=1,NUMA(NM,LLVS)  !AT S-TUBE
!          I=I_B(NN,NM,LLVS)
!          J=J_B(NN,NM,LLVS)
!          K=K_B(NN,NM,LLVS)  
!       IF (MASK_LVS(I,J,K) .EQ. 1) THEN
!         DO K2=K-1,K+1
!            KK=K2
!            IF (IPZ .EQ. 1) THEN
!             IF (KK .EQ. 0) KK=N3FM
!             IF (KK .EQ. N3F) KK=1  
!            ELSE
!             IF(K2 .EQ. 1) THEN        !IF INTERFACE TOUCHES THE BOUNDARY
!              goto 1005
!             ELSE IF(K2 .EQ. N3F) THEN
!              goto 1005
!             ENDIF
!            ENDIF
!         DO J2=J-1,J+1
!            JJ=J2
!            IF (IPY .EQ. 1) THEN
!             IF (JJ .EQ. 0) JJ=N2FM
!             IF (JJ .EQ. N2F) JJ=1  
!            ELSE
!             IF(J2 .EQ. 1) THEN 
!              goto 1006
!             ELSE IF(J2 .EQ. N2F) THEN
!              goto 1006
!             ENDIF
!            ENDIF
!          DO I2=I-1,I+1
!            II=I2
!            IF (IPX .EQ. 1) THEN
!            ELSE
!             IF(I2 .EQ. 1) THEN 
!              goto 1007
!             ELSE IF (I2 .EQ. N1FM) THEN
!              goto 1007
!             ENDIF
!            ENDIF
!            
!       IF ( MASK_S(II,JJ,KK) .NE. 1 ) THEN
!
!        !IF ( MASK_S(II,JJ,KK) .EQ. 2 )  -> ALREADY IN I_B,J_B,K_B & HAVE ALPHI VALUE 
!
!!C-----NEW_CELL ASSIGN
!      IF( (MASK_S(II,JJ,KK).EQ. 0) .OR. (N .GE. NN2) ) THEN
!
!      IF ( K2 .EQ. 0) THEN
!        IF (J2 .EQ. 0) THEN
!          IF (I2 .EQ. 0) THEN
!            DX=XPF(II)-(XPF(N1FM)+SDXF)
!            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
!            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
!          ELSE IF (I2 .EQ. N1F) THEN
!            DX=XPF(II)-(XPF(1)-SDXF)
!            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
!            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
!          ELSE
!            DX=XPF(II)-XPF(I)
!            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
!            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
!          ENDIF
!        ELSE IF (J2 .EQ. N2F) THEN
!          IF (I2 .EQ. 0) THEN
!            DX=XPF(II)-(XPF(N1FM)+SDXF)
!            DY=YPF(JJ)-(YPF(1)-SDYF)
!            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
!          ELSE IF (I2 .EQ. N1F) THEN
!            DX=XPF(II)-(XPF(1)-SDXF)
!            DY=YPF(JJ)-(YPF(1)-SDYF)
!            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
!          ELSE
!            DX=XPF(II)-XPF(I)
!            DY=YPF(JJ)-(YPF(1)-SDYF)
!            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
!          ENDIF
!        ELSE
!          IF (I2 .EQ. 0) THEN
!            DX=XPF(II)-(XPF(N1FM)+SDXF)
!            DY=YPF(JJ)-YPF(J)
!            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
!          ELSE IF (I2 .EQ. N1F) THEN
!            DX=XPF(II)-(XPF(1)-SDXF)
!            DY=YPF(JJ)-YPF(J)
!            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
!          ELSE
!            DX=XPF(II)-XPF(I)
!            DY=YPF(JJ)-YPF(J)
!            DZ=ZPF(KK)-(ZPF(N3FM)+SDZF)
!          ENDIF
!        ENDIF
!       ELSE IF (K2 .EQ. N3F) THEN
!         IF (J2 .EQ. 0) THEN
!          IF (I2 .EQ. 0) THEN
!            DX=XPF(II)-(XPF(N1FM)+SDXF)
!            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
!            DZ=ZPF(KK)-(ZPF(1)-SDZF)
!          ELSE IF (I2 .EQ. N1F) THEN
!            DX=XPF(II)-(XPF(1)-SDXF)
!            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
!            DZ=ZPF(KK)-(ZPF(1)-SDZF)
!          ELSE
!            DX=XPF(II)-XPF(I)
!            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
!            DZ=ZPF(KK)-(ZPF(1)-SDZF)
!          ENDIF
!        ELSE IF (J2 .EQ. N2F) THEN
!          IF (I2 .EQ. 0) THEN
!            DX=XPF(II)-(XPF(N1FM)+SDXF)
!            DY=YPF(JJ)-(YPF(1)-SDYF)
!            DZ=ZPF(KK)-(ZPF(1)-SDZF)
!          ELSE IF (I2 .EQ. N1F) THEN
!            DX=XPF(II)-(XPF(1)-SDXF)
!            DY=YPF(JJ)-(YPF(1)-SDYF)
!            DZ=ZPF(KK)-(ZPF(1)-SDZF)
!          ELSE
!            DX=XPF(II)-XPF(I)
!            DY=YPF(JJ)-(YPF(1)-SDYF)
!            DZ=ZPF(KK)-(ZPF(1)-SDZF)
!          ENDIF
!        ELSE
!          IF (I2 .EQ. 0) THEN
!            DX=XPF(II)-(XPF(N1FM)+SDXF)
!            DY=YPF(JJ)-YPF(J)
!            DZ=ZPF(KK)-(ZPF(1)-SDZF)
!          ELSE IF (I2 .EQ. N1F) THEN
!            DX=XPF(II)-(XPF(1)-SDXF)
!            DY=YPF(JJ)-YPF(J)
!            DZ=ZPF(KK)-(ZPF(1)-SDZF)
!          ELSE
!            DX=XPF(II)-XPF(I)
!            DY=YPF(JJ)-YPF(J)
!            DZ=ZPF(KK)-(ZPF(1)-SDZF)
!          ENDIF
!        ENDIF
!        ELSE
!        IF (J2 .EQ. 0) THEN
!          IF (I2 .EQ. 0) THEN
!            DX=XPF(II)-(XPF(N1FM)+SDXF)
!            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
!            DZ=ZPF(KK)-ZPF(K)
!          ELSE IF (I2 .EQ. N1F) THEN
!            DX=XPF(II)-(XPF(1)-SDXF)
!            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
!            DZ=ZPF(KK)-ZPF(K)
!          ELSE
!            DX=XPF(II)-XPF(I)
!            DY=YPF(JJ)-(YPF(N2FM)+SDYF)
!            DZ=ZPF(KK)-ZPF(K)
!          ENDIF
!        ELSE IF (J2 .EQ. N2F) THEN
!          IF (I2 .EQ. 0) THEN
!            DX=XPF(II)-(XPF(N1FM)+SDXF)
!            DY=YPF(JJ)-(YPF(1)-SDYF)
!            DZ=ZPF(KK)-ZPF(K)
!          ELSE IF (I2 .EQ. N1F) THEN
!            DX=XPF(II)-(XPF(1)-SDXF)
!            DY=YPF(JJ)-(YPF(1)-SDYF)
!            DZ=ZPF(KK)-ZPF(K)
!          ELSE
!            DX=XPF(II)-XPF(I)
!            DY=YPF(JJ)-(YPF(1)-SDYF)
!            DZ=ZPF(KK)-ZPF(K)
!          ENDIF
!        ELSE
!          IF (I2 .EQ. 0) THEN
!            DX=XPF(II)-(XPF(N1FM)+SDXF)
!            DY=YPF(JJ)-YPF(J)
!            DZ=ZPF(KK)-ZPF(K)
!          ELSE IF (I2 .EQ. N1F) THEN
!            DX=XPF(II)-(XPF(1)-SDXF)
!            DY=YPF(JJ)-YPF(J)
!            DZ=ZPF(KK)-ZPF(K)
!          ELSE
!            DX=XPF(II)-XPF(I)
!            DY=YPF(JJ)-YPF(J)
!            DZ=ZPF(KK)-ZPF(K)
!          ENDIF
!        ENDIF
!      ENDIF
!
!       IF (N3FM .EQ. 1) DZ=0. !2D CASE
!       IF (ALPHI_EXT(I,J,K,1) .GE. 0.)  THEN
!       ALPHI_EXT(II,JJ,KK,1)=ALPHI_EXT(I,J,K,1)+SQRT(DX**2+DY**2+DZ**2) 
!       ELSE
!       ALPHI_EXT(II,JJ,KK,1)=ALPHI_EXT(I,J,K,1)-SQRT(DX**2+DY**2+DZ**2) 
!       ENDIF
!
!       ENDIF !IF( (MASK_S(II,JJ,KK).EQ. 0) .OR. (N .GE. 11) ) THEN
!
!!C-----NEW_CELL ASSIGN
!          NUMS=NUMS+1
!          I_B(NUMS,N,LLVS)=II
!          J_B(NUMS,N,LLVS)=JJ
!          K_B(NUMS,N,LLVS)=KK
!          MASK_S(II,JJ,KK)=1  
!
!      ENDIF !IF ( MASK_S(II,JJ,KK) .NE. 1 )
!
! 1007   CONTINUE
!       ENDDO
! 1006   CONTINUE
!       ENDDO
! 1005   CONTINUE
!       ENDDO
!       ENDIF
! 1004  CONTINUE
!           NUMA(N,LLVS)=NUMS
!!C==========================makes a I_B,J_B,K_B=========================C
! 1003  CONTINUE
!
!       NUMA_MAX=0
!       DO N=1,N_MAX
!         NUMA_MAX=max(NUMA_MAX,NUMA(N,LLVS))
!       ENDDO
!!       write(*,*) 'numa_max=',numa_max,mf_band       
!!        WRITE(*,*) 'NUMA_MAX/MF_BAND=',FLOAT(NUMA_MAX)/FLOAT(MF_BAND)
!       IF (NUMA_MAX .GT. MF_BAND) THEN
!        WRITE(*,*) 'NUMA_MAX/MF_BAND=',FLOAT(NUMA_MAX)/FLOAT(MF_BAND)
!        WRITE(*,*) 'MF_BAND IS NOT ENOUGH!!'
!        WRITE(*,*) 'SIMULATION IS STOPED!!'
!        STOP
!       ENDIF
!       
!      
!        DO N=1,N_MAX
!!$OMP PARALLEL DO private(I,J,K)
!        DO NN=1,NUMA(N,LVS1)
!         I=I_B(NN,N,LVS1)
!         J=J_B(NN,N,LVS1)
!         K=K_B(NN,N,LVS1)
!          ALPHI(NN,N,LVS1)=ALPHI_EXT(I,J,K,1)
!        ENDDO
!       ENDDO
!!----------------------------------------------------------------------------!
!!-----------------EMPTY LVS2 ALPHI AND NUMA(LVS2)----------------------------!       
!       DO N=1,N_MAX
!!$OMP PARALLEL DO private(I,J,K)
!        DO NN=1,NUMA(N,LVS2)
!         I=I_B(NN,N,LVS2)
!         J=J_B(NN,N,LVS2)
!         K=K_B(NN,N,LVS2)
!          ALPHI(NN,N,LVS2)=0.
!        ENDDO
!       ENDDO     
!       DO N=1,N_MAX
!        NUMA(N,LVS2)=0
!       ENDDO       
!!----------------------------------------------------------------------------!
!!----------------------------------------------------------------------------!
!       LVSON(LVS1)=1 
!       LVSON(LVS2)=0 
! 
!!!2017-08-19_ to see how alphi changes according to band#...
!      ! DO N=1,11!(change)
!        ! DO NN=1,NUMA(N,LVS1)
!         ! I=I_B(NN,N,LVS1)
!         ! J=J_B(NN,N,LVS1)
!         ! K=K_B(NN,N,LVS1)
!        ! alphi_ext(i,j,k,1)=ALPHI(NN,N,LVS1)
!        ! ENDDO
!       ! ENDDO
!       
!       ! OPEN(146,FILE='LVS1alphi_contour_AFTER_BANDG.DAT')
!       ! WRITE(146,*) 'VARIABLES="X","Y","Z","alphi"'
!      ! WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!      ! DO K=1,N3FM
!      ! DO J=1,N2FM
!      ! DO I=1,N1FM
!        ! WRITE(146,147) XPF(I),YPF(J),ZPF(K),ALPHI_EXT(i,j,k,1)
!      ! ENDDO
!      ! ENDDO
!      ! ENDDO          
!       ! CLOSE(146)
! ! 147  FORMAT(4F15.8)
!      ! stop
!       
!       
!       call cal_psif_mtp(1,nn1,lvs1,psif_mtp)
!       
!       DO K=1,N3FM
!       DO J=1,N2FM
!       DO I=1,N1FM
!         PSIF(I,J,K)=PSIF_MTP(I,J,K)        
!       ENDDO
!       ENDDO
!       ENDDO
!       
!       
!       ! DO N=1,nn1
!         ! DO NN=1,NUMA(N,LVS1)
!          ! I=I_B(NN,N,LVS1)
!          ! J=J_B(NN,N,LVS1)
!          ! K=K_B(NN,N,LVS1)
!         ! psif3(i,j,k)=PSIF(I,J,K)
!         ! ENDDO
!        ! ENDDO         
!        
!        ! OPEN(148,FILE='03PSIF.DAT')
!        ! WRITE(148,*) 'VARIABLES="X","Y","Z","psif"'
!       ! WRITE(148,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
!       ! DO K=1,N3FM
!       ! DO J=1,N2FM
!       ! DO I=1,N1FM
!         ! WRITE(148,149) XPF(I),YPF(J),ZPF(K),psif(i,j,k)
!       ! ENDDO
!       ! ENDDO
!       ! ENDDO
!        ! CLOSE(148)
!  ! 149  FORMAT(4F15.8)
!     
!  
!      !CALL CAL_PSIF_MTP(0,NN1,LVS1,PSIF_MTP)  
!
!
!      
!       IF (N3FM .EQ. 1) THEN
!        DVOL=SDXF*SDYF
!       ELSE
!        DVOL=SDXF*SDYF*SDZF
!       ENDIF
!        VOL_TOT_ORI(LVS1)=0.
!        DO K=1,N3FM
!        DO J=1,N2FM
!        DO I=1,N1FM
!         VOL_TOT_ORI(LVS1)=VOL_TOT_ORI(LVS1)+(1.-PSIF_MTP(I,J,K))*DVOL
!        ENDDO
!        ENDDO
!        ENDDO
!      WRITE(*,*) 'TOTAL BUBBLE VOL = ',LVS1,VOL_TOT_ORI(LVS1)
!      
!      
!      CALL CAL_PSI_TMP(0,NN1,LVS1,PSI_TMP,PSIF_MTP,PSI_CN)
!
!!$OMP PARALLEL DO private(I,J)      
!       DO K=1,N3FM
!       DO J=1,N2FM
!       DO I=1,N1FM
!          IF (PSIF(I,J,K) .LT. 1.E-3) THEN !FOR NUMERICAL ERROR
!            PSIF(I,J,K)=0. 
!          ELSE IF (PSIF(I,J,K) .GT. 0.999) THEN
!            PSIF(I,J,K)=1.
!          ENDIF
!       ENDDO
!       ENDDO
!       ENDDO  
!       DEALLOCATE(SUR_X,SUR_Y,SUR_Z)       
!       CALL CAL_CUR(4,LVS1,PSI_TMP,PSIF_MTP) 
!
!       
!      
!      RETURN
!      END           
!
!!***********************************************************************
!      SUBROUTINE CAL_TCON(L,LL)
!!***********************************************************************
!      USE PARAM_VAR
!      USE FLOW_GEOM_VAR
!      USE FLOW_VAR
!
!      USE LVS_VAR
!      USE LVS_GEOM_VAR
!      USE TWO_PHASE_PROPERTY
!      
!      IMPLICIT NONE
!      
!      REAL UREL
!      REAL UU1,UU2,VV1,VV2,WW1,WW2
!      INTEGER L,LL
!            
!      UU1=VELVS(L,1)
!      VV1=VELVS(L,2)
!      WW1=VELVS(L,3)
!             
!      UU2=VELVS(LL,1)
!      VV2=VELVS(LL,2)
!      WW2=VELVS(LL,3)
!      
!      
!      UREL=((UU1-UU2)**2+(VV1-VV2)**2+(WW1-WW2)**2)**0.5
!      
!      T_CON(L,LL)=DENR*UREL*(0.5)**2*SURF_J
!      
!      ! write(*,*) uu1,uu2
!      ! write(*,*) vv1,vv2
!      ! write(*,*) ww1,ww2
!      ! write(*,*) urel, t_con(l,ll)
!      
!      RETURN
!      END
!      
      