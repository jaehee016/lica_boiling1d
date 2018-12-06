!***********************************************************************
      SUBROUTINE SHARP_DIS(DH)
!***********************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE HEAT_VAR
      
      IMPLICIT NONE
      INTEGER*8   I,J,K
      REAL*8      U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)    
      REAL*8      PSI_CN(0:M1,0:M2,0:M3)  
      REAL*8      ANPSO(M1M,M2M,M3M),DENSC,DENSCF  
      REAL*8      DH(0:M1,0:M2,0:M3,3)    

     
      
      
      
      RETURN
      END
      
!***********************************************************************
      SUBROUTINE LSHPS2(PSI_CN,PSI_C,DF_C)     
!***********************************************************************
      IMPLICIT NONE
      INTEGER*8   I,J,K
      REAL*8      U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)    
      REAL*8      PSI_CN(0:M1,0:M2,0:M3)  
      REAL*8      ANPSO(M1M,M2M,M3M),DENSC,DENSCF  
      REAL*8      DH(0:M1,0:M2,0:M3,3)    

     ! NEW LHSINIT COEFFICIENTS ARE NEEDE HERE (EVERY TIMESTEP)
      
      
      
      RETURN
      END







!***********************************************************************
      SUBROUTINE TEMP_GRAD(PSI_CN)
!***********************************************************************
! PARAMETERS USED
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
       use flow_var ! this is to save files at certain steps                                                          

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      USE HEAT_VAR
      USE PHASE_BOILING
      use lvs_coupling
      
      IMPLICIT NONE
      
      REAL*8      PSI_CN(0:M1,0:M2,0:M3)  
      REAL*8      PSI_TMP(0:M1,0:M2,0:M3)
      
      REAL*8      NPHI(0:M1,0:M2,0:M3)             
      REAL*8      PSIH(0:M1,0:M2,0:M3) 
!      REAL*8      TEMP(0:M1,0:M2,0:M3,0:1) 
      INTEGER*8   IC
      REAL*8      NSPHI,NPHI_PRE,DX,DY,DZ
      
      REAL*8      TF(MF_BAND,M_MAX,MLVS,0:1)
      REAL*8      PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)
      REAL*8      check0(M1L:M1U,M2L:M2U,M3L:M3U)
      REAL*8      check1(M1L:M1U,M2L:M2U,M3L:M3U)
      INTEGER*8   NN1,NN2,N,NN,I,J,K,LLVS,IM1,JM1,KM1
      REAL*8      T_INTP,MDOT
      INTEGER*8   NGRAD,kk,jj,ii
      
       REAL*8      MFDOT(M1L:M1U,M2L:M2U,M3L:M3U)  
       
      ALLOCATE(MFLVS(MF_BAND,M_MAX,MLVS))
      ALLOCATE(SGLVS(MF_BAND,M_MAX,MLVS))
      ! ALLOCATE(IMASK_PH(0:M1,0:M2,0:M3))
      
      check0=0d0
      check1=0d0
      NGRAD=N_MAX-6
      
      nn1=n_max-3
      nn2=n_max-1
      IMASK_PH=0
      
!INITIALISE     
       DO LLVS=1,NLVS 
       DO N=1,N_MAX
!$OMP PARALLEL DO
       DO NN=1,NUMA(N,LLVS)
       MFLVS(NN,N,LLVS)=0D0
       SGLVS(NN,N,LLVS)=0D0
       ENDDO
       ENDDO
       ENDDO
       
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3
      DO J=0,N2
      DO I=0,N1
      TEMP(I,J,K,0)=T(I,J,K)
      TEMP(I,J,K,1)=T(I,J,K)
      NPHI(I,J,K)=0D0
      PSIH(I,J,K)=0D0
      ENDDO
      ENDDO
      ENDDO
      
 !!=====SAVE==============================================================
! TWO_PHASE.F LVSINIT 부분 ALPHI_EXT2 AFTER BANDGENERATION 주석 부분 제거하고
! 여기서 출력되는 것과 어떻게 다른지 비교해보기 (1D CONFIGURATION에서는 LINEAR해서 문제 없음)
        OPEN(198,FILE='0Tprofile.DAT')
       !WRITE(198,*) 'VARIABLES="X","Y","Z","t"'
       WRITE(198,*) 'VARIABLES="Z","t"'
      !WRITE(198,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
      WRITE(198,*)'ZONE K=',N3M,',F=POINT'
      DO K=1,N3M
      !DO J=1,N2M
      !DO I=1,N1M
      ! WRITE(198,197)XP(I),YP(J),ZP(K),t(I,J,K)
      WRITE(198,197)ZP(K),t(1,1,K)
      !ENDDO
      ! ENDDO
      ENDDO
       CLOSE(198)
! 197  FORMAT(3F16.8,1D20.8)
 197  FORMAT(1F16.8,1D20.8)

       !STOP
!!=======================================================================
       DO 774 LLVS=1,NLVS  
      
!(0) ALLOCATE ALPHI_EXT OF THE CURRENT LVS++++++++++++++++++++++++++++++
       ALLOCATE(ALPHI_EXT(-2:M1F+3,-2:M2F+3,-2:M3F+3,1))
!$OMP PARALLEL DO            
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       ALPHI_EXT(I,J,K,1)=0D0
      ENDDO
      ENDDO
      ENDDO
       DO N=1,N_MAX
!$OMP PARALLEL DO private(i,j,k)
       DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)       
       ALPHI_EXT(I,J,K,1)=ALPHI(NN,N,LLVS)
       ENDDO
       ENDDO
       
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!(1) DETERMINE THE PHASE OF THE NAVIER STOKES CELLS: PSIH(II,JJ,KK)+++++     
!nn2-1, 2배 refined 일때: NSgrid 왼5개, 오른4개가 제대로 계산됐고,
!nn2,   2배 refined 일때: NSgrid 왼6개, 오른5개가 계산됐는데, 
!                      왼1개는 보간에 0을 가져다 사용해서 더 작은 값이 나옴. 
       DO N=1,NN2      ! LVS AT THE BRINK CAN'T USE ALPHI OUTSIDE THE BANDS          
!$OMP PARALLEL DO private(I,J,K,T_INTP,ii,jj,kk,NSPHI,IC)
       DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         II=ICOUMP_VEL(I)
         JJ=JCOUMP_VEL(J)
         KK=KCOUMP_VEL(K)
        IF (NPHI(II,JJ,KK) .EQ. 0D0) THEN 
         CALL PHINSINTP(II,JJ,KK,IC,NSPHI) 
          NPHI(II,JJ,KK)=1D-7 ! TO USE THE NPHI AS AN ACTIVATION MARKER
          NPHI(II,JJ,KK)=NSPHI      
         IF (DABS(PSI_CN(II,JJ,KK)-0.5D0) .NE. 0.5D0) THEN 
         PSIH(II,JJ,KK)=DBLE(IC) ! 1 or 0
         ELSE
         PSIH(II,JJ,KK) = PSI_CN(II,JJ,KK)
         ENDIF              
        ENDIF        
       ENDDO
       ENDDO   

! 바깥 경계에 걸리는 Nphi가 없도록 Nphi 에 해당하는 부분을 새롭게 정의하도록 한다.--------------- 
!$OMP PARALLEL DO private(I,J,K,ii,jj,kk,DX,DY,DZ,NPHI_PRE,IC)
        DO NN=1,NUMA(N_MAX,LLVS)
           I =I_B(NN,N_MAX,LLVS) ! ALPHI INDEX
           J =J_B(NN,N_MAX,LLVS)
           K =K_B(NN,N_MAX,LLVS)
           II=ICOUMP_VEL(I)      ! NS INDEX
           JJ=JCOUMP_VEL(J)
           KK=KCOUMP_VEL(K)           
        DX=XPF(I)-XP(II)
        DY=YPF(J)-YP(JJ)
        DZ=ZPF(K)-ZP(KK)                
        NPHI_PRE=ALPHI_EXT(I,J,K,1)+DSQRT(DX**2D0+DY**2D0+DZ**2D0)
        NPHI(II,JJ,KK)=DMIN1(NPHI(II,JJ,KK),NPHI_PRE)
        PSIH(II,JJ,KK) = PSI_CN(II,JJ,KK)                  
        ENDDO
        
! 새롭게 정의한 NPHI (N_MAX)에 대해 REINITIALIZATION을 못하네....마지막 그리드라...
       
!!=====SAVE==============================================================
! TWO_PHASE.F LVSINIT 부분 ALPHI_EXT2 AFTER BANDGENERATION 주석 부분 제거하고
! 여기서 출력되는 것과 어떻게 다른지 비교해보기 (1D CONFIGURATION에서는 LINEAR해서 문제 없음) 
       IF ((MOD(NTIME,10).EQ.0) .and. (msub .eq. 3)) THEN                                                         
      OPEN(142,FILE='0TEST_1.DAT')
       !WRITE(142,*) 'VARIABLES="X","Y","Z","PHI0"'
       WRITE(142,*) 'VARIABLES= "Z","PHI0"'
      !WRITE(142,*)'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
      WRITE(142,*)'ZONE  K=',N3FM,',F=POINT'
      DO K=1,N3FM
      !DO J=1,N2FM
      !DO I=1,N1FM
      !WRITE(142,241)XPF(I),YPF(J),ZPF(K),ALPHI_EXT(I,J,K,1)
      WRITE(142,241) ZPF(K),ALPHI_EXT(1,1,K,1)
      !ENDDO
      !ENDDO
      ENDDO
       CLOSE(142)
! 241  FORMAT(3F16.8,1D20.8)
 241  FORMAT(1F16.8,1D20.8)

       OPEN(148,FILE='0TEST_2.DAT')
       !WRITE(148,*) 'VARIABLES="X","Y","Z","H0","NSPHI"'
       WRITE(148,*) 'VARIABLES= "Z","H0","NSPHI"'
      !WRITE(148,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
      WRITE(148,*)'ZONE  K=',N3M,',F=POINT'
      DO K=1,N3M
      !DO J=1,N2M
      !DO I=1,N1M
!      WRITE(148,147)XP(I),YP(J),ZP(K),PSIH(I,J,K),NPHI(I,J,K)
      WRITE(148,147) zP(K),PSIH(1,1,K),NPHI(1,1,K)
      !ENDDO
      !ENDDO
      ENDDO
       CLOSE(148)
! 147  FORMAT(3F16.8,2D20.8)
 147  FORMAT(1F16.8,2D20.8)
       !STOP
       endif            
!!=======================================================================
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++       
!(2) EXTRAPOLATE THETA_NAVIERSTOKES W.R.T PHASE+++++++++++++++++++++++++
! USE PSIH TO UPDATE ONLY THE TEMP OF THE SAME PHASE.
      CALL TEMP_EXTP(PSI_CN,temp,psih,nphi,T_SAT)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++        
!(3) INTERPOLATE TF USING THE EXTRAPOLATED THETA_NS+++++++++++++++++++++
       DO N=1,NN1 
!$OMP PARALLEL DO private(I,J,K,T_INTP)
       DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         ! OBTAIN TF(GAS) USE TEMP(I,J,K,0)
         CALL TEMPASSIGN(I,J,K,T_INTP,0,TEMP)
         TF(NN,N,LLVS,0)=T_INTP
         check0(i,j,k)=t_intp
         ! OBTAIN TF(LIQ) USE TEMP(I,J,K,1)
         CALL TEMPASSIGN(I,J,K,T_INTP,1,TEMP)
         TF(NN,N,LLVS,1)=T_INTP 
         check1(i,j,k)=t_intp
       ENDDO
       ENDDO  
!!=====SAVE============================================================== 
       IF ((MOD(NTIME,10).EQ.0) .and. (msub .eq. 3)) THEN                                                         
      OPEN(142,FILE='0tf_assign.dat')
       WRITE(142,*) 'VARIABLES="X","Y","Z","tf0","tf1"'
      WRITE(142,*)'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
      DO K=1,N3FM
      DO J=1,N2FM
      DO I=1,N1FM
      WRITE(142,141)XPF(I),YPF(J),ZPF(K),check0(i,j,k),check1(i,j,k)
      ENDDO
      ENDDO
      ENDDO
       CLOSE(142)
 141  FORMAT(3F16.8,2D20.8)       
      endif           
       ! stop
!!=======================================================================         
!(4) OBTAIN DTF/DN(0~1)+++++++++++++++++++++++++++++++++++++++++++++++++            
!(5) EXTRAPOLATE DTF/DN IF NECESSARY++++++++++++++++++++++++++++++++++++            
!(6) OBTAIN MFDOT+++++++++++++++++++++++++++++++++++++++++++++++++++++++  
       CALL TF_GRAD(TF,NGRAD,LLVS,PSIF_MTP,mfdot)  
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++               
!(7) OBTAIN MDOT++++++++++++++++++++++++++++++++++++++++++++++++++++++++              
       IF (LLVS .EQ. 1) THEN       
       CALL CAL_PSI_TMP(1,NN1,LLVS,PSI_TMP,PSIF_MTP,PSI_CN)  
       ALLOCATE(SGAS_X(0:M1,0:M2,0:M3),
     &          SGAS_Y(0:M1,0:M2,0:M3),
     &          SGAS_Z(0:M1,0:M2,0:M3))        
       ELSE       
       CALL CAL_PSI_TMP(0,NN1,LLVS,PSI_TMP,PSIF_MTP,PSI_CN)
       ENDIF

!-----FOR MULTIPLE LVS
!$OMP PARALLEL DO private(I,J,IM1,MDOT)
       DO K=1,N3M
       DO J=1,N2M
       DO I=IBG,N1M
         IM1=IMV(I)
      IF (ABS(PSI_TMP(I,J,K)-PSI_TMP(IM1,J,K)) .NE. 0D0 ) THEN
        CALL CAL_CUR_INTEGRAL(I,J,K,1,MDOT,PSIF_MTP,MFDOT)
        SGAS_X(I,J,K)=SGAS_X(I,J,K)
     &    +(MDOT**2d0)*(PSI_TMP(I,J,K)-PSI_TMP(IM1,J,K))*VVDX(I)
      IMASK_PH(I,J,K)=1
      ENDIF
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J,JM1,MDOT)
       DO K=1,N3M
       DO J=JBG,N2M
         JM1=JMV(J)
       DO I=1,N1M
      IF (ABS(PSI_TMP(I,J,K)-PSI_TMP(I,JM1,K)) .NE. 0D0 ) THEN
        CALL CAL_CUR_INTEGRAL(I,J,K,2,MDOT,PSIF_MTP,MFDOT)
        SGAS_Y(I,J,K)=SGAS_Y(I,J,K)
     &    +(MDOT**2d0)*(PSI_TMP(I,J,K)-PSI_TMP(I,JM1,K))*VVDY(J)
      IMASK_PH(I,J,K)=1
      ENDIF
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J,KM1,MDOT)
       DO K=KBG,N3M
         KM1=KMV(K)
       DO J=1,N2M
       DO I=1,N1M
      IF (ABS(PSI_TMP(I,J,K)-PSI_TMP(I,J,KM1)) .NE. 0D0 ) THEN
        CALL CAL_CUR_INTEGRAL(I,J,K,3,MDOT,PSIF_MTP,sfgas) 
        writE(*,*) '111111',mdot
                if (mdot .lt. 0d0) write(*,*) 'minus mdot?????'
                
        SGAS_Z(I,J,K)=mdot
!        SGAS_Z(I,J,K)
!     &      +mdot*(PSI_TMP(I,J,K)-PSI_TMP(I,J,KM1))*VVDZ(K)
!     &     +(MDOT**2d0)*(PSI_TMP(I,J,K)-PSI_TMP(I,J,KM1))*VVDZ(K)
      IMASK_PH(I,J,K)=1
      ENDIF
       ENDDO
       ENDDO
       ENDDO
              
!      DEALLOCATE(MFDOT)! SFGAS 는 DIVGS 까지 구하고 DEALLOCATE.
      DEALLOCATE(ALPHI_EXT)
      
      
  774 continue

      
      RETURN
      END


!C*******************************************************************      
      SUBROUTINE TEMPASSIGN(I,J,K,T_INTP,IMAT,TEMP)      
! THE DIFFERENT PHASES ARE NOT CONSIDERED IN THE INTERPOLATION (YET)      
!C*******************************************************************      
      USE PARAM_VAR      
      USE FLOW_GEOM_VAR      
      
      USE LVS_GEOM_VAR      
      USE LVS_COUPLING      
    
      IMPLICIT NONE      
      
      REAL*8       TEMP(0:M1,0:M2,0:M3,0:1)             
      INTEGER*8    I,J,K,IMAT       
      REAL*8       T_INTP      
      
      INTEGER*8    II,JJ,KK,II1,JJ1,KK1,II2,JJ2,KK2        
      REAL*8       XD,YD,ZD,C00,C10,C01,C11,C0,C1    
      
      ! IMAT=0 : GAS, IMAT=1 : LIQ
      
         II=ICOUMP_VEL(I) ! INDEX OF THE NS CELL CONTAINING THE LVS CELL
         JJ=JCOUMP_VEL(J)
         KK=KCOUMP_VEL(K)
      !DETERMINE X(NS) INDICES
      IF (XP(II) .GT. XPF(I)) then 
      II1=IMV(II)!짝궁은 II-1
      II2=II
      if (ii1 .eq. ii2) then 
      xd=0d0
      else
      XD=(XPF(I)-(XP(II)-VDX(II)))/(XP(II2)-XP(II1))
      endif
      else ! (XP(II) .LT. XPF(I))
      II1=II
      II2=IPV(II)!짝궁은 II+1
      if (ii1 .eq. ii2) then 
      xd=0d0
      else
       XD=(XPF(I)-XP(II1))/(XP(II2)-XP(II1))
      endif
      endif
      
      !DETERMINE Y(NS) INDICES
      IF (YP(II) .GT. YPF(I)) then 
      JJ1=JMV(JJ)!짝궁은 JJ-1
      JJ2=JJ
      if (jj1 .eq. jj2) then 
      yd=0d0
      else
       YD=(YPF(J)-(YP(JJ)-VDY(JJ)))/(YP(JJ2)-YP(JJ1))      
      endif
      else !IF (YP(II) .LT. YPF(I)) 
      JJ1=JJ  
      JJ2=JPV(JJ)!짝궁은 JJ+1     
      if (jj1 .eq. jj2) then 
      yd=0d0
      else      
       YD=(YPF(J)-YP(JJ1))/(YP(JJ2)-YP(JJ1))      
      endif 
      endif
      
      !DETERMINE Z(NS) INDICES
      IF (ZP(KK) .GT. ZPF(K)) then
      KK1=KMV(KK)!짝궁은 KK-1
      KK2=KK
      if (kk1 .eq. kk2) then 
      zd=0d0
      else
       ZD=(ZPF(K)-(ZP(KK)-VDZ(KK)))/(ZP(KK2)-ZP(KK1))   
      endif 
      else !IF (ZP(KK) .LT. ZPF(K)) 
      KK1=KK
      KK2=KPV(KK)!짝궁은 KK+1    
      if (kk1 .eq. kk2) then 
      zd=0d0
      else
       ZD=(ZPF(K)-ZP(KK1))/(ZP(KK2)-ZP(KK1))   
      endif 
      endif
       
        C00=TEMP(II ,JJ ,KK ,IMAT)*(1-XD)+TEMP(II2,JJ ,KK ,IMAT)*XD
        C10=TEMP(II ,JJ2,KK ,IMAT)*(1-XD)+TEMP(II2,JJ2,KK ,IMAT)*XD
        C01=TEMP(II ,JJ ,KK2,IMAT)*(1-XD)+TEMP(II2,JJ ,KK2,IMAT)*XD
        C11=TEMP(II ,JJ2,KK2,IMAT)*(1-XD)+TEMP(II2,JJ2,KK2,IMAT)*XD 

        C0=C00*(1-YD)+C10*YD
        C1=C01*(1-YD)+C11*YD

        T_INTP=(C0*(1.-ZD)+C1*ZD)
      
      return
      end


!C*******************************************************************      
      SUBROUTINE TEMP_EXTP(PSI_CN,temp,psih,nphi,T_SAT)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR

      IMPLICIT NONE

      REAL*8      NPHI(0:M1,0:M2,0:M3)             
      REAL*8      PSIH(0:M1,0:M2,0:M3) 
      REAL*8      TEMP(0:M1,0:M2,0:M3,0:1)          
      REAL*8      TEMPO(0:M1,0:M2,0:M3)          
      REAL*8     TEMP_RK3(0:M1,0:M2,0:M3) 
      REAL*8     DTEMPN(0:M1,0:M2,0:M3) 
      REAL*8     DTEMPN_RK3(0:M1,0:M2,0:M3) 
      REAL*8     RAT,EPS,SM_GRID,REDTL
      INTEGER*8  ITER_RE,ILIMIT
      INTEGER*8  I,J,K,IM
      REAL*8      DPX,DPY,DPZ,DALPHI,AMDPI,NX,NY,NZ
      REAL*8      DTDX,DTDY,DTDZ
      INTEGER*8   IM1,IP1,JM1,JP1,KM1,KP1 ,ND
      REAL*8      PSI_CN(0:M1,0:M2,0:M3),T_SAT
      
      IF (N3FM .EQ. 1 ) THEN
        SM_GRID=MIN(2*SDXF,2*SDYF)
      ELSE
        SM_GRID=MIN(2*SDXF,2*SDYF,2*SDZF)
      ENDIF
      
        REDTL =0.2D0*SM_GRID   !CFL=<0.5 because of stability of Godunov scheme
        ILIMIT=50              !MAX_ITER=11/CFL (2008 herrmann)
        EPS=1e-4
            DO IM=0,1

        dpx=0d0
        dpy=0d0
        dpz=0d0      
!$OMP PARALLEL DO private(I,J,IM1,IP1,JM1,JP1,KM1,KP1)
!$OMP&private(DPX,DPY,DPZ,DALPHI,NX,NY,NZ,AMDPI)
!$OMP&private(DTDX,DTDY,DTDZ)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
       DTEMPN(I,J,K)=0D0
       DTEMPN_RK3(I,J,K)=0D0
         IM1=IMV(I);JM1=JMV(J);KM1=KMV(K)
         IP1=IPV(I);JP1=JPV(J);KP1=KPV(K)   
!       IF (NPHI(I,J,K) .NE. 0D0) THEN ! McCaslin et al (2014) 참고 
       if (((DABS(nphi(i,j,k)) .GT. 1.5*SM_GRID)
     & .AND. (NPHI(I,J,K).NE. 0D0))
     &      .OR. ((DABS(PSI_CN(I,J,K)-0.5D0) .EQ. 0.5D0)
     &      .AND. (NPHI(I,J,K).NE. 0D0))) THEN 
       
       IF ( NPHI(IM1,J,K) .EQ. 0D0) THEN 
       DPX=(NPHI(IP1,J,K)-NPHI(I,J,K))*VVDX(IP1)
       DTDX=(TEMP(IP1,J,K,IM)-TEMP(I,J,K,IM))*VVDX(IP1)
       ELSEIF (NPHI(IP1,J,K) .EQ. 0D0) THEN 
       DPX=(NPHI(I,J,K)-NPHI(IM1,J,K))*VVDX(I)
       DTDX=(TEMP(I,J,K,IM)-TEMP(IM1,J,K,IM))*VVDX(I)
       ELSE
       DPX=0.5D0*(NPHI(IP1,J,K)-NPHI(IM1,J,K))*VVDX(I)
       DTDX=0.5D0*(TEMP(IP1,J,K,IM)-TEMP(IM1,J,K,IM))*VVDX(I)
       ENDIF
       
       IF ( NPHI(I,JM1,K) .EQ. 0D0) THEN 
       DPY=(NPHI(I,JP1,K)-NPHI(I,J,K))*VVDY(JP1)
       DTDY=(TEMP(I,JP1,K,IM)-TEMP(I,J,K,IM))*VVDY(JP1)
       ELSEIF (NPHI(I,JP1,K) .EQ. 0D0) THEN 
       DPY=(NPHI(I,J,K)-NPHI(I,JM1,K))*VVDY(J)
       DTDY=(TEMP(I,J,K,IM)-TEMP(I,JM1,K,IM))*VVDY(J)
       ELSE
       DPY=0.5D0*(NPHI(I,JP1,K)-NPHI(I,JM1,K))*VVDY(J)
       DTDY=0.5D0*(TEMP(I,JP1,K,IM)-TEMP(I,JM1,K,IM))*VVDY(J)
       ENDIF
       
       IF ( NPHI(I,J,KM1) .EQ. 0D0) THEN 
       DPZ=(NPHI(I,J,KP1)-NPHI(I,J,K))*VVDZ(KP1)
       DTDZ=(TEMP(I,J,KP1,IM)-TEMP(I,J,K,IM))*VVDZ(KP1)
       ELSEIF (NPHI(I,J,KP1) .EQ. 0D0) THEN 
       DPZ=(NPHI(I,J,K)-NPHI(I,J,KM1))*VVDZ(K)
       DTDZ=(TEMP(I,J,K,IM)-TEMP(I,J,KM1,IM))*VVDZ(K)
       ELSE
       DPZ=0.5D0*(NPHI(I,J,KP1)-NPHI(I,J,KM1))*VVDZ(K)
       DTDZ=0.5D0*(TEMP(I,J,KP1,IM)-TEMP(I,J,KM1,IM))*VVDZ(K)
       ENDIF     

       DALPHI=DPX**2D0+DPY**2D0+DPZ**2D0
        if (dalphi .ne. 0d0) then 
        AMDPI=1D0/DSQRT(DALPHI)
        NX=0d0!DPX*AMDPI
        NY=0d0!DPY*AMDPI
        NZ=1d0!DPZ*AMDPI
        
        !if (nz .eq. -1d0) then 
        !write(*,*) 'nz eq -1'
        !write(*,*) nphi(i,j,km1)
        !write(*,*) nphi(i,j,kp1)
        !writE(*,*) nphi(i,j,k  )
        !writE(*,*) dpz, amdpi
        !! stop
        !endif
        
!ONLY THE POINTS WITHIN EACH DOMAIN SHOULD BE USED (MCCASLIN ET AL. 2014)
       IF (IM .EQ. 1) THEN 
       DTEMPN(I,J,K)=NX*DTDX+NY*DTDY+NZ*DTDZ*PSIH(I,J,K)     
       ELSE
       DTEMPN(I,J,K)=NX*DTDX+NY*DTDY+NZ*DTDZ*(1-PSIH(I,J,K))      
       ENDIF      

       ! DTEMPN_RK3(I,J,K)=DTEMPN(I,J,K) 밑에 있음
       ENDIF
       
       !  near the interface
        ELSEIF (((DABS(NPHI(I,J,K)) .LE. 1.5*SM_GRID) .OR.
     &          (DABS(PSI_CN(I,J,K)-0.5D0) .NE. 0.5D0))  .AND.
     &               ((NPHI(I,J,K) .NE. 1D-7 ) .AND.
     &                (NPHI(I,J,K).NE. 0D0)) ) THEN
      IF (IM .EQ. 1) THEN 
      DTEMPN(I,J,K)= (TEMP(I,J,K,IM)-T_SAT)/NPHI(I,J,K)*PSIH(I,J,K)       
      ELSE
      DTEMPN(I,J,K)=(TEMP(I,J,K,IM)-T_SAT)/NPHI(I,J,K)*(1d0-PSIH(I,J,K))   
      ENDIF
       
       ENDIF
       ENDDO
       ENDDO       
       ENDDO  
!!=====SAVE==============================================================    
       IF ((MOD(NTIME,10).EQ.0) .and. (msub .eq. 3)) THEN                                                         
      IF (IM .EQ. 0) THEN 
        OPEN(118,FILE='0DTDN0.DAT') 
      else
        open(118,FILE='0DTDN1.DAT')
      endif
       WRITE(118,*) 'VARIABLES="X","Y","Z","DTMPN"'
       WRITE(118,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      WRITE(118,117)XP(I),YP(J),ZP(K),DTEMPN(I,J,K)
      ENDDO
      ENDDO
      ENDDO
       CLOSE(118)
 117  FORMAT(3F16.8,1D20.8)
      endif           
!!=======================================================================
 998   continue
! for linear extrapolation, extrapolating the gradients should be preceded..    
      DO 999 ITER_RE=1,100       
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
         DTEMPN_RK3(I,J,K)=DTEMPN(I,J,K)
       ENDDO
       ENDDO       
       ENDDO  

        CALL TEMP_EXTP_RK3(DTEMPN,PSIH,NPHI,IM,REDTL,0,DTEMPN) 
        CALL TEMP_EXTP_RK3(DTEMPN,PSIH,NPHI,IM,REDTL,0,DTEMPN)         

!$OMP PARALLEL DO private(I,J)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         DTEMPN(I,J,K)=0.75D0*DTEMPN_RK3(I,J,K)+0.25D0*DTEMPN(I,J,K)
        ENDDO
        ENDDO
        ENDDO
       
       CALL TEMP_EXTP_RK3(DTEMPN,PSIH,NPHI,IM,REDTL,0,DTEMPN)
      
!$OMP PARALLEL DO private(I,J)             
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         DTEMPN(I,J,K)=1d0/3d0*DTEMPN_RK3(I,J,K)+2d0/3d0*DTEMPN(I,J,K)   
        ENDDO
        ENDDO
        ENDDO
  
        !EPS=1D-3
        RAT=0D0
      !STEADY STATE 판별 식. ㅇㅅㅇ A(N+1)/A(N)-1D0< 0.001 이런 느낌?
!$OMP PARALLEL DO private(I,J)
!$OMP&reduction(MAX:RAT)  
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         if (DTEMPN_RK3(i,j,k) .ne. 0d0) then                                                                 
         RAT=DMAX1(DABS(DTEMPN(I,J,K)/DTEMPN_RK3(I,J,K)-1.),RAT)
         elseif (DTEMPN(I,J,K) .ne. 0d0) then 
         RAT=DMAX1(DABS(DTEMPN_RK3(I,J,K)/DTEMPN(I,J,K)-1.),RAT)
         endif
                 
        ENDDO
        ENDDO
        ENDDO      
         !WRITE(*,*) "STEADY?:",RAT!-1D0
      IF (DABS(RAT) .LT. EPS) THEN
      WRITE(*,*) 'CONVERGED',ITER_RE
      GOTO 9991
      ENDIF
 
  999 CONTINUE 
      WRITE(*,*) 'NOT CONVERGED-CHECK DT', rat
      REDTL =0.5D0*redtl
      goto 998
 9991 CONTINUE
  
!!=====SAVE==============================================================  
       IF ((MOD(NTIME,10).EQ.0) .and. (msub .eq. 3)) THEN                                                         
      IF (IM .EQ. 0) THEN 
      OPEN(128,FILE='0DTDN0_3.DAT')
      ELSEIF (IM .EQ. 1) THEN 
      OPEN(128,FILE='0DTDN1_4.DAT')
      ENDIF
      
       WRITE(128,*) 'VARIABLES="X","Y","Z","DTMPN"'
       WRITE(128,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      WRITE(128,127)XP(I),YP(J),ZP(K),DTEMPN(I,J,K)
      ENDDO
      ENDDO
      ENDDO
       CLOSE(128)
 127  FORMAT(3F16.8,1D20.8)   
      endif           
!!=======================================================================      
      
!C==================MAIN REINITIALIZATION CAL._START=====================C
      REDTL =0.2D0*SM_GRID
 1002 continue    
      DO 1000 ITER_RE=1,100
      
!C-----3RD_ORDER TVD RK
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
         TEMP_RK3(I,J,K)=TEMP(I,J,K,im)
         TEMPO(I,J,K)=TEMP(I,J,K,IM)
         ! TEMP_RK3_2(I,J,K)=0d0
       ENDDO
       ENDDO       
       ENDDO   

        CALL TEMP_EXTP_RK3(TEMPO,PSIH,NPHI,IM,REDTL,1,DTEMPN) 
        CALL TEMP_EXTP_RK3(TEMPO,PSIH,NPHI,IM,REDTL,1,DTEMPN)         

!$OMP PARALLEL DO private(I,J)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         TEMPO(I,J,K)=0.75D0*TEMP_RK3(I,J,K)+0.25D0*TEMPO(I,J,K)
        ENDDO
        ENDDO
        ENDDO

       CALL TEMP_EXTP_RK3(TEMPO,PSIH,NPHI,IM,REDTL,1,DTEMPN)
       
!$OMP PARALLEL DO private(I,J)             
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         TEMPO(I,J,K)=1d0/3d0*TEMP_RK3(I,J,K)+2d0/3d0*TEMPO(I,J,K)  
         TEMP(I,J,K,IM)=TEMPO(I,J,K)
        ENDDO
        ENDDO
        ENDDO

        
        !EPS=1D-3
        RAT=0D0
      !STEADY STATE 판별 식. ㅇㅅㅇ A(N+1)/A(N)-1D0< 0.001 이런 느낌?
!$OMP PARALLEL DO private(I,J)
!$OMP&reduction(MAX:RAT)  
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         !TEMP(I,J,K,IM) VS. TEMP_RK3(I,J,K)
         if (temp_rk3(i,j,k) .ne. 0d0) then
         RAT=DMAX1(DABS(TEMP(I,J,K,IM)/TEMP_RK3(I,J,K)-1.),RAT)
        endif
        ENDDO
        ENDDO
        ENDDO      
         !WRITE(*,*) "STEADY?:",RAT!-1D0
      IF (DABS(RAT) .LT. EPS) THEN
      WRITE(*,*) 'CONVERGED',ITER_RE
      GOTO 1001
      ENDIF
 1000 CONTINUE
       WRITE(*,*) 'NOT CONVERGED-CHECK DT', rat
       REDTL =0.5D0*redtl
       goto 1002
 1001 CONTINUE
      
      ENDDO ! do im=0,1
!!=====SAVE==============================================================      
        IF ((MOD(NTIME,10).EQ.0) .and. (msub .eq. 3)) THEN                                                               
!        OPEN(138,FILE='0TEST_1.DAT')
!       WRITE(138,*) 'VARIABLES="X","Y","Z","H0","NSPHI"'
!       WRITE(138,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!      WRITE(138,137)XP(I),YP(J),ZP(K),PSIH(I,J,K),NPHI(I,J,K)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(148)
! 137  FORMAT(3F16.8,2D20.8)      
!
       OPEN(148,FILE='0TEST_temp.DAT')
       WRITE(148,*) 'VARIABLES="X","Y","Z","t0","t1"'
      WRITE(148,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M        
      WRITE(148,147)XP(I),YP(J),ZP(K),TEMP(I,J,K,0),TEMP(I,J,K,1)        
      ENDDO        
      ENDDO        
      ENDDO        
       CLOSE(148)        
 147  FORMAT(3F16.8,2D20.8)             
!      STOP      
      endif           
!!=======================================================================      
      RETURN
      END
      
!C*******************************************************************      
      SUBROUTINE TEMP_EXTP_RK3(TEMPO,PSIH,NPHI,IM,REDTL,ND,DTEMPN)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE   
      INTEGER*8   I,J,K,MM,N,NN,LLVS,NN1,NGRAD,NN2,
     &                               IM1,IM2,IM3,IP1,IP2,
     &                               JM1,JM2,JM3,JP1,JP2,
     &                               KM1,KM2,KM3,KP1,KP2
      REAL*8      TEMP_NEW(0:M1,0:M2,0:M3)
      REAL*8      NPHI(0:M1,0:M2,0:M3)             
      REAL*8      PSIH(0:M1,0:M2,0:M3) 
      REAL*8      TEMPO(0:M1,0:M2,0:M3)  
      REAL*8      DTEMPN(0:M1,0:M2,0:M3)       

      REAL*8      EPS,DPX,DPY,DPZ,DALPHI,AMDPI,ASIGN,NX,NY,NZ
      REAL*8      WENO_X,WENO_Y,WENO_Z,REDTL
      REAL*8      V1,V2,V3,V4,V5,DRHS
      REAL*8      APXM,APXP,ALPHI_X,
     &            APYM,APYP,ALPHI_Y,
     &            APZM,APZP,ALPHI_Z
      integer*8   im,IMAT,ND
      INTEGER*8   COLOR(MF_BAND,M_MAX)
      real*8      nxo,nyo,nzo
      real*8 sm_grid
        SM_GRID=MIN(2*SDXF,2*SDYF,2*SDZF) 
      drhs=0d0
        IF (N3F .EQ. 2 ) THEN
         EPS=MIN(2*SDXF,2*SDYF)**2
        ELSE
         EPS=MIN(2*SDXF,2*SDYF,2*SDZF)**2
        ENDIF
       

        dpx=0d0        ;ALPHI_X=0d0
        dpy=0d0        ;ALPHI_Y=0d0
        dpz=0d0        ;ALPHI_Z=0d0
! 일반 그리드 사용하는 경우에는 LIB가 아니라서 그냥 모든 점에 대해 구해주면 된다...
!$OMP PARALLEL DO private(I,J,IM1,JM1,KM1) 
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M            
        TEMP_NEW(I,J,K)=tempO(i,j,k)
       ENDDO
       ENDDO
       ENDDO

!$OMP PARALLEL DO private(I,J,IM1,IP1,JM1,JP1,KM1,KP1)
!$OMP&private(DPX,DPY,DPZ,DALPHI,ASIGN,DRHS,NX,NY,NZ,AMDPI)
!$OMP&private(APXM,APXP,ALPHI_X,APYM,APYP,ALPHI_Y)
!$OMP&private(APZM,APZP,ALPHI_Z)
!$OMP&private(NXo,NYo,NZo)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         IM1=IMV(I);JM1=JMV(J);KM1=KMV(K)
         IP1=IPV(I);JP1=JPV(J);KP1=KPV(K)
         
         if (nphi(i,j,k) .ne. 0d0) then     
       IF ( NPHI(IM1,J,K) .EQ. 0D0) THEN 
       DPX=(NPHI(IP1,J,K)-NPHI(I,J,K))*VVDX(IP1)
       ELSEIF (NPHI(IP1,J,K) .EQ. 0D0) THEN 
       DPX=(NPHI(I,J,K)-NPHI(IM1,J,K))*VVDX(I)   
       ELSE
       DPX=0.5D0*(NPHI(IP1,J,K)-NPHI(IM1,J,K))*VVDX(I)
       ENDIF
       
       IF ( NPHI(I,JM1,K) .EQ. 0D0) THEN 
       DPY=(NPHI(I,JP1,K)-NPHI(I,J,K))*VVDY(JP1)
       ELSEIF (NPHI(I,JP1,K) .EQ. 0D0) THEN 
       DPY=(NPHI(I,J,K)-NPHI(I,JM1,K))*VVDY(J)
       ELSE
       DPY=0.5D0*(NPHI(I,JP1,K)-NPHI(I,JM1,K))*VVDY(J)    
       ENDIF
       
       IF ( NPHI(I,J,KM1) .EQ. 0D0) THEN 
       DPZ=(NPHI(I,J,KP1)-NPHI(I,J,K))*VVDZ(KP1)   
       ELSEIF (NPHI(I,J,KP1) .EQ. 0D0) THEN 
       DPZ=(NPHI(I,J,K)-NPHI(I,J,KM1))*VVDZ(K)     
       ELSE
       DPZ=0.5D0*(NPHI(I,J,KP1)-NPHI(I,J,KM1))*VVDZ(K)  
       ENDIF        

       DALPHI=DPX**2D0+DPY**2D0+DPZ**2D0
       if (dalphi .ne. 0d0) then 
       AMDPI=1D0/DSQRT(DALPHI)
        NX=0d0! DPX*AMDPI
        NY=0d0! DPY*AMDPI
        NZ=1d0! DPZ*AMDPI
       
        ASIGN=NPHI(I,J,K)/dSQRT(NPHI(I,J,K)**2d0)!+DALPHI*EPS)        
        APXM=(TEMPO(I  ,J  ,K  )-TEMPO(IM1,J  ,K  ))*VVDX(I)     
        APXP=(TEMPO(IP1,J  ,K  )-TEMPO(I  ,J  ,K  ))*VVDX(I)     
        APYM=(TEMPO(I  ,J  ,K  )-TEMPO(I  ,JM1,K  ))*VVDY(J)     
        APYP=(TEMPO(I  ,JP1,K  )-TEMPO(I  ,J  ,K  ))*VVDY(J)     
        APZM=(TEMPO(I  ,J  ,K  )-TEMPO(I  ,J  ,KM1))*VVDZ(K)     
        APZP=(TEMPO(I  ,J  ,KP1)-TEMPO(I  ,J  ,K  ))*VVDZ(K)     
        

        IF (IM .EQ. 1) THEN ! 반대방향으로 채워야해서 그래... 
        NX=(-1D0)*NX
        NY=(-1D0)*NY
        NZ=(-1D0)*NZ
        ENDIF       
        ALPHI_X=DMIN1(NX,0D0)*APXP+DMAX1(NX,0D0)*APXM
        ALPHI_Y=DMIN1(NY,0D0)*APYP+DMAX1(NY,0D0)*APYM
        ALPHI_Z=DMIN1(NZ,0D0)*APZP+DMAX1(NZ,0D0)*APZM  

        
       
        
!C***** EULER METHOD
      IF (IM .EQ. 1) THEN  ! extrapolating the liquid T to gas region
       DRHS=(ALPHI_X +ALPHI_Y +ALPHI_Z)+ND*DTEMPN(I,J,K)
       TEMP_NEW(I,J,K)=TEMPO(I,J,K)-REDTL*DRHS*(1-PSIH(I,J,K))       
      ELSE                 ! extrapolating the gas T to liquid region
       DRHS=(ALPHI_X +ALPHI_Y +ALPHI_Z)-ND*DTEMPN(I,J,K)
       TEMP_NEW(I,J,K)=TEMPO(I,J,K)-REDTL*DRHS*PSIH(I,J,K)      
      ENDIF

! linear extrapolation is done with less grids: empty slots at the boundary 
      IF ((ND .EQ. 0).AND.(DABS(NPHI(I,J,K)) .LE. SM_GRID)) THEN 
         DRHS= (ALPHI_X +ALPHI_Y +ALPHI_Z)
         TEMP_NEW(I,J,K)=TEMPO(I,J,K)-REDTL*DRHS 
      ENDIF

       endif      !if (dalphi .ne. 0d0) then 
       endif      !if (nphi(i,j,k) .ne. 0d0) then    

       ENDDO
       ENDDO
       ENDDO
              
!C=================== ONLY DONE NEAR THE INTERFACE=======================C
!$OMP PARALLEL DO private(I,J)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
          TEMPO(I,J,K)=TEMP_NEW(I,J,K)
        if (nd .eq. 0) DTEMPN(I,J,K)=tempo(i,j,k)
        ENDDO
        ENDDO
       ENDDO

      RETURN
      END 

!C*******************************************************************      
      SUBROUTINE TF_GRAD(TF,NGRAD,LLVS,PSIF_MTP,mfdot)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      USE TWO_PHASE_PROPERTY
      USE PHASE_BOILING
      USE HEAT_VAR
      
      IMPLICIT NONE
      INTEGER*8   I,J,K,N,NN
      INTEGER*8   NGRAD,LLVS,NN1,NN2
      ! REAL*8      ALPHI_TMP(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      REAL*8      DTDN(-2:M1F+3,-2:M2F+3,-2:M3F+3,0:1)
      REAL*8      TEMP_EXT(-2:M1F+3,-2:M2F+3,-2:M3F+3,0:1)
      INTEGER*8   MASK(-2:M1F+3,-2:M2F+3,-2:M3F+3) 
      
      INTEGER*8   iMASK(-2:M1F+3,-2:M2F+3,-2:M3F+3) 
      REAL*8      TF(MF_BAND,M_MAX,MLVS,0:1)      
      
       REAL*8      MFDOT(M1L:M1U,M2L:M2U,M3L:M3U)  
       REAL*8      PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U) 
       REAL*8      GPSI
       REAL*8      GX,GY,GZ,AGI,NORM_X,NORM_Y,NORM_Z
       REAL*8      DTX0,DTY0,DTZ0,DTX1,DTY1,DTZ1
       INTEGER*8   IM1,IP1,JM1,JP1,KM1,KP1
       
       
      REAL*8      CURF_TMP(M1L:M1U,M2L:M2U,M3L:M3U)
      REAL*8      GXX,GYY,GZZ,GXY,GXZ,GYZ,dvol
      real*8      place

!       ALLOCATE( ALPHI_EXT(-2:M1F+3,-2:M2F+3,-2:M3F+3,1))
!       ALLOCATE(MFDOT(M1L:M1U,M2L:M2U,M3L:M3U))
!       ALLOCATE(SFGAS(M1L:M1U,M2L:M2U,M3L:M3U))
       
      dvol=sdxf*sdyf*sdzf 
       
      nn1=n_max-3
      nn2=n_max-1

      IF (LLVS .EQ. 1) THEN
!$OMP PARALLEL DO            
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       DTDN(I,J,K,0)     =0D0
       DTDN(I,J,K,1)     =0d0
       MASK(I,J,K)       =0
       iMASK(I,J,K)      =0
       TEMP_EXT(I,J,K,0) =0d0
       TEMP_EXT(I,J,K,1) =0d0
       ! tf(i,j,k)  =0d0
      ENDDO
      ENDDO
      ENDDO      
      ENDIF
      
	  !EMPTY SFGAS(I,J,K)
!$OMP PARALLEL DO private(I,J,K)	
        DO K=1,N3FM
        DO J=1,N2FM
        DO I=1,N1FM
         SFGAS(I,J,K)=0D0
         SFGAS2(I,J,K)=0D0
         MFDOT(I,J,K)=0D0
        ENDDO
        ENDDO
        ENDDO
        
      ! SAVE EACH LEVEL SET FUNCT TO ALPHI_TMP      
      ! FOR EACH LEVEL-SET, SAVE TO ALPHI_TMP(I,J,K)
       DO N=1,NN1
!$OMP PARALLEL DO private(I,J,K)
       DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS) 
         TEMP_EXT(I,J,K,0)=TF(NN,N,LLVS,0)         
         TEMP_EXT(I,J,K,1)=TF(NN,N,LLVS,1)         
       ENDDO
       ENDDO  

       IF (LLVS .EQ. 1) THEN
       CALL CAL_PSIF_MTP(1,NN1,LLVS,PSIF_MTP)     
       ELSE
       CALL CAL_PSIF_MTP(0,NN1,LLVS,PSIF_MTP)
       ENDIF
  
!================CALCULATE THE TEMPERATURE GRADIENT (REFINED)===========
      DO N=1,Nn1-2
!$OMP PARALLEL DO private(I,J,K,GX,GY,GZ,AGI)
!$OMP&private(DTX0,DTY0,DTZ0,DTX1,DTY1,DTZ1)
!$OMP&private(IM1,IP1,JM1,JP1,Km1,KP1)
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
       
      GX=0.5D0*( ALPHI_EXT(ip1,J,K,1)-ALPHI_EXT(im1,J,K,1) )*SSDXF
      GY=0.5D0*( ALPHI_EXT(I,jp1,K,1)-ALPHI_EXT(I,jm1,K,1) )*SSDYF
      GZ=0.5D0*( ALPHI_EXT(I,J,kp1,1)-ALPHI_EXT(I,J,km1,1) )*SSDZF     
      AGI=1D0/DSQRT(GX**2D0 + GY**2D0 + GZ**2D0)          
      ! if (n .lt. nn1) then 
      DTX0=0.5D0*(TEMP_EXT(IP1,J,K,0)-TEMP_EXT(IM1,J,K,0))*SSDXF
      DTY0=0.5D0*(TEMP_EXT(I,JP1,K,0)-TEMP_EXT(I,JM1,K,0))*SSDYF
      DTZ0=0.5D0*(TEMP_EXT(I,J,KP1,0)-TEMP_EXT(I,J,KM1,0))*SSDZF
      DTDN(I,J,K,0)=(DTX0*GX+DTY0*GY+DTZ0*GZ)*AGI
      
      DTX1=0.5D0*(TEMP_EXT(IP1,J,K,1)-TEMP_EXT(IM1,J,K,1))*SSDXF
      DTY1=0.5D0*(TEMP_EXT(I,JP1,K,1)-TEMP_EXT(I,JM1,K,1))*SSDYF
      DTZ1=0.5D0*(TEMP_EXT(I,J,KP1,1)-TEMP_EXT(I,J,KM1,1))*SSDZF
      DTDN(I,J,K,1)=(DTX1*GX+DTY1*GY+DTZ1*GZ)*AGI
      
      ! 왜 이렇게 하면 뻑나는지 모르겟다. (oscillate) ! alphi가 작아질 때 불안정해지는 것 같음
!      IF ( (DABS(PSIF_MTP(I,J,K)-0.5D0) .NE. 0.5D0) .and. 
!     &     (dabs(alphi_ext(i,j,k,1)) .gt. 1d-3))THEN 
!      DTDN(I,J,K,1)=(TEMP_EXT(I,J,K,1)-T_SAT)/ALPHI_EXT(I,J,K,1)
!      DTDN(I,J,K,0)=(TEMP_EXT(I,J,K,0)-T_SAT)/ALPHI_EXT(I,J,K,1)
!      endif     
      ENDDO
      ENDDO
      
! 이게 먼저 다 끝나고 그다음에 계산해야 차분이 가능..  
!=====SAVE==============================================================
       IF ((MOD(NTIME,10).EQ.0) .and. (msub .eq. 3)) THEN                                                         
       OPEN(148,FILE='0dtdn_pre.DAT')
       WRITE(148,*) 'VARIABLES="X","Y","Z","dtdn0","dtdn1"'
      WRITE(148,*)'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
      DO K=1,N3FM
      DO J=1,N2FM
      DO I=1,N1FM
      WRITE(148,147)XPF(I),YPF(J),ZPF(K),dtdn(i,j,k,0),dtdn(i,j,k,1)
      ENDDO
      ENDDO
      ENDDO
       CLOSE(148)
!  147  FORMAT(3F16.8,2D20.8)
       !stop
        endif           
!=======================================================================
!      CALL DTDNBC(DTDN) ! 이게 필요한가..?

      ! M_DOT_REFINED를 구한다.      
      DO N=1,nn1-2
!!$OMP PARALLEL DO private(I,J,K)
!!$OMP&private(IM1,IP1,JM1,JP1,Km1,KP1)
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
!"For efficient implementation of U=uf+mdot*n/rho,            
! the mass flux, which is evaluated at the interface, is extrapolated into the entire domain."
! Son & Dhir, 2007 (여기서는 lvs band 까지만 구한다.)            
      MFDOT(I,J,K)=JA_AIR*PRAI*REI!*TCM*dtdn(i,j,k,0)
     &                   *(-TCP*DTDN(I,J,K,1)+TCM*DTDN(I,J,K,0)) ! GRAD(T)*NORM 
                          !!the signs are changed as the temp is normalized in such way
      MFLVS(NN,N,LLVS)=MFDOT(I,J,K)
        ENDDO
      ENDDO            
 
!(*) ADJUST VOL_TOT_ORI WITH MFDOT!+++++++++++++++++++++++++++++++++++++
       CALL VOL_TOT_MFDOT(PSIF_MTP,LLVS)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ 
!!=====SAVE=============================================================
        IF ((MOD(NTIME,10).EQ.0) .and. (msub .eq. 3)) THEN                                                        
       OPEN(146,FILE='0mfdots.DAT')
       WRITE(146,*) 'VARIABLES="X","Y","Z","mfdots","sfgas"'
      WRITE(146,*)'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
      DO K=1,N3FM
      DO J=1,N2FM
      DO I=1,N1FM
      WRITE(146,147)XPF(I),YPF(J),ZPF(K),mfdot(i,j,k),sfgas(i,j,k)
      ENDDO
      ENDDO
      ENDDO
       CLOSE(146)
 147  FORMAT(3F16.8,2D20.8)
       endif     
      !if (m .eq. 1) stop
      !if (m .eq. 10) stop
!======================================================================= 

      RETURN
      END

!C*******************************************************************      
      SUBROUTINE DTDNBC(DTDN)
!C*******************************************************************      
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE
      REAL*8      DTDN(-2:M1F+3,-2:M2F+3,-2:M3F+3,0:1) 
      INTEGER*8   I,J,K


      !BOUNDARY CONDITION
      IF ( IPX .EQ. 1 ) THEN
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F    !START FROM ZERO BECAUSE OF EDGE OF DOMAIN.(NEED FOR INTER_PSI)
      DO J=0,N2F
      DO I=0,N1F_BD
        DTDN(-I,J,K,0)=DTDN(N1FM-I,J,K,0)
        DTDN(-I,J,K,1)=DTDN(N1FM-I,J,K,1)
        DTDN(N1F+I,J,K,0)=DTDN(1+I,J,K,0)
        DTDN(N1F+I,J,K,1)=DTDN(1+I,J,K,1)
      ENDDO
      ENDDO
      ENDDO
      ELSE
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F    !START FROM ZERO BECAUSE OF EDGE OF DOMAIN.(NEED FOR INTER_PSI)
      DO J=0,N2F
      DO I=0,N1F_BD
        DTDN(-I,J,K,0)=DTDN(1,J,K,0)
        DTDN(N1F+I,J,K,0)=DTDN(N1FM,J,K,0)
        DTDN(-I,J,K,1)=DTDN(1,J,K,1)
        DTDN(N1F+I,J,K,1)=DTDN(N1FM,J,K,1)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      IF ( IPY .EQ. 1 ) THEN
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F
      DO J=0,N2F_BD
      DO I=0,N1F
        DTDN(I,-J,K,0)=DTDN(I,N2FM-J,K,0)
        DTDN(I,N2F+J,K,0)=DTDN(I,1+J,K,0)
        DTDN(I,-J,K,1)=DTDN(I,N2FM-J,K,1)
        DTDN(I,N2F+J,K,1)=DTDN(I,1+J,K,1)
      ENDDO
      ENDDO
      ENDDO
      ELSE
!$OMP PARALLEL DO private(I,J)
      DO K=0,N3F
      DO J=0,N2F_BD
      DO I=0,N1F
        DTDN(I,-J,K,0)=DTDN(I,1,K,0)
        DTDN(I,N2F+J,K,0)=DTDN(I,N2FM,K,0)
        DTDN(I,-J,K,1)=DTDN(I,1,K,1)
        DTDN(I,N2F+J,K,1)=DTDN(I,N2FM,K,1)
      ENDDO
      ENDDO
      ENDDO
      ENDIF
      IF ( IPZ .EQ. 1 ) THEN
      DO K=0,N3F_BD
!$OMP PARALLEL DO private(I)
      DO J=0,N2F
      DO I=0,N1F
        DTDN(I,J,-K,0)=DTDN(I,J,N3FM-K,0)
        DTDN(I,J,N3F+K,0)=DTDN(I,J,1+K,0)
        DTDN(I,J,-K,1)=DTDN(I,J,N3FM-K,1)
        DTDN(I,J,N3F+K,1)=DTDN(I,J,1+K,1)
      ENDDO
      ENDDO
      ENDDO
      ELSE
      DO K=0,N3F_BD
!$OMP PARALLEL DO private(I)
      DO J=0,N2F
      DO I=0,N1F
        DTDN(I,J,-K,0)=DTDN(I,J,1,0)
        DTDN(I,J,N3F+K,0)=DTDN(I,J,N3FM,0)
        DTDN(I,J,-K,1)=DTDN(I,J,1,1)
        DTDN(I,J,N3F+K,1)=DTDN(I,J,N3FM,1)
      ENDDO
      ENDDO
      ENDDO
      ENDIF


      RETURN
      END
  
!***********************************************************************
      SUBROUTINE PHINSINTP(IO,JO,KO,IC,NSPHI)
!***********************************************************************     
! IO,JO,KO ARE NS INDEX
      USE PARAM_VAR      
      USE FLOW_GEOM_VAR            
      USE LVS_GEOM_VAR 
      USE LVS_VAR
      USE LVS_COUPLING                  

      implicit none
      REAL*8      NSPHI
      INTEGER*8   IC
      integer*8   I1,I2,J1,J2,K1,K2
      INTEGER*8   IO,JO,KO
      REAL*8      C00,C10,C01,C11,C0,C1,XD,YD,ZD
      integer*8   i
      
      I1=ICOU2(IO)
      I2=ICOUMP2(IO)
      J1=JCOU2(JO)
      J2=JCOUMP2(JO)
      K1=KCOU2(KO)
      K2=KCOUMP2(KO)
      
      if (i1 .eq. i2) then 
      xd=0d0
      else
      XD=(XP(IO)-XPF(I1))/(XPF(I2)-XPF(I1))
      endif
      
      if (j1 .eq. j2) then 
      yd=0d0
      else
       YD=(YP(JO)-YPF(J1))/(YPF(J2)-yPF(J1))
      endif
      
      if (k1 .eq. k2) then
      zd=0d0
      else
      ZD=(ZP(KO)-ZPF(K1))/(ZPF(K2)-zPF(K1))
      endif
        C00=ALPHI_EXT(I1,J1,K1,1)*(1d0-XD)+ALPHI_EXT(I2,J1,K1,1)*XD
        C10=ALPHI_EXT(I1,J2,K1,1)*(1d0-XD)+ALPHI_EXT(I2,J2,K1,1)*XD
        C01=ALPHI_EXT(I1,J1,K2,1)*(1d0-XD)+ALPHI_EXT(I2,J1,K2,1)*XD
        C11=ALPHI_EXT(I1,J2,K2,1)*(1d0-XD)+ALPHI_EXT(I2,J2,K2,1)*XD 

        C0=C00*(1d0-YD)+C10*YD
        C1=C01*(1d0-YD)+C11*YD

        NSPHI=(C0*(1d0-ZD)+C1*ZD)
       
        IF (NSPHI .GT. 0D0) IC=1 
        IF (NSPHI .LT. 0D0) IC=0    
      !?? nn2 brink에 있는 녀석들 중 alphi 값이 없는 cell을 참조할 때가 있어서,,,
        IF (ALPHI_EXT(I1,J1,K1,1) .EQ. 0D0 .OR.
     &      ALPHI_EXT(I1,J2,K1,1) .EQ. 0D0 .OR.
     &      ALPHI_EXT(I1,J1,K2,1) .EQ. 0D0 .OR.
     &      ALPHI_EXT(I1,J2,K2,1) .EQ. 0D0 .OR.
     &      ALPHI_EXT(I2,J1,K1,1) .EQ. 0D0 .OR.
     &      ALPHI_EXT(I2,J2,K1,1) .EQ. 0D0 .OR.
     &      ALPHI_EXT(I2,J1,K2,1) .EQ. 0D0 .OR.
     &      ALPHI_EXT(I2,J2,K2,1) .EQ. 0D0 ) THEN 
       NSPHI =0D0 ! 괜한 어중간한 값대신 
       ! WRITE(*,*) XPF(I1),YPF(J1),ZPF(K1)
        ENDIF
        
      return
      end

!***********************************************************************
      SUBROUTINE VOL_TOT_MFDOT(PSIF_MTP,LLVS)
!***********************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      USE TWO_PHASE_PROPERTY
      USE PHASE_BOILING
      USE HEAT_VAR
      
      IMPLICIT NONE
      INTEGER*8    IM1,IP1,JM1,JP1,KM1,KP1       
      INTEGER*8    NN1,NN2,I,J,K,N,NN,LLVS,ISUBCYCLE,IMASS
      
      REAL*8      PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U) 
      REAL*8      GPSI
      REAL*8      GX,GY,GZ,AG

      REAL*8       DTLVS,TIME_TRAN,CFLM,CFLM_OMP,CFLX,CFLY,CFLZ,CFLI
      REAL*8       UUVVWW
      REAL*8       DALPHIMAX,DALPHIMIN      

      REAL*8      ALPHI_S(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      REAL*8      ALPHI_RK3(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      REAL*8      PSIFMTP_S(M1L:M1U,M2L:M2U,M3L:M3U) 
      REAL*8      PSIF_MTP2(M1L:M1U,M2L:M2U,M3L:M3U) 
       REAL*8     MFD(-2:M1F+3,-2:M2F+3,-2:M3F+3) ,DVOL,VOL_TOT
       real*8     place,place2
       DVOL=SDXF*SDYF*SDZF
       
      NN1=N_MAX-3
      NN2=N_MAX-1
      
!$OMP PARALLEL DO private(I,J)         
      DO K=0,N3F
      DO J=0,N2F
      DO I=0,N1F
      ALPHI_S(I,J,K)=ALPHI_EXT(I,J,K,1)
      PSIFMTP_S(I,J,K)=PSIF_MTP(I,J,K)
      ENDDO
      ENDDO
      ENDDO

!C-----CAL. TRANSPORT DT
      DO N=1,NN1-2 
!$OMP PARALLEL DO private(I,J,K)
!$OMP&private(IP1,IM1,JP1,JM1,KP1,KM1)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         IM1=IMF(I) ;JM1=JMF(J) ;  KM1=KMF(K)
         IP1=IPF(I) ;JP1=JPF(J) ;  KP1=KPF(K)             

        IF (ALPHI_EXT(I,J,K,1) .GT. 0D0) THEN
        MFD(I,J,K) = MFLVS(NN,N,LLVS)/DENR!*AG
        ELSE                              
        MFD(I,J,K) = MFLVS(NN,N,LLVS)!/denr     !*AG
        ENDIF
   
        ENDDO
       ENDDO       


!provisional 하게 vol_tot을 구하려는 step이라서 굳이 dtlvs으로 쪼개지 않음. 
!추후 사용될 dtconst를 바로 적용. 
!C-----3RD_ORDER RK-----------------------------------------------------
       DO N=1,NN1-2                                                       
!$OMP PARALLEL DO PRIVATE(I,J,K)                                        
        DO NN=1,NUMA(N,LLVS)                                            
         I=I_B(NN,N,LLVS)                                               
         J=J_B(NN,N,LLVS)                                               
         K=K_B(NN,N,LLVS)
         
      ALPHI_RK3(I,J,K)=ALPHI_EXT(I,J,K,1)                                                
      ALPHI_EXT(I,J,K,1)=ALPHI_EXT(I,J,K,1)-2*DTconst*MFD(I,J,K)     
      ALPHI_EXT(I,J,K,1)=0.75*ALPHI_RK3(I,J,K)+0.25*ALPHI_EXT(I,J,K,1)                   
      ALPHI_EXT(I,J,K,1)=ALPHI_EXT(I,J,K,1)-DTconst*MFD(I,J,K)  
      ALPHI_EXT(I,J,K,1)=1./3.*ALPHI_RK3(I,J,K)+2./3.*ALPHI_EXT(I,J,K,1)        
      !ALPHI_EXT(I,J,K,1)=ALPHI_EXT(I,J,K,1)-dtconst*mfd(i,j,k)     
      
        ENDDO                                                           
       ENDDO                                                                                                                        
!C-----3RD_ORDER RK-----------------------------------------------------
       call reinitialization(nn1,nn2,nn2,llvs) !?해줘야하나? 해줘야하네..(해줘야 전후 vol이 비슷해짐)

       
!!=====SAVE=============================================================  
      place=0d0
       do k=1,n3fm
       if (ALPHI_EXT(2,2,k,1)*ALPHI_EXT(2,2,k+1,1) .lt. 0d0) then
       place=zpf(k) + DABS(ALPHI_EXT(2,2,K,1))
       place2=mfd(2,2,k)
!       sdzf *dabs(alphi_ext(2,2,k,1))
!     & /(dabs(alphi_ext(2,2,k,1))+dabs(alphi_ext(2,2,k+1,1)))     
       write(*,*) 'interface',place
       endif
       enddo      
       
       if (msub .eq. 3) then
      OPEN(100,FILE='0interface.dat',POSITION='APPEND')
      OPEN(102,FILE='0interface2.dat',POSITION='APPEND')
      WRITE(100,101)TIME,place
      WRITE(102,103)TIME,place,place2
      CLOSE(100)
      CLOSE(102)
 101   FORMAT(2ES15.7) 
 103   FORMAT(3ES15.7) 
       endif
!=======================================================================

! OBTAIN NEW PSIF_MTP WITH THE NEW ALPHI_EXT, INOUT=0; NOT BIG CHANGE
      CALL CAL_PSIF_MTP(1,NN1,LLVS,PSIF_MTP2)   
        VOL_TOT=0D0
!$OMP PARALLEL DO PRIVATE(I,J)
!$OMP&REDUCTION(+:VOL_TOT)
        DO K=1,N3FM
        DO J=1,N2FM
        DO I=1,N1FM
         VOL_TOT=VOL_TOT+(1D0-PSIF_MTP2(I,J,K))*DVOL
        ENDDO
        ENDDO
        ENDDO       

!여기까지 VOL_TOT_ORI(LLVS)를 끌고 다닐 수 없으니까, 
!PHASE_CHANGE에 종속되는 변수를 정의해서, IPHS가 1이면 적용해주는 걸로 하자.        
      VOL_TOT_IPHS(LLVS)=VOL_TOT
      
      !SFGAS IS DIFFERENCE BETWEEN PSIF_MTP AND PSIF_MTP2 * DVOL
      DO N=1,NN1-2
!$OMP PARALLEL DO PRIVATE(I,J,K)                                        
      DO NN=1,NUMA(N,LLVS)   
       IF (mflvs(nn,n,llvs) .NE. 0D0) THEN 
         I=I_B(NN,N,LLVS)                                               
         J=J_B(NN,N,LLVS)                                               
         K=K_B(NN,N,LLVS)       
      IF  (DABS(PSIF_MTP(I,J,K) -0.5D0) .NE. 0.5D0)  THEN   
            if (alphi(nn,n,llvs) .gt. 0d0) then 
      SFGAS(I,J,K)=-(PSIF_MTP2(I,J,K)-PSIF_MTP(I,J,K))*denr !*DVOL ! increased gas volume fraction
      SFGAS2(I,J,K)=-(PSIF_MTP2(I,J,K)-PSIF_MTP(I,J,K))!*DVOL ! increased gas volume fraction
      else
      sfgas(i,j,k)=-(PSIF_MTP2(I,J,K)-PSIF_MTP(I,J,K))
      sfgas2(i,j,k)=-(PSIF_MTP2(I,J,K)-PSIF_MTP(I,J,K))
      endif
      ! SGLVS(NN,N,LLVS)=SFGAS(I,J,K)*DVOL
      ENDIF
      ENDIF
      ENDDO
      ENDDO
 
! restore the original values.  
!$OMP PARALLEL DO private(I,J)   
      DO K=0,N3F
      DO J=0,N2F
      DO I=0,N1F
      ALPHI_EXT(I,J,K,1)=ALPHI_S(I,J,K)
      PSIF_MTP(I,J,K)=PSIFMTP_S(I,J,K)
      ENDDO
      ENDDO
      ENDDO
      
      RETURN
      END
      
!***********************************************************************
      SUBROUTINE TEST_TFIELD(XX,YY,ZZ,sigma1,TEMP)
!***********************************************************************     
! to make the initial/analytic temperature field of 
! growing bubble in superheated liquid, without gravity.

      implicit none
      
      integer*8   N,i,M,k
      real*8      x,sigma1,sigma2,AA,rad,R0
      REAL*8      TEMP,XX,YY,ZZ

      N=50000

      R0=1d0 ! NONDIME WITH R0

      rad=DSQRT(XX**2+YY**2+ZZ**2) !R0+R0*k*0.0005d0
      AA=R0/rad !1d0-R0/rad
      sigma2=0d0
    
      do i=1,N-2
      x=(1d0-R0/rad)+i*AA/N
      sigma2=sigma2+(AA/N)*exp(
     & -((15.2465d0)**2d0)*((1d0-x)**(-2d0)-2d0*0.9994d0*x-1d0))
      enddo

      TEMP=sigma2/sigma1 

      return
      end




      
      


      
      