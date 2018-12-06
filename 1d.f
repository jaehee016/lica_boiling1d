!***********************************************************************
!***********************************************************************
!***********************************************************************

      PROGRAM MAIN
      USE PARAM_VAR     
      USE FLOW_VAR
      USE FLOW_GEOM_VAR
      USE FLD_AVG      
      USE MG_OPTION
      USE HEAT_VAR
   
      
      USE TWO_PHASE_PROPERTY
      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING
      
      USE TIME_VAR
      USE PHASE_BOILING
      
      IMPLICIT NONE
      
      CHARACTER*10 DUMMY
      CHARACTER*30 gridfile
      CHARACTER*30 fileprevel

      INTEGER*8   I,J,K,L,N,NN,LLVS,IP1,JP1                                                     
      INTEGER*8   ICFL_CHANGE
      REAL*8      FUNCBODY,PI
      INTEGER*8   IHIST      

      REAL*8      DENP_OLD,DENM_OLD,VISP_OLD,VISM_OLD,SURF_OLD,FR_OLD
      REAL*8      DENP_NEW,DENM_NEW,VISP_NEW,VISM_NEW,SURF_NEW,FR_NEW
      REAL*8      TIME_OLD,TIME_INTERVAL,T_TAR
      INTEGER*8   ITIME_VARING_PROPERTIES,NV,NAV,L1,L2
      
      INTEGER*8   ITER_NS
      REAL*8      DEN_DIFFI,U_MAX,U_MAX_OLD,E_SUM,VOL
      REAL*8      CFLM

      REAL*8      U_BUB_AVG,VOL_BUB,DVOL,DEN_CN
      REAL*8      UU,VV,WW,VEL_KIM,E_AVG,DVMAX,QMMAX
      REAL*8      E_SUM_TMP,U_BUB_AVG_TMP,U_MAX_TMP,VOL_TMP,VOL_BUB_TMP
      REAL*8      FTRTIME1,FTRTIME2,FTRTIME3,FTRTIME4,FTRTIMPS,FTRTIME5
      INTEGER*8   IERR
      INTEGER*8   NNVS,NTII
      REAL*8       PHI_MAX,PHI_MIN,PHI_MAX_TMP,PHI_MIN_TMP,T_PERIOD             

      !FOR INITIAL MEANPGR
      REAL*8, DIMENSION (:,:,:), ALLOCATABLE :: DENF_Z
      !RK3_OLD
      REAL*8, ALLOCATABLE :: RK3XO(:,:,:),RK3YO(:,:,:),RK3ZO(:,:,:)       
      !RK3_OLD_ENERGY
      REAL*8, ALLOCATABLE :: ANPSO(:,:,:) 
      !TRACE
      INTEGER*8 ITR(10000),JTR(10000),KTR(10000)
      INTEGER*8 ITR2(10000),JTR2(10000),KTR2(10000)      
            
      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 P(0:M1,0:M2,0:M3)
      REAL*8 PSI_XN(0:M1,0:M2,0:M3)
      REAL*8 PSI_YN(0:M1,0:M2,0:M3)
      REAL*8 PSI_ZN(0:M1,0:M2,0:M3)
      REAL*8 PSI_CN(0:M1,0:M2,0:M3)
      REAL*8 VOL_TOT_ORI(MLVS)
      REAL*8 QVOL_ORI
      REAL*8 VOL1 !TOTAL FLUID VOLUME 
      
      PI=ACOS(-1.)

       CALL TIMESTAMP
       CALL REAL_TIME(TOTAL_TIME_B)
       CALL INPUT_SIMULATION(gridfile,fileprevel)
!------------------MEMORY ALLOCATAION & INITIALIZATION-----------------C      
!C----MAIN VARIABLES----------------------------------------------------
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
        U(I,J,K)=0d0
        V(I,J,K)=0d0
        W(I,J,K)=0d0
        PSI_XN(I,J,K)=0d0
        PSI_YN(I,J,K)=0d0
        PSI_ZN(I,J,K)=0d0
        PSI_CN(I,J,K)=0d0
       ENDDO
       ENDDO
       ENDDO
!$OMP PARALLEL DO
       DO LLVS=1,MLVS
       VOL_TOT_ORI(LLVS)=0.
       ENDDO      
!$OMP PARALLEL DO
       DO L=1,10000
        ITR(L)=0 ;JTR(L)=0 ;KTR(L)=0 ;ITR2(L)=0 ;JTR2(L)=0 ;KTR2(L)=0
       ENDDO
!C----MODULES-----------------------------------------------------------       
       CALL ALLOCINIT
!----------------------------------------------------------------------C
       CALL GEOM(gridfile)

      vol1=0.
!$OMP PARALLEL DO private(I,J)
!$OMP&reduction(+:VOL1)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       vol1=vol1+SDX(I)*SDY(J)*VDZ(K)
      ENDDO
      ENDDO
      ENDDO


!=====read or make initial field
      IF (IREAD.NE.0) THEN
       OPEN(12,FILE=fileprevel)
       CALL PREFLD(IHIST,fileprevel,U,V,W,P,VOL_TOT_ORI,QVOL_ORI)
       CLOSE(12)
      ELSE
       CALL MAKEFLD(IHIST,NTII,U,V,W,P)
      ENDIF

! LEVEL-SET INITIALISATION PROCESS      
      IF (ITRACKING .EQ. 0) THEN
      WRITE(*,*) 'LVS IS NOT SOLVED, ALL LVS FUNCTION = ZERO'
      ALLOCATE(SUR_X(M1M,M2M,M3M),SUR_Y(M1M,M2M,M3M),SUR_Z(M1M,M2M,M3M))
!$OMP PARALLEL DO private(I,J)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       SUR_X(I,J,K)=0.; SUR_Y(I,J,K)=0.;  SUR_Z(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO 
!-----IF THERE IS A SECOND PHASE_ITRACKING=/=0'     

      ELSE
       !NLVS=15 !MLVS
        N_MAX=M_MAX          
        CALL LVSINIT(U,V,W,PSI_CN,VOL_TOT_ORI)
        CALL GRID_COUPLING(1,PSI_XN,PSI_YN,PSI_ZN,PSI_CN)
      ENDIF ! IF (ITRACKING .EQ. 0) THEN
      
      
       IF (IPHS.EQ.1) THEN
        ALLOCATE(VOL_TOT_IPHS(MLVS))
        DO I=1,NLVS
        VOL_TOT_IPHS(I)=VOL_TOT_ORI(I)
        ENDDO
       ENDIF
      
       IF ( ICH .EQ. 0 ) THEN
         PMI=0.
       ELSE
         IF (IRESET .EQ. 1 ) THEN
         OPEN(96,FILE='0MEANP.dat')
         ELSE
         OPEN(96,FILE='0MEANP.dat',POSITION='APPEND')
         ENDIF
       ENDIF

       ALLOCATE (DENF_Z(0:M1,0:M2,0:M3))
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        DENF_Z(I,J,K)=DENM+DEN_DIFF*PSI_ZN(I,J,K)
       ENDDO
       ENDDO
       ENDDO

       DEALLOCATE (DENF_Z)

!=====initialize others
      IF (IRESET .EQ. 1) THEN
         IHIST=0
         TIME=0d0!0.01  !0.04d0 !0.02!69 !0.00861 ! 1DPART
      END IF
      NTIME=0
      
      CALL RK3COEF
      CALL LHSINIT


      NV=101                  ! INITIAL FILE NUMBER(BINARY)
      NAV=2001
      L1=3001
      L2=5001
           

!=======================================================================
!=======================================================================
!=====start time dependent calculation==================================
      CFLM=CFLMAX
      DO 3000 M=1,NTST              ! TOTAL ITERATION
      CALL real_time(TIME_BEGIN)

      NTIME=NTIME+1
!-----determine time step (DT)
      CALL DETERMINE_DT(CFLM,u,v,w,PSI_CN)           
!-----MAIN SOLVER-------------------------------------------------------------
! the IMPL input condition has been taken out to avoid confusion (2017.06.10)
        ITER_NS=3
        
      ALLOCATE(RK3XO(M1M,M2M,M3M),RK3YO(M1M,M2M,M3M),RK3ZO(M1M,M2M,M3M))
      ALLOCATE(ANPSO(M1M,M2M,M3M))
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        RK3XO(I,J,K)=0.
        RK3YO(I,J,K)=0.
        RK3ZO(I,J,K)=0.
        ANPSO(I,J,K)=0.
       ENDDO
       ENDDO
       ENDDO

      DO 2000 MSUB=1,3                    ! SUB-ITERATION(K=1,2,3)
      ALPHA_RK3=0.5*(GAMMA(MSUB)+RO(MSUB))
      ACOEF=ALPHA_RK3*DT
      ACOEFI=1./ACOEF
      DTCONST=2.*ALPHA_RK3*DT
      DTCONSTI=1./DTCONST
  
      CALL RHSNLHS(ITRACKING,U,V,W,P,PSI_XN,PSI_YN,PSI_ZN,PSI_CN
     &,RK3XO,RK3YO,RK3ZO,ANPSO,QVOL_ORI,VOL1)

      CALL real_time(ALVS_B(MSUB))
      
      ! iphs 할 때는 psi_field 가 바뀌면, sfgas 필드 구하는게 달라져서, 
      ! 이게 transport 뒤에 있으면 당연히 수렴 안한걸로 나옴...
      CALL CONVRGE1(DVMAX,U,V,W)          ! CONVERGENCY CHECK for FLOW
      IF (DVMAX.GT.RESID_POI)  WRITE(*,300) NTIME,DVMAX,TIME
 300  FORMAT(I15,'  LOCAL RESIDUE FOR MASS FLOW CONSV.:',ES18.5,
     &             ' AT ',ES11.4)      
      ! 원래는 lvs 뒤에 있었음!-----------------------------------
      
      IF (ITRACKING.NE.0) THEN
       WRITE(*,*)'>>>>>>>>TRACKING CACULATION START>>>>>>>>'  
          CALL TRANSPORT(U,V,W,PSI_CN,VOL_TOT_ORI)
          IF (IPHS .EQ. 1) DEALLOCATE(MFLVS,SGLVS)         
          CALL GRID_COUPLING(0,PSI_XN,PSI_YN,PSI_ZN,PSI_CN)
       WRITE(*,*)'<<<<<<<<<TRACKING CACULATION END<<<<<<<<<'    
       ENDIF    
      CALL real_time(ALVS_E(MSUB))
      
      
 2000  CONTINUE
      DEALLOCATE(RK3XO,RK3YO,RK3ZO)
      DEALLOCATE(ANPSO)

      IF (DENR .EQ. 1.) THEN
      DEN_DIFFI=0.
      ELSE
      DEN_DIFFI=1./DEN_DIFF
      ENDIF      

      TIME=TIME+DT

      IF (MOD(NTIME,NPIN).EQ.0) THEN
         IHIST=IHIST+1
         CALL WRITEHISTORY(CFLM,DVMAX,QMMAX)
      ENDIF
       
      IF (MOD(NTIME,NPRINT).EQ.0) THEN
         CALL WRITEFIELD(NV,IHIST,U,V,W,P,VOL_TOT_ORI,QVOL_ORI)  ! binary file
      ENDIF
      
      CALL real_time(TIME_END)

      FTRTIME1=(CONTINUITY_E(1)-CONTINUITY_B(1))
     &        +(CONTINUITY_E(2)-CONTINUITY_B(2))
     &        +(CONTINUITY_E(3)-CONTINUITY_B(3))
      FTRTIME2=0.
      FTRTIME3=0.
      FTRTIMPS=0.
      DO I=1,ITER_NS
      FTRTIME2=FTRTIME2+(ANSETIME_E(I)-ANSETIME_B(I))
      FTRTIME3=FTRTIME3+(POISSTIME_E(I)-POISSTIME_B(I))
      FTRTIMPS=FTRTIMPS+(TEMP_E(I)-TEMP_B(I))
      ENDDO
      FTRTIME4=(ALVS_E(1)-ALVS_B(1))
     &        +(ALVS_E(2)-ALVS_B(2))
     &        +(ALVS_E(3)-ALVS_B(3))

      FTRTIME5=TIME_END-TIME_BEGIN

      WRITE(*,206) FTRTIME1,FTRTIME2,FTRTIME3,FTRTIME4,FTRTIMPS,FTRTIME5
 206   FORMAT('CONTI',F6.2,'  NSE',F6.2,'  POI',F6.2,'  LVS',F6.2,
     & '  TMP',F6.2,'  TOT',F6.2)
      WRITE(87,211)TIME,FTRTIME1,FTRTIME2,FTRTIME3,FTRTIME4,FTRTIME5

 3000 CONTINUE

      CLOSE(100)
      CLOSE(13)
      IF (ICH   .NE. 0) CLOSE(96)
      IF (IPS   .EQ. 1) CLOSE(97)
      CLOSE(87)

      IF (MOD(NTIME,NPRINT).NE.0) THEN
         CALL WRITEFIELD(NV,IHIST,U,V,W,P,VOL_TOT_ORI,QVOL_ORI)
      ENDIF

      CALL real_time(TOTAL_TIME_E)
      WRITE(*,*) ' '
      WRITE(*,210) TOTAL_TIME_E-TOTAL_TIME_B
 210  FORMAT('TOTAL TIME      : ',F12.2,' SECONDS')
 211  FORMAT(F13.5,5F12.4)

      call timestamp
      write(*,*) 'End of simulation'
      STOP
      END

!===============================================================
!     SUBROUTINS FOR INPUT AND OUTPUT
!===============================================================

!**********************************************************************
      SUBROUTINE INPUT_SIMULATION(gridfile,fileprevel)
!**********************************************************************      
      
      USE FLOW_VAR
      USE FLOW_GEOM_VAR
      USE FLD_AVG
      USE MG_OPTION
      USE HEAT_VAR
      USE LES_VAR
      USE IBM_VAR
      USE LVS_VAR
      USE TWO_PHASE_PROPERTY 
      
      USE PHASE_BOILING
      IMPLICIT NONE
      
      INTEGER*8   IERR
      
      CHARACTER*10 DUMMY
      CHARACTER*30 gridfile,fileprevel
      REAL*8       PI
      PI=ACOS(-1.)
      
!C=====read input file
      OPEN(10,FILE='lica.in')
      READ(10,*) DUMMY
      READ(10,*) IRESET,IREAD,IAVG,IPZERO,IDTOLD
      READ(10,*) DUMMY
      READ(10,*) NTST,NPRINT,NPRIAVG,NPIN
      READ(10,*) DUMMY
      READ(10,*) IMPL,IDTOPT,DT,CFLMAX
      READ(10,*) DUMMY
      READ(10,*) RESID_POI,NLEV,IWC,NBLI,MGITR,IOLDV,IMGSOR,WWSOR,IPCG
      READ(10,*) DUMMY
      READ(10,*) IPS,IPHS
      ! READ(10,*) DUMMY
      ! READ(10,*) ILES,ISGSMOD,DYNMON,CSGS,CSGSPS,IFILTER
      READ(10,*) DUMMY
      READ(10,*) IBMON,MASSON
      READ(10,*) DUMMY
      READ(10,*) IPOISS
      READ(10,*) DUMMY
      READ(10,*) ICH,IPX,IPY,IPZ
      READ(10,*) DUMMY
      READ(10,*) DUMMY
      READ(10,*) ITRACKING,MCLS,MGLOBAL
      READ(10,*) DUMMY
      READ(10,*) RE_AIR,VISR,DENR
      READ(10,*) DUMMY
      READ(10,*) SURF_J,FR,FK
      READ(10,*) DUMMY
      READ(10,*) PRM,SCR,TCR
      READ(10,*) DUMMY
      READ(10,*) T_SAT,JA_AIR
      READ(10,*) DUMMY
      READ(10,*) DUMMY
      READ(10,'(a)') gridfile
      READ(10,'(a)') fileprevel
      READ(10,*) DUMMY
      READ(10,*) NTEST,NTPRINT
            READ(10,*) DUMMY
      READ(10,*) NLVS
            READ(10,*) DUMMY
      READ(10,*) ilag

      WRITE(*,101) IRESET,IREAD,IAVG,IDTOLD
      WRITE(*,103) NTST,NPRINT,NPRIAVG,NPIN
      WRITE(*,105) IMPL,IDTOPT,DT,CFLMAX
      WRITE(*,109) RESID_POI,NLEV,IWC,NBLI,MGITR,IOLDV,IMGSOR,WWSOR
      WRITE(*,110) IPS
      WRITE(*,111) ILES,ISGSMOD,DYNMON,CSGS,CSGSPS,IFILTER
      WRITE(*,113) IBMON,MASSON
      WRITE(*,115) IPOISS
      WRITE(*,117) ICH,IPX,IPY,IPZ
      WRITE(*,119) ITRACKING,MCLS,MGLOBAL
      WRITE(*,121) RE_AIR,VISR,DENR
      WRITE(*,123) SURF_J,FR,FK
      WRITE(*,125) gridfile
      IF (IREAD .EQ. 1) WRITE(*,127) fileprevel
     
 101  FORMAT('IRESET=',I8,'  IREAD=',I9,'  IAVG=',I10,
     &       '  IPZERO=',I8,'  IDTOLD=',F7.2)
 103  FORMAT('NTST=',I10,'  NPRINT=',I8,'  NPRIAVG=',I7,'  NPIN=',I10)
 105  FORMAT('IMPL=',I10,'  IDTOPT=',I8,'  DT=',ES12.4,
     &       '  CFLMAX=',F8.2)
!107  FORMAT('IUD=',I11,'  IEND=',I10,'  IENDM=',I9,'  ALPZ=',F10.4)
 109  FORMAT('RESID=',ES9.2,'  NLEV=',I10,'  IWC=',I11,'  NBLI=',I10
     &      ,'  MGITR=',I9,'  IOLDV=',I9,'  IMGSOR=',I8,'  WWSOR=',F9.3)
 110  FORMAT('IPS=',I11)
 111  FORMAT('ILES=',I10,'  ISGSMOD=',I7,'  DYNMON=',L8,
     &       '  CSGS=',F10.4,'  CSGSPS=',F8.4,'  IFILTER=',I7)
 113  FORMAT('IBMON=',I9,'  MASSON=',I8)
 115  FORMAT('IPOISS=',I9)
 117  FORMAT('ICH=',I9,'  IPX=',I8,'  IPY=',I8,'  IPZ=',I8)
 119  FORMAT('ITRACKING=',I9,'  MCLS=',I8,'  MGLOBAL=',I8)
 121  FORMAT('RE_AIR=',F12.5,'  VISR=',F12.5,'  DENR=',F12.5)
 123  FORMAT('SUR_J=',F12.5,'  FR=',F12.5,'  FK=',F12.5)
 125  FORMAT('GRID=   ',A30)
 127  FORMAT('FIELD_READ=   ',A30) 

!100  FORMAT('INITIAL FIELD TREATMNT:                READING DONE')

       IF (ICH .EQ. 0 ) THEN
       WRITE(*,*) 'NO-MEAN_PRESSURE GRADIENT'
       ELSE IF (ICH .EQ. 1 ) THEN
       WRITE(*,*) 'MASS-FLOW RATE CONST'
       ELSE IF (ICH .EQ. 2 ) THEN
       WRITE(*,*) 'PRESSURE GRADIENT CONST'
       ENDIF

!============================FOR TWO-PHASE=============================C
      RE=RE_AIR
      REI=1./RE

      DENM=1.        !NON-DIMENSIONALIZED VIS
      DENP=DENR
      DEN_DIFF=DENP-DENM

      VISM=1.*REI    !NON-DIMENSIONALIZED VIS
      VISP=VISR*REI
      VIS_DIFF=VISP-VISM

      REM=RE_AIR
      REP=REM*DENR/VISR
      
      DENSCM=DENM*1.   !SPECIFIC HEAT CAPACITY
      DENSCP=DENSCM*(DENR*SCR)
      DENSC_DIFF=DENSCP-DENSCM
      
      TCM=1.   !THERMAL CONDUCTIVITY
      TCP=TCR
      TC_DIFF=TCP-TCM
      
      PRM=PRM !PRANDTL NUMBER
      PRP=PRM*VISR*SCR/TCR
      
      PRA=PRM
      PRAI=1./PRM
      IF ( SURF_J .EQ. 0. ) THEN
       SURF=0.
      ELSE
       SURF=1./SURF_J  !SURF_J IS WEBER NUMBER
      ENDIF

      IF (ICH .EQ. 2) THEN
       PMI_CONST=-FK
      ELSE
        PMI_CONST=0.
      ENDIF

      IF(FR .EQ. 0.) THEN
        GRAVITY=0.
      ELSE
        GRAVITY=1./FR**2
      ENDIF

      TIME_ST=SURF*PI/(DENM+DENP)
      TIME_GV=DEN_DIFF/(DENM+DENP)*GRAVITY/PI

      PRINT*,'MEAN PRESSURE GRADIENT :',FK
      PRINT*,'SURF_J :',SURF_J
      PRINT*,'FR :',FR
      PRINT*,'DENSITY_RATIO :',DENR
      PRINT*,'VISCOSITY_RATIO :',VISR
      PRINT*,'RE_AIR :',REM,'RE_WATER :',REP
      PRINT*,'PR_AIR :',PRM,'PR_WATER :',PRP
!============================FOR TWO-PHASE=============================C     
     
      RETURN      
      END
      
      
!*******************************************************************
      SUBROUTINE ALLOCINIT
!*******************************************************************

      USE PARAM_VAR     
      USE FLOW_VAR
      USE FLOW_GEOM_VAR
      USE FLD_AVG      
      USE HEAT_VAR
      USE IBM_VAR
      !USE LES_VAR            
      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING      
      USE TIME_VAR
      USE PHASE_BOILING
      
      IMPLICIT NONE  
      INTEGER*8   I,J,K,L,N,NN,NNVS,LLVS,IP1,JP1          
      
!------------------MEMORY ALLOCATAION & INITIALIZATION-----------------C
!-----------------------------------------------------------------------      
!      MODULE FLOW_VAR
!-----------------------------------------------------------------------   
      ALLOCATE(UO(0:M1,0:M2,0:M3),VO(0:M1,0:M2,0:M3),WO(0:M1,0:M2,0:M3))
      ALLOCATE(AIU(M1),CIU(M1),AIVW(M1),CIVW(M1)
     &        ,AJV(M2),CJV(M2),AJUW(M2),CJUW(M2)
     &        ,AKW(M3),CKW(M3),AKUV(M3),CKUV(M3))

!$OMP PARALLEL DO private(I,J)
       DO K=0,N3 ;
       DO J=0,N2
       DO I=0,N1
        UO(I,J,K)=0. ;VO(I,J,K)=0.; WO(I,J,K)=0.
       ENDDO ;ENDDO ;ENDDO
!$OMP PARALLEL DO
       DO I=1,N1
        AIU(I)=0. ;CIU(I)=0. ;AIVW(I)=0. ;CIVW(I)=0.
       ENDDO
!$OMP PARALLEL DO
       DO J=1,N2
        AJV(J)=0. ;CJV(J)=0. ;AJUW(J)=0. ;CJUW(J)=0.
       ENDDO
!$OMP PARALLEL DO
       DO K=1,N3
        AKW(K)=0. ;CKW(K)=0. ;AKUV(K)=0. ;CKUV(K)=0.
       ENDDO
       
!-----------------------------------------------------------------------      
!      MODULE LVS_VAR
!----------------------------------------------------------------------- 
      ALLOCATE( ALPHI(MF_BAND,M_MAX,MLVS))
      ALLOCATE( NUMA(M_MAX,MLVS))
      ALLOCATE( I_B(MF_BAND,M_MAX,MLVS),J_B(MF_BAND,M_MAX,MLVS)
     &                           ,K_B(MF_BAND,M_MAX,MLVS))
      ALLOCATE( PSIF(M1L:M1U,M2L:M2U,M3L:M3U))

      ALLOCATE(MASK_BUB(M1M,M2M,M3M))
      ALLOCATE(MASK_GLOBAL(0:M1F,0:M2F,0:M3F))     
      ALLOCATE( ALPHI_COL(M1L:M1U,M2L:M2U,M3L:M3U))                    ! 2017-08-20
      ALLOCATE( MASK_COL(0:M1F,0:M2F,0:M3F))                           ! 2017-08-20
      ALLOCATE( CONTACT(MLVS,MLVS),CONTACT_OLD(MLVS,MLVS))             ! 2017-08-16 
      ALLOCATE( CONTACT_START(MLVS,MLVS),CONTACT_TIME(MLVS,MLVS))      ! 2017-08-16 
      ALLOCATE( VELVS(MLVS,3))                                         ! 2017-09-24
      ALLOCATE( T_CON(MLVS,MLVS))                                      ! 2017-09-24
      ALLOCATE( LVSON(MLVS))
      ALLOCATE( XLAGPRT(10000,4), uLAGPRT(10000,6), vlagprt(10000,3)) !20171222, 2018-01-07
      ALLOCATE( DVDT(10000,6),DUDX(6,10000,3))
            
       XLAGPRT=0.                                   !20171222
       uLAGPRT=0.                                   !20171222
       VLAGPRT=0.                                   !20171222
       NPRT=0                                       !20171222
       DVDT=0.                                      !20180107
      
!$OMP PARALLEL DO private(NNVS)
       DO LLVS=1,MLVS
       DO NNVS=1,MLVS
       CONTACT(NNVS,LLVS)=0                          ! 2017-08-16   INT
       CONTACT_OLD(NNVS,LLVS)=0                      ! 2017-08-16   INT
       CONTACT_START(NNVS,LLVS)=0D0                  ! 2017-08-16   DBLE
       CONTACT_TIME(NNVS,LLVS)=0D0                   ! 2017-08-16   DBLE
       t_con(nnvs,llvs)=100.                         ! 2017-09-24
       ENDDO ;       ENDDO 
       
      DO LLVS=1,NLVS!  MLVS
      LVSON(LLVS)=1    ! this is set to MLVS 2017-10-08
      ENDDO

       DO LLVS=1,MLVS ! this is set to MLVS 2017-10-08
       DO N=1,M_MAX
        NUMA(N,LLVS)=0
!$OMP PARALLEL DO
       DO NN=1,MF_BAND
        ALPHI(NN,N,LLVS)=0.
        I_B(NN,N,LLVS)=0 ;J_B(NN,N,LLVS)=0 ;K_B(NN,N,LLVS)=0
       ENDDO;       ENDDO;       ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=M3L,M3U
       DO J=M2L,M2U
       DO I=M1L,M1U
        PSIF(I,J,K)=0.
        ALPHI_COL(I,J,K)=0.                       !2017-08-20
       ENDDO;       ENDDO;       ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        MASK_BUB(I,J,K)=0
       ENDDO;       ENDDO;       ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=0,N3F
       DO J=0,N2F
       DO I=0,N1F
        MASK_GLOBAL(I,J,K)=0
        MASK_COL(I,J,K)=0                         !2017-08-20
       ENDDO;       ENDDO;       ENDDO
       
!-----------------------------------------------------------------------  
!      MODULE IBM_VAR
!-----------------------------------------------------------------------
      ALLOCATE(QMASS(0:M1,0:M2,0:M3))
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
        QMASS(I,J,K)=0.
       ENDDO
       ENDDO
       ENDDO

!-----------------------------------------------------------------------        
!      MODULE HEAT_VAR
!----------------------------------------------------------------------- 
      IF (IPS .EQ. 1) THEN
      ALLOCATE(T(0:M1,0:M2,0:M3))
      ALLOCATE(TALPH(0:M1,0:M2,0:M3))
      ALLOCATE(RHSPS(M1M,M2M,M3M))
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3 ;DO J=0,N2 ;DO I=0,N1
        T(I,J,K)=0D0
        TALPH(I,J,K)=0.
       ENDDO ;ENDDO ;ENDDO

!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M ;DO J=1,N2M ;DO I=1,N1M
        RHSPS(I,J,K)=0.
       ENDDO ;ENDDO ;ENDDO
      ENDIF
      
!-----------------------------------------------------------------------
!      MODULE FLOW_GEOM_VAR
!-----------------------------------------------------------------------  
      ALLOCATE(!SGSTIME_B(3),SGSTIME_E(3),
     &    ANSETIME_B(100),ANSETIME_E(100),POISSTIME_B(3),POISSTIME_E(3)
     &   ,CONTINUITY_B(3),CONTINUITY_E(3),ALVS_B(3),ALVS_E(3))
      allocate(temp_b(3),temp_E(3))
      
       DO L=1,3
        !SGSTIME_B(L)=0. ;SGSTIME_E(L)=0.
        POISSTIME_B(L)=0d0    ;POISSTIME_E(L)=0d0
        CONTINUITY_B(L)=0d0   ;CONTINUITY_E(L)=0d0
        ALVS_B(L)=0d0         ;ALVS_E(L)=0d0
        temp_B(L)=0d0         ; temp_E(L)=0d0
       ENDDO
      
       DO L=1,100
        ANSETIME_B(L)=0. ;ANSETIME_B(L)=0
       ENDDO
!-----------------------------------------------------------------------        
!      MODULE PHASE_BOILING
!-----------------------------------------------------------------------
      IF (IPHS.EQ. 1) THEN 
      ALLOCATE(TEMP(0:M1,0:M2,0:M3,0:1))
      ALLOCATE(SFGAS(M1L:M1U,M2L:M2U,M3L:M3U))
      ALLOCATE(SFGAS2(M1L:M1U,M2L:M2U,M3L:M3U))
      ALLOCATE(IMASK_PH(0:M1,0:M2,0:M3))
      
      ENDIF       
!C------------------MEMORY ALLOCATAION & INITIALIZATION-----------------C

      RETURN
      END
      
!*******************************************************************
      SUBROUTINE DETERMINE_DT(CFLM,u,v,w,PSI_CN)
!*******************************************************************
      USE PARAM_VAR
      USE FLOW_VAR
      USE FLOW_GEOM_VAR
      USE TWO_PHASE_PROPERTY 
      

      IMPLICIT NONE
      INTEGER*8   I,J,K
      REAL*8      U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)    
      REAL*8      PSI_CN(0:M1,0:M2,0:M3)      
      REAL*8      DT_OLD,DT_SOURCE,DT_SOURCE_TMP,DX1,DX2,DT_TMP1,DT_TMP2
      REAL*8      CFL_SURF,CFL_SUM,CFL_MAX_VAR,
     &            CFLM,CFL_MAX,CFLCONVMAX,CFLSURFMAX


      DT_OLD=DT
      !this is coded refering 'capillary wave' at wikipedia
      IF( TIME_ST .EQ. 0. ) THEN
      DT_SOURCE=1000000.  !VERY SMALL NUMBER FOR AVODING DEVIDED BY ZERO
      ELSE
      DT_SOURCE=1000.
      
      
       IF (N3M .EQ. 1 ) THEN
      K=1
!$OMP PARALLEL DO private(I,DX1,DX2,DT_TMP1,DT_TMP2)
!$OMP&firstprivate(K)
!$OMP&reduction(MIN:DT_SOURCE)
      DO J=1,N2M
      DO I=1,N1M
        IF ( DABS(PSI_CN(I,J,K)-0.5) .LT. 0.5 ) THEN
           DX1=DMIN1(SDX(I),SDY(J))
           DT_TMP1=0.5*DX1/DSQRT(TIME_ST/DX1+TIME_GV*DX1)
           DX2=DMAX1(SDX(I),SDY(J))
           DT_TMP2=0.5*DX2/DSQRT(TIME_ST/DX2+TIME_GV*DX2)
           DT_SOURCE=DMIN1(DT_SOURCE,DT_TMP1,DT_TMP2)
        ENDIF
      ENDDO
      ENDDO
      
       ELSE IF (N1M .EQ. 1 ) THEN
      I=1
!$OMP PARALLEL DO private(J,DX1,DX2,DT_TMP1,DT_TMP2)
!$OMP&firstprivate(I)
!$OMP&reduction(MIN:DT_SOURCE)
      DO K=1,N3M
      DO J=1,N2M
        IF ( DABS(PSI_CN(I,J,K)-0.5) .LT. 0.5 ) THEN
           DX1=DMIN1(SDY(J),SDZ(K))
           DT_TMP1=0.5*DX1/DSQRT(TIME_ST/DX1+TIME_GV*DX1)
           DX2=DMAX1(SDY(J),SDZ(K))
           DT_TMP2=0.5*DX2/DSQRT(TIME_ST/DX2+TIME_GV*DX2)
           DT_SOURCE=DMIN1(DT_SOURCE,DT_TMP1,DT_TMP2)
        ENDIF
      ENDDO
      ENDDO
      ELSE IF ((N2M .EQ. 1 ) .and. (n3m .eq. 1)) THEN  ! 1d equation
      k=1
      j=1
      DO i=1,N1M
        IF ( DABS(PSI_CN(I,J,K)-0.5) .LT. 0.5 ) THEN
           DX1=DMIN1(SDY(J),SDZ(K))
           DT_TMP1=0.5*DX1/DSQRT(TIME_ST/DX1+TIME_GV*DX1)
           DX2=DMAX1(SDY(J),SDZ(K))
           DT_TMP2=0.5*DX2/DSQRT(TIME_ST/DX2+TIME_GV*DX2)
           DT_SOURCE=DMIN1(DT_SOURCE,DT_TMP1,DT_TMP2)
        ENDIF
      ENDDO
      eLSE
!$OMP PARALLEL DO private(I,J,DX1,DX2,DT_TMP1,DT_TMP2)
!$OMP&reduction(MIN:DT_SOURCE)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
        IF ( DABS(PSI_CN(I,J,K)-0.5) .LT. 0.5 ) THEN
             DX1=DMIN1(SDX(I),SDY(J),SDZ(K))
           DT_TMP1=0.5*DX1/DSQRT(TIME_ST/DX1+TIME_GV*DX1)
             DX2=DMAX1(SDX(I),SDY(J),SDZ(K))
           DT_TMP2=0.5*DX2/DSQRT(TIME_ST/DX2+TIME_GV*DX2)
           DT_SOURCE=DMIN1(DT_SOURCE,DT_TMP1,DT_TMP2)
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDIF

      ENDIF

      CALL CFL(CFLM,U,V,W)
      
      CFLCONVMAX=CFLMAX
      CFLSURFMAX=CFLMAX
      CFL_SURF=DT/DT_SOURCE
      
      
      IF (IDTOPT .EQ. 1) THEN
!-----DT_DETERMINE 1ND TYPE
        CFL_SUM=0.5*(CFLM+DSQRT(CFLM**2+4*CFL_SURF))
        DT=DT*CFLCONVMAX/CFL_SUM
!-----DT_DETERMINE 2ND TYPE
!        CFL_SUM=CFLM+CFL_SURF
!        DT=DT*AMIN1(CFLCONVMAX/CFLM,CFLSURFMAX/CFL_SURF)
      else 
      cfl_sum = 0.5d0*(cflm+cfl_surf) ! temporary: no background
      ENDIF

      IF ( DT/DT_OLD .GT. 2. ) THEN
        DT=DT_OLD*1.2    !CHECK LATTER, SMOOTH STARTING WHEN INITIAL VEL. IS ZERO.
        WRITE(*,*) 'DT CHANGE IS LARGE ALGORITHM IS ACTIVATED'
      ENDIF

      CFL_MAX=DMAX1(CFLM,CFL_SURF)
      IF (CFL_MAX .GT.(CFLMAX*1.1)) THEN
      PRINT*,' '
!      WRITE(*,310) NTIME,CFLM,TIME
      WRITE(*,320) NTIME,TIME,DT
      WRITE(*,*) 'CFL NUMBER IS EXCEEDED!!'
      ELSE
      PRINT*,' '
      WRITE(*,320) NTIME,TIME,DT
      ENDIF
      !WRITE(*,153) CFLM,CFL_SURF,CFL_SUM
      WRITE(*,*) CFLM,CFL_SURF,CFL_SUM
      
      
 310  FORMAT(I15,'  CFL# EXCEED GIVEN CFL LIMIT :',ES18.5,
     &             ' AT TIME',F12.5)
 320  FORMAT('--------------------------',I6,'  TIME=',F10.5,
     &             '  DT=',F12.8)
 153  FORMAT('CFL_CON= ',F7.4,'  CFL_SURF= ',F7.4,'  CFLSUM= ',F7.4)
      
      RETURN
      END


!*******************************************************************
      SUBROUTINE MAKEFLD(IHIST,NTII,U,V,W,P)
!*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      USE HEAT_VAR
      use phase_boiling
      
      IMPLICIT NONE

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 P(0:M1,0:M2,0:M3)
      
      INTEGER*8 IHIST,NTII      
      
      INTEGER*8 I,J,K,ICHANNEL,IDUM
      REAL*8 PI,RCOEF,RE_WATER,EPS_PTR,MFD_TURB,UMAX,RR,Y_PLUS,W1,W2
      REAL*8 u2,v2
      REAL*8 RAN1,RAN_NUM1,EVMM,EVM,EVM_DIVIDE
      REAL*8 FUNCBODY,PBODY,ASSIGN_LVS
      REAL*8 UU,VV,WW,U_BULK
      REAL*8 VOL,DVOL
      REAL*8 DTDX,COEF,AA,BB,CC
      
      real*8 xc,sigma1
      integer*8 N
      
      real*8 xzero

       INTEGER*8 LLVS
      PI=ACOS(-1.)

        U=0d0
        V=0d0
        W=0d0
        
        W2=0.
        u2=0.
        v2=0.
        IHIST=0
        NTII=0
        TIME=0.
        
!$OMP PARALLEL DO private(I,J)
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         P(I,J,K)=0.  
        ENDDO
        ENDDO
        ENDDO
        
!      
!      xzero=0.4998d0
!       !xzero=0.5002d0
!         !IF (IPS.EQ.1) CALL MAKEFLD_TEMP
!         if (IPHS .EQ. 1) THEN 
!!$OMP PARALLEL DO private(PBODY,J,I)         
!          DO K=1,N3M
!          DO J=1,N2M
!          DO I=1,N1M
!           PBODY=ASSIGN_LVS(XP(I),YP(J),ZP(K),1,PI)
!            if (pbody .lt. 0d0) then 
!            T(i,j,k) =1d0 ! 1d0/xzero*zp(k) ! 
!            ELSEIF (PBODY .gt. 0D0) THEN 
!            T(I,J,K) =-1d0/(1d0-xzero)*(zp(k)-xzero)+1d0  !1d0             ! 
!                                             
!                       
!            ENDIF
!          ENDDO
!          ENDDO
!          ENDDO
!         ENDIF
       if (IPHS .EQ. 1) THEN 
!$OMP PARALLEL DO private(PBODY,J,I)         
          DO K=1,N3M
          DO J=1,N2M
          DO I=1,N1M
           PBODY=ASSIGN_LVS(XP(I),YP(J),ZP(K),1,PI)
           if (pbody .gt. 0d0) then 
         T(i,j,k) =1d0 !-zp(k)!
         w(i,j,k) =1d0 !76175d0
         w(i,j,n3)=1d0 ! 0.76175d0
         ELSEIF (PBODY .Le. 0D0) THEN 
         ! T(I,J,K)=1d0-18.39588d0*zp(k)
           T(I,J,K)=16.47d0*zp(k) !23.2d0*zp(k)!1d0-  23.2d0*zp(k)!
           T(I,J,0)=0D0
           w(i,j,k)=0d0
           w(i,j,1)=0d0
           ENDIF
          ENDDO
          ENDDO
          ENDDO
         ENDIF
 
 
 
        OPEN(146,FILE='0INITIAL_VEL.DAT')
        IF (IPS.EQ. 1) THEN 
        WRITE(146,*) 'VARIABLES="X","Y","Z","U","V","W","T"'
        ELSEIF (IPS .EQ. 0) THEN 
        WRITE(146,*) 'VARIABLES="X","Y","Z","U","V","W"'
        ENDIF
       WRITE(146,*) 'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
          UU=0.5*(U(I,J,K)+U(IPV(I),J,K))
          VV=0.5*(V(I,J,K)+V(I,JPV(J),K))
          WW=0.5*(W(I,J,K)+W(I,J,KPV(K)))
      IF (IPS.EQ. 1) THEN 
      WRITE(146,147) XP(I),YP(J),ZP(K),UU,VV,WW,T(I,J,K)
      ELSE
      WRITE(146,147) XP(I),YP(J),ZP(K),UU,VV,WW
      ENDIF
      ENDDO
       ENDDO
       ENDDO
        CLOSE(146)
  147  FORMAT(6F15.8)
      !STOP
        RETURN
        END    
!C*******************************************************************
      SUBROUTINE MAKEFLD_TEMP
!C*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE HEAT_VAR

      IMPLICIT NONE

      INTEGER*8 I,J,K
      REAL*8 FUNCBODY,RR
      REAL*8 DTDX,U_BULK,COEF,AA,BB,CC

       T=0.

       
        RETURN
        END
        
 
!*******************************************************************
      SUBROUTINE PREFLD(IHIST,fileprevel,U,V,W,P,
     & VOL_TOT_ORI,QVOL_ORI)
!*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR

      USE HEAT_VAR

      CHARACTER*30 fileprevel

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 P(0:M1,0:M2,0:M3)

      REAL*8 VOL_TOT_ORI(MLVS)
      
      REAL*8 QVOL_ORI

      INTEGER*8 LVS,LVS2
      INTEGER*8 ihist
      
!     dum for future use
      READ(12,1001) N1,N2,N3
      READ(12,1002) IHIST,M,TIME,DT_O
      READ(12,1003) IPSS,IIPX,IIPY,IIPZ,PRAA
      READ(12,1004) QVOL_ORI
      READ(12,1005) ((( U(I,J,K) ,I=1,N1),J=0,N2),K=0,N3)
      READ(12,1005) ((( V(I,J,K) ,I=0,N1),J=1,N2),K=0,N3)
      READ(12,1005) ((( W(I,J,K) ,I=0,N1),J=0,N2),K=1,N3)
      READ(12,1005) ((( P(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)

      !LVS
      READ(12,1011) N1F,N2F,N3F
      READ(12,1012) DENR,VISR,FR,SURF_J,FK,RE_AIR

      READ(12,1013) NLVS,N_MAX
      READ(12,1014) (VOL_TOT_ORI(LLVS), LLVS=1,NLVS)
      READ(12,1015) ((NUMA(N,LLVS) ,N=1,N_MAX),LLVS=1,NLVS)
      DO LLVS=1,NLVS
       DO N=1,N_MAX
        DO NN=1,NUMA(N,LLVS)
        READ(12,1016) I,J,K,ALPHI_TMP
         I_B(NN,N,LLVS)=I
         J_B(NN,N,LLVS)=J
         K_B(NN,N,LLVS)=K
         ALPHI(NN,N,LLVS)=ALPHI_TMP
        ENDDO
       ENDDO
      ENDDO
      
      
      !2017-09-03 CONTACT-RELATED SAVE      !2017-09-03 CONTACT-RELATED SAVE
       READ(12,1007) ((CONTACT(LVS,LVS2)      ,LVS=1,NLVS),LVS2=1,NLVS)
       READ(12,1007) ((CONTACT_OLD(LVS,LVS2)  ,LVS=1,NLVS),LVS2=1,NLVS)
       READ(12,1005) ((CONTACT_START(LVS,LVS2),LVS=1,NLVS),LVS2=1,NLVS)
       READ(12,1005) ((CONTACT_TIME(LVS,LVS2) ,LVS=1,NLVS),LVS2=1,NLVS)
       READ(12,1007) ( LVSON(LVS) ,LVS=1,NLVS)
       READ(12,1005) ((T_CON(LVS,LVS2)        ,LVS=1,NLVS),LVS2=1,NLVS)
      !2017-09-03 CONTACT-RELATED SAVE      !2017-09-03 CONTACT-RELATED SAVE   

      
      IF (IPS.EQ.1) THEN
      IF (IPSS.EQ.0) THEN
       CALL MAKEFLD_TEMP
      ELSE IF (IPSS.EQ.1) THEN
       READ(12,1006) PRM,SCR,TCR
       READ(12,1005) ((( T(I,J,K) ,I=0,N1),J=0,N2),K=0,N3)
      ENDIF
      ENDIF

 1001   FORMAT(3I8)
 1002   FORMAT(2I8,2ES20.12)
 1003   FORMAT(4I8,1ES20.12)
 1004   FORMAT(1ES20.12)
 1005   FORMAT(5ES20.12)
 1006   FORMAT(3ES20.12)
 1007   FORMAT(4I8)
 
 1011   FORMAT(3I8)
 1012   FORMAT(6ES20.12)
 1013   FORMAT(2I8)
 1014   FORMAT(8ES20.12)
 1015   FORMAT(8I8)
 1016   FORMAT(3I8,1ES20.12)


      IF (IPZERO .EQ. 1) P=0.

      WRITE(*,*)' '
      WRITE(*,100)
      WRITE(*,101)
      WRITE(*,102) fileprevel
      WRITE(*,104) N1,N2,N3
      WRITE(*,105) IHIST,M,TIME,DT_O
      WRITE(*,106) IPSS,PRAA

  100 FORMAT('----------- INITIAL FIELD INFORMATION -----------')
  101 FORMAT('INITIAL FIELD      : READING DONE')
  102 FORMAT('INITIAL FIELD NAME : ',A30)
  104 FORMAT('N1=',I12,'  N2=',I12,'  N3=',I12)
  105 FORMAT('IHIST=',I9,'  M=',I13,'  TIME=',F10.5,'  DT=',F12.8)
  106 FORMAT('IPS=',I11,'  PRA=',F11.3)

       IF (IDTOLD .EQ. 1) THEN
        DT=DT_O
        WRITE(*,*) 'DT_OLD IS USED'
       ELSE
        WRITE(*,*) 'DT_NEW IS USED'
       ENDIF

      DO LLVS=1,NLVS
       NUMA_MAX=0
       DO N=1,N_MAX
         NUMA_MAX=MAX(NUMA_MAX,NUMA(N,LLVS))
       ENDDO
        WRITE(*,*) 'NUMA_MAX/MF_BAND=',FLOAT(NUMA_MAX)/FLOAT(MF_BAND)
      ENDDO

      
      RETURN
      END
!*******************************************************************
        FUNCTION FUNCBODY(X_IBM,Y_IBM,Z_IBM)
!*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      USE IBM_VAR



        IF ( IBMON .EQ. 1 ) THEN

       FUNCBODY=-(X_IBM**2.+Y_IBM**2.-(0.5*YL)**2.)


        ELSE
          FUNCBODY=100 !FOR NO IBM CASE.

        ENDIF

        RETURN
        END
!C*******************************************************************
      SUBROUTINE WRITEFIELD(NV,IHIST,U,V,W,P,VOL_TOT_ORI,QVOL_ORI)
!C*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR

      USE HEAT_VAR
      
      USE IBM_VAR
      
      IMPLICIT NONE

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 P(0:M1,0:M2,0:M3)

      REAL*8 VOL_TOT_ORI(MLVS)
      REAL*8 QVOL_ORI
      
      INTEGER*8 NV,IHIST
      
      INTEGER*8 IDUM
      REAL*8 DUM
      INTEGER*8 I,J,K,N,NN,LLVS,L,LVS,LVS2 ! LVS1,LVS2 ADDED 2017-09-03
      CHARACTER*9 tname
      CHARACTER*3 tfn1
      INTEGER*8 IDG1,IDG2,IDG3,IDG4,IDG5,IDG6

      idum=0
      dum=0.

      tfn1='fld'
      idg1=ihist/100000
      idg2=(ihist-idg1*100000)/10000
      idg3=(ihist-idg1*100000-idg2*10000)/1000
      idg4=(ihist-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6=ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
      tname=tfn1//char(idg1+48)//char(idg2+48)//
     &      char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)

      OPEN(NV,FILE=tname)
!ccc      IVER=-2
!ccc      WRITE(NV) IVER
      WRITE(NV,1001) N1,N2,N3
      WRITE(NV,1002) IHIST,M,TIME,DT
      WRITE(NV,1003) IPS,IPX,IPY,IPZ,PRA
      WRITE(NV,1004) QVOL_ORI
      WRITE(NV,1005) ((( U(I,J,K) ,I=1,N1),J=0,N2),K=0,N3)
      WRITE(NV,1005) ((( V(I,J,K) ,I=0,N1),J=1,N2),K=0,N3)
      WRITE(NV,1005) ((( W(I,J,K) ,I=0,N1),J=0,N2),K=1,N3)
      WRITE(NV,1005) ((( P(I,J,K) ,I=1,N1M),J=1,N2M),K=1,N3M)

      !LVS
      WRITE(NV,1011) N1F,N2F,N3F
      WRITE(NV,1012) DENR,VISR,FR,SURF_J,FK,RE_AIR

      WRITE(NV,1013) NLVS,N_MAX
      ! do llvs=1,nlvs
      ! WRITE(*,*)'vol_tot_ori=',VOL_TOT_ORI(LLVS)
      ! enddo
      WRITE(NV,1014) (VOL_TOT_ORI(LLVS), LLVS=1,NLVS)
      WRITE(NV,1015) ((NUMA(N,LLVS) ,N=1,N_MAX),LLVS=1,NLVS)
      DO LLVS=1,NLVS
       DO N=1,N_MAX
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
        WRITE(NV,1016) I,J,K,ALPHI(NN,N,LLVS)
        ENDDO
       ENDDO
      ENDDO
      
      !2017-09-03 CONTACT-RELATED SAVE      !2017-09-03 CONTACT-RELATED SAVE
      WRITE(NV,1007) ((CONTACT(LVS,LVS2)      ,LVS=1,NLVS),LVS2=1,NLVS)
      WRITE(NV,1007) ((CONTACT_OLD(LVS,LVS2)  ,LVS=1,NLVS),LVS2=1,NLVS)
      WRITE(NV,1005) ((CONTACT_START(LVS,LVS2),LVS=1,NLVS),LVS2=1,NLVS)
      WRITE(NV,1005) ((CONTACT_TIME(LVS,LVS2) ,LVS=1,NLVS),LVS2=1,NLVS)
      WRITE(NV,1007) ( LVSON(LVS) ,LVS=1,NLVS)
      WRITE(NV,1005) ((T_CON(LVS,LVS2)        ,LVS=1,NLVS),LVS2=1,NLVS)
      !2017-09-03 CONTACT-RELATED SAVE      !2017-09-03 CONTACT-RELATED SAVE
      IF (IPS.EQ.1) THEN
      WRITE(NV,1006) PRM,SCR,TCR
      WRITE(NV,1005) ((( T(I,J,K) ,I=0,N1),J=0,N2),K=0,N3)
      ENDIF
      WRITE(*,1007) ( LVSON(LVS) ,LVS=1,NLVS) 
      CLOSE(NV)

 1001   FORMAT(3I8)
 1002   FORMAT(2I8,2ES20.12)
 1003   FORMAT(4I8,1ES20.12)
 1004   FORMAT(1ES20.12)
 1005   FORMAT(5ES20.12)
 1006   FORMAT(3ES20.12)
 1007   FORMAT(4I8)
 
 1011   FORMAT(3I8)
 1012   FORMAT(6ES20.12)
 1013   FORMAT(2I8)
 1014   FORMAT(8ES20.12)
 1015   FORMAT(8I8)
 1016   FORMAT(3I8,1ES20.12)

      NV=NV+1

      RETURN
      END

!C*******************************************************************
      SUBROUTINE WRITEHISTORY(CFLM,DVMAX,QMMAX)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR

      IMPLICIT NONE

      REAL*8 CFLM,DVMAX,QMMAX

      WRITE(13,130) TIME,DT,CFLM,DVMAX,QMMAX
 130  FORMAT(F13.5,4ES15.7)

      IF (ICH.NE.0) THEN
      WRITE(96,140) TIME,PMI_DUDY,PMI
 140  FORMAT(F13.5,3ES20.12)
      ENDIF

      RETURN
      END
      
!C*******************************************************************
      SUBROUTINE WRITEHISTORY_HEAT(W,PSI_C)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR

      USE TWO_PHASE_PROPERTY
      USE HEAT_VAR

      IMPLICIT NONE

      REAL*8 W(0:M1,0:M2,0:M3)
      REAL*8 PSI_C(0:M1,0:M2,0:M3)

      INTEGER*8 I,J,K
      REAL*8 AA,UU,FUNCBODY,DVOL,UU_TMP,DENSCF,AMEAN_TEMP,ANUSSET
      
      AA=0.
      UU=0.
!$OMP PARALLEL DO private(I,J,DVOL,UU_TMP,DENSCF)
!$OMP&reduction(+:AA,UU)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
         DVOL=SDX(I)*SDY(J)*SDZ(K)

         UU_TMP=0.5*(W(I,J,KPV(K))+W(I,J,K))
         DENSCF=DENSCM+DENSC_DIFF*PSI_C(I,J,K)

         AA=AA+DENSCF*UU_TMP*T(I,J,K)*DVOL
         UU=UU+DENSCF*UU_TMP*DVOL
      ENDDO
      ENDDO
      ENDDO

       IF (AA .NE. 0. .AND. UU.NE. 0.) THEN
        AMEAN_TEMP=AA/UU
        ANUSSET=2./AMEAN_TEMP
        WRITE(*,*) 'THETA_MEAN=',AMEAN_TEMP,'NUSSET=',ANUSSET
       ELSE
        AMEAN_TEMP=0.
        ANUSSET=0.
       ENDIF
       
      WRITE(97,140)TIME,ANUSSET,AMEAN_TEMP
 140  FORMAT(F13.5,3ES20.12)
      

      RETURN
      END


!C*******************************************************************
      SUBROUTINE RK3COEF
!C*******************************************************************
      USE FLOW_VAR

      GAMMA(1)=8./15.
      GAMMA(2)=5./12.
      GAMMA(3)=3./4.
      RO(1)=0.
      RO(2)=-17./60.
      RO(3)=-5./12.

      RETURN
      END

!C******************************************************************
      SUBROUTINE CFL(CFLM,U,V,W)
!C******************************************************************
!C     CALCULATE THE MAXIMUM CFL NUMBER OF FLOW FIELD
      USE FLOW_VAR
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL*8 CFLM
      
      INTEGER*8 I,J,K
      REAL*8 CFLI1,CFLI2,CFLI3,CFLI

      CFLM=0d0

       IF ( N3M .EQ. 1 ) THEN
!$OMP PARALLEL DO private(I,CFLI1,CFLI2,CFLI)
!$OMP&reduction(MAX:CFLM)
      DO 12 J=1,N2M
      DO 12 I=1,N1M
        CFLI1=DABS(U(I,J,1)+U(IPV(I),J,1))*SSDX(I)
        CFLI2=DABS(V(I,J,1)+V(I,JPV(J),1))*SSDY(J)
        CFLI=0.5d0*(CFLI1+CFLI2)*DT
        CFLM=DMAX1(CFLI,CFLM)
 12    CONTINUE
          WRITE(*,*) '2D CFL IS CAL.'

        ELSE
!$OMP PARALLEL DO private(I,J,CFLI1,CFLI2,CFLI3,CFLI)
!$OMP&reduction(MAX:CFLM)
      DO 13 K=1,N3M
      DO 13 J=1,N2M
      DO 13 I=1,N1M
         CFLI1=DABS(U(I,J,K)+U(IPV(I),J,K))*SSDX(I)
         CFLI2=DABS(V(I,J,K)+V(I,JPV(J),K))*SSDY(J)
         CFLI3=DABS(W(I,J,K)+W(I,J,KPV(K)))*SSDZ(K)
         CFLI=0.5*(CFLI1+CFLI2+CFLI3)*DT
         CFLM=DMAX1(CFLI,CFLM)
!         IF (CFLI .GE. CFLM) THEN
!              write(*,*) i,j,k,cflm
!         ENDIF
   13 CONTINUE
       ENDIF

      RETURN
      END
       
!C------------------------------------------------
! RANDOM_NUMBER_GENERATOR

      FUNCTION RAN1(IDUM)
      
      INTEGER*8    IDUM,IA,IM,IQ,IR,NTAB,NDIV
      REAL*8       RAN1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     &           NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2E-7,RNMX=1.-EPS)
      INTEGER*8    J,K,IV(NTAB),IY
      
      SAVE IV,IY
      DATA IV /NTAB*0/, IY /0/
      IF (IDUM.LE.0 .OR. IY.EQ.0) THEN
        IDUM=MAX(-IDUM,1)
        DO 11 J=NTAB+8,1,-1
          K=IDUM/IQ
          IDUM=IA*(IDUM-K*IQ)-IR*K
          IF (IDUM.LT.0) IDUM=IDUM+IM
          IF (J.LE.NTAB) IV(J)=IDUM
11      CONTINUE
        IY=IV(1)
      ENDIF
      K=IDUM/IQ
      IDUM=IA*(IDUM-K*IQ)-IR*K
      IF (IDUM.LT.0) IDUM=IDUM+IM

      J=1+IY/NDIV
      IY=IV(J)
      IV(J)=IDUM
      RAN1=DMIN1(AM*IY,RNMX)

      RETURN
      END


!C******************************************************************
!C     SECTION FOR SOLVING UNCOUPLED MOMENTUM EQS. OF NAVIER-STOKES EQS.
!C     - RHSNLHS-EI.F
!C     - USES 2ND. ORDER CRANK-NICHOLSON METHOD FOR VISCOUS TERMS
!C     - USES 3RD. ORDER RUNGE-KUTTA METHOD FOR CONVECTIVE TERMS
!C     - ADI SCHEME IS USED TO SOLVE ELLIPTIC MOMENTUM EQS. IN Z->Y->X
!C       ORDER
!C     - CAN HANDLE OPEN( 2 TYPES ) OR CLOSED TOP ACCORDING TO ITOPEN.
!C     - FOR THESE SUB-ROUTINES, 'COMMON/WRK2~4/~ ' VARIABLES WORK
!C       AS TEMPORARY STORAGE. THEY STORE NUt AT EXCEPT CELL CENTER AND
!C       SHOULD NOT BE MODIFIED FOR WHOLE SGS & RHSNLHS PROCESS
!C     - TO REDUCE REQUIRED MEMORY FOR DNS, ERASE COMMON STATEMENTS
!C       HAVING TNU,TNUX,TNUY,TNUZ WHICH STORE SGS VISCOSITY FOR LES
!C       IN WHOLE RHSNLHS ROUTINES
!C     - MUST BE ACCOMPANIED BY LHS-EI.F OR LHS-EI-L.F(IN CASE OF LES)
!C
!C                                  LAST REVISED BY JUNGWOO KIM
!C                                            1999.08.04.   11:00
!c     -code for Interpolation scehme(bilinear interploation)
!C******************************************************************

!***********************************************************************
!       FULLY IMPLICIT & FULLY EXPLICIT DIRECT NUMERICAL SIMULATION CODE
!        -LARGE EDDY SIMULATION
!        -IMMERSED BOUNDARY METHOD
!        -LEVEL SET METHOD    
!
!       ORININAL INCOMPRESSIBLE NAVIER-STOKES EQUATION 
!                                IS CODED BY JUNIL LEE(SEPTEMBER 2007)
!
!                                                  2013.10. KIYOUNH KIM
!***********************************************************************
!C*******************************************************************
      SUBROUTINE RHSNLHS(ITRACKING,U,V,W,P,PSI_XN,PSI_YN,PSI_ZN
     &,PSI_CN,RK3XO,RK3YO,RK3ZO,ANPSO,QVOL_ORI,VOL1)
!*******************************************************************
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      
      USE HEAT_VAR
      
      USE TIME_VAR
      
      USE PHASE_BOILING
      
      IMPLICIT NONE

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 P(0:M1,0:M2,0:M3)

      REAL*8 PSI_XN(0:M1,0:M2,0:M3)
      REAL*8 PSI_YN(0:M1,0:M2,0:M3)
      REAL*8 PSI_ZN(0:M1,0:M2,0:M3)
      REAL*8 PSI_CN(0:M1,0:M2,0:M3)

      REAL*8 RK3XO(M1M,M2M,M3M),RK3YO(M1M,M2M,M3M),RK3ZO(M1M,M2M,M3M)
      REAL*8 ANPSO(M1M,M2M,M3M)

      REAL*8 DENF_X(0:M1,0:M2,0:M3)
      REAL*8 DENF_Y(0:M1,0:M2,0:M3)
      REAL*8 DENF_Z(0:M1,0:M2,0:M3)
      REAL*8 DENF_C(0:M1,0:M2,0:M3)
      
      REAL*8 QVOL_ORI
      REAL*8 VOL1
      
      REAL*8 DF_X(0:M1,0:M2,0:M3,3),VF_X(0:M1,0:M2,0:M3,3)
      REAL*8 DF_Y(0:M1,0:M2,0:M3,3),VF_Y(0:M1,0:M2,0:M3,3)
      REAL*8 DF_Z(0:M1,0:M2,0:M3,3),VF_Z(0:M1,0:M2,0:M3,3)
      REAL*8 DF_C(0:M1,0:M2,0:M3,3)

      REAL*8 RHS1(M1M,M2M,M3M,3)
      REAL*8 DENF_XI(0:M1,0:M2,0:M3),DENF_YI(0:M1,0:M2,0:M3)
     &                                         ,DENF_ZI(0:M1,0:M2,0:M3) 
     
      REAL*8 PSI_X(0:M1,0:M2,0:M3),PSI_Y(0:M1,0:M2,0:M3)
     &                     ,PSI_Z(0:M1,0:M2,0:M3),PSI_C(0:M1,0:M2,0:M3)

      INTEGER*8 I,J,K,ITRACKING
      REAL*8 QM1

      REAL*8 WRITE_ON
      REAL*8 RHS11,RHS22,RHS33,UU,VV,WW
      REAL*8 DENSC,DENSCF

!$OMP PARALLEL DO private(i,j)
       do k=1,n3m
       do j=1,n2m
       do i=1,n1m
       rhs1(i,j,k,1)=0d0
       rhs1(i,j,k,2)=0d0
       rhs1(i,j,k,3)=0d0
       enddo
       enddo
       enddo
       
      write_on=0

!      IF(ICH.EQ.1)CALL QVOLCALC(QVOL1,1,W,PSI_ZN,DENF_Z,VOL1)
      CALL REAL_TIME(CONTINUITY_B(MSUB))

      CALL CAL_CONTINUITY(DF_X,DF_Y,DF_Z,DF_C,VF_X,VF_Y,VF_Z
     &,U,V,W,PSI_XN,PSI_YN,PSI_ZN,PSI_CN,PSI_X,PSI_Y,PSI_Z,PSI_C)

      CALL CAL_CONT_BC(PSI_X,PSI_Y,PSI_Z,PSI_C,
     &     DENF_X,DENF_XI,DENF_Y,DENF_YI,DENF_Z,DENF_ZI,DENF_C)

      CALL REAL_TIME(CONTINUITY_E(MSUB))
      
      !ENERGY_EQUATION      
      IF (IPS .EQ. 1 ) THEN 
      CALL REAL_TIME(temp_B(MSUB))            
      CALL PASSIVE_SCALAR   
      CALL REAL_TIME(temp_E(MSUB))     
      ENDIF
      
!C-----MOMENTUM EQUATION SOLVER
      !CAL RHS
      CALL REAL_TIME(ANSETIME_B(MSUB))
      CALL RHSZ(RHS1,U,V,W,P,DF_Z,VF_Z,PSI_ZN,DENF_Z,RK3ZO)
!!=====SAVE==============================================================    
        OPEN(119,FILE='0VEL_rhs1.DAT') 
       WRITE(119,*) 'VARIABLES="X","Y","Z","rhs1"'
       WRITE(119,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      WRITE(119,117)XP(I),YP(J),Z(K),rhs1(i,j,k,3)
      ENDDO
      ENDDO
      ENDDO
       CLOSE(119)
 117  FORMAT(3F16.8,1D20.8)
      ! if (m .eq. 2) stop
!!=======================================================================    
       IF (ITRACKING .NE. 0) DEALLOCATE(SUR_X,SUR_Y,SUR_Z)
       IF (IPHS      .EQ. 1) DEALLOCATE(SGAS_X,SGAS_Y,SGAS_Z)
             
          I=1
          J=1
          
!$OMP PARALLEL DO 
         DO K=KBG,N3M
          RHS1(I,J,K,3)=RHS1(I,J,K,3)*DENF_ZI(I,J,K)-W(I,J,K)         
         ENDDO
!!=====SAVE==============================================================    
!        OPEN(119,FILE='0VEL_rhs1.DAT') 
!       WRITE(119,*) 'VARIABLES="X","Y","Z","rhs1"'
!       WRITE(119,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!      WRITE(119,117)XP(I),YP(J),Z(K),rhs1(i,j,k,3)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(119)
! 117  FORMAT(3F16.8,1D20.8)
!      stop
!!=======================================================================    
      !CAL LHS
        CALL LHSW(RHS1,U,V,W,DF_Z,VF_Z,DENF_ZI)

    
      !!BOUNDARY UPDATE BEFORE FSM_CHOI FOR 2ND-ORDER ACCURACY.
        CALL BC(RHS1,U,V,W)
        CALL FSM_CHOI(DENF_XI,DENF_YI,DENF_ZI,U,V,W,P,PSI_ZN,DENF_Z
     &,QVOL_ORI,VOL1)

      CALL REAL_TIME(ANSETIME_E(MSUB))
      CALL REAL_TIME(POISSTIME_B(MSUB))
  

      CALL PRESSURE_CORRECTION(U,V,W,P,
     &  DENF_Y,DENF_C,DENF_XI,DENF_YI,DENF_ZI)
 
       ! DEALLOCATED HERE FOR REQUIREMENT IN DIVGS
       !IF (IPHS .EQ. 1) DEALLOCATE(SFGAS) 
       ! if (iphs.eq. 1) DEALLOCATE(MFDOT)! SFGAS 는 DIVGS 까지 구하고 DEALLOCATE.  
       
      CALL real_time(POISSTIME_E(MSUB))
      

      RETURN
      END
      
!C******************************************************************
      SUBROUTINE RHSZ(RHS1,U,V,W,P,DF_Z,VF_Z,PSI_ZN,DENF_Z
     &,RK3ZO)
!C******************************************************************
!C     CALCULATION OF RHS(NAVIER-STOKES)
      USE PARAM_VAR
      USE FLOW_VAR

      USE FLOW_GEOM_VAR

      USE TWO_PHASE_PROPERTY
      USE PHASE_BOILING
      IMPLICIT NONE

      REAL*8 RHS1(M1M,M2M,M3M,3)
      REAL*8 RK3ZO(M1M,M2M,M3M)

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 P(0:M1,0:M2,0:M3)

      REAL*8 DENF_Z(0:M1,0:M2,0:M3)
      REAL*8 PSI_ZN(0:M1,0:M2,0:M3)
      REAL*8 DF_Z(0:M1,0:M2,0:M3,3),VF_Z(0:M1,0:M2,0:M3,3)


      INTEGER*8 I,J,K,IP1,IM1,JP1,JM1,KP1,KM1
      REAL*8 UE,UW,VN,VS,WE,WW,WN,WS,WT,WB
      REAL*8 ANZX,ANZY,ANZZ,ANZ
      REAL*8 ALZ1,ALZ2,ALZ3,ALZ4,ALZ5,ALZ6,ALZ7,ALZ8,ALZ9,ALZ10
      REAL*8 ALZX,ALZY,ALZZ,TALZX,TALZY,TALZZ,ALZ_SUM1,ALZ_SUM2
      REAL*8 PGR_Z,RK3Z,CN2Z,DEN_ZN,BUOYANCY

      REAL*8 FUNCBODY

!C-----RHS1 calculation for w momentum -----------------
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1)
!$OMP&private(UE,UW,VN,VS,WE,WW,WN,WS,WT,WB)
!$OMP&private(ANZX,ANZY,ANZZ,ANZ,ALZ1,ALZ2,ALZ3,ALZ4,ALZ5,ALZ6,ALZ7)
!$OMP&private(ALZ8,ALZ9,ALZ10,ALZX,ALZY,ALZZ,TALZX,TALZY,TALZZ)
!$OMP&private(ALZ_SUM1,ALZ_SUM2,PGR_Z,RK3Z,CN2Z,DEN_ZN,BUOYANCY)
       DO 1000 K=KBG,N3M
         KP1=KPV(K)
         KM1=KMV(K)
       DO 1000 J=1,N2M
         JP1=JPV(J)
         JM1=JMV(J)
       DO 1000 I=1,N1M
         IP1=IPV(I)
         IM1=IMV(I)

      !CONVECTION
      UE=0.5*(SDZ(KM1)*U(IP1,J,K)+SDZ(K)*U(IP1,J,KM1))*VVDZ(K)
      UW=0.5*(SDZ(KM1)*U(I,J,K)+SDZ(K)*U(I,J,KM1))*VVDZ(K)

      VN=0.5*(SDZ(KM1)*V(I,JP1,K)+SDZ(K)*V(I,JP1,KM1))*VVDZ(K)
      VS=0.5*(SDZ(KM1)*V(I,J,K)+SDZ(K)*V(I,J,KM1))*VVDZ(K)

      WE=(0.5*(SDX(IP1)*W(I,J,K)+SDX(I)*W(IP1,J,K))*VVDX(IP1))
     &*(1.-FIXIU(I))+FIXIU(I)*W(N1,J,K)
      WW=(0.5*(SDX(I)*W(IM1,J,K)+SDX(IM1)*W(I,J,K))*VVDX(I))
     &*(1.-FIXIL(I))+FIXIL(I)*W(0,J,K)

      WN=(0.5*(SDY(JP1)*W(I,J,K)+SDY(J)*W(I,JP1,K))*VVDY(JP1))
     & *(1.-FIXJU(J))+FIXJU(J)*W(I,N2,K)
      WS=(0.5*(SDY(J)*W(I,JM1,K)+SDY(JM1)*W(I,J,K))*VVDY(J))
     & *(1.-FIXJL(J))+FIXJL(J)*W(I,0,K)

      WT=0.5*(W(I,J,KP1)+W(I,J,K))
      WB=0.5*(W(I,J,K)+W(I,J,KM1))

!      ELSE !2ND-ORDER CENTRAL DIFFERENCE
      ANZX=(DF_Z(I,J,K,1)*UE*WE-DF_Z(IM1,J,K,1)*UW*WW)*SSDX(I)
      ANZY=(DF_Z(I,J,K,2)*VN*WN-DF_Z(I,JM1,K,2)*VS*WS)*SSDY(J)
      ANZZ=(DF_Z(I,J,K,3)*WT**2-DF_Z(I,J,KM1,3)*WB**2)*VVDZ(K)
!      ENDIF
      ANZ=-ANZX-ANZY-ANZZ                     

      !DIFFUSION
      ALZ1=(W(IP1,J,K)-W(I,J,K))*VVDX(IP1)
      ALZ2=(W(I,J,K)-W(IM1,J,K))*VVDX(I)
      ALZ3=(W(I,JP1,K)-W(I,J,K))*VVDY(JP1)
      ALZ4=(W(I,J,K)-W(I,JM1,K))*VVDY(J)
      ALZ5=(W(I,J,KP1)-W(I,J,K))*SSDZ(K)
      ALZ6=(W(I,J,K)-W(I,J,KM1))*SSDZ(KM1)
      ALZX=(VF_Z(I,J,K,1)*ALZ1-VF_Z(IM1,J,K,1)*ALZ2)*SSDX(I)
      ALZY=(VF_Z(I,J,K,2)*ALZ3-VF_Z(I,JM1,K,2)*ALZ4)*SSDY(J)
      ALZZ=(VF_Z(I,J,K,3)*ALZ5-VF_Z(I,J,KM1,3)*ALZ6)*VVDZ(K)

      ALZ7=(U(IP1,J,K)-U(IP1,J,KM1))
      ALZ8=(U(I,J,K)-U(I,J,KM1))
      ALZ9=(V(I,JP1,K)-V(I,JP1,KM1))
      ALZ10=(V(I,J,K)-V(I,J,KM1))
      TALZX=(VF_Z(I,J,K,1)*ALZ7-VF_Z(IM1,J,K,1)*ALZ8)*SSDX(I)*VVDZ(K)
      TALZY=(VF_Z(I,J,K,2)*ALZ9-VF_Z(I,JM1,K,2)*ALZ10)*SSDY(J)*VVDZ(K)
      TALZZ=ALZZ

      ALZ_SUM1=ALZX+ALZY+ALZZ+TALZZ
      ALZ_SUM2=TALZX+TALZY

      !OTHER TERM
       PGR_Z=(P(I,J,K)-P(I,J,KM1))*VVDZ(K)

      !BUOYANCY
      DEN_ZN=(DENM+DEN_DIFF*PSI_ZN(I,J,K))
      BUOYANCY=0d0! 0.5*(DEN_ZN+DENF_Z(I,J,K) )*GRAVITY

      RK3Z=ANZ+ALZ_SUM2
      CN2Z=ALZ_SUM1
!      RK3Z=ANZ+ALZ_SUM1+ALZ_SUM2
!      CN2Z=0.


      RHS1(I,J,K,3)=DEN_ZN*W(I,J,K)+
     &       ( GAMMA(MSUB)*RK3Z+RO(MSUB)*RK3ZO(I,J,K)
     &         +2.*ALPHA_RK3*CN2Z
     &         +2.*ALPHA_RK3*(SUR_Z(I,J,K)-PGR_Z)
     &         -2.*ALPHA_RK3*BUOYANCY
     &         -2.*ALPHA_RK3*(PMI+PMI_GRAV) )*DT
      if (sur_z(i,j,k) .ne. 0d0) then 
      write(*,*) ')(*&^',sur_z(i,j,k)  
      endif
      
       IF (IPHS .EQ. 1)then 
        if (IMASK_PH(I,J,K) .eq. 1) then 
        if (dabs(psi_zn(i,j,k)-0.5d0) .ne. 0.5d0) then 
       RHS1(I,J,K,3) = RHS1(I,J,K,3) + GAMMA(MSUB)*DT*den_zn*w(i,j,k)*
     & ((UE-UW)*SSDX(I)+(VN-VS)*SSDY(J)+(WT-WB)*VVDZ(K))
!     & + 2d0*ALPHA_RK3*DT*SGAS_Z(I,J,K)*(DEN_DIFF/DENR)*den_zn*w(i,j,k)
!     &  + 2D0*ALPHA_RK3*DT*den_zn*w(i,j,k)*(wt-wb)*vvdz(k)
      if (sgas_z(i,j,k) .ne. 0d0) then 
      write(*,*) '!@#$%^&*',sgas_z(i,j,k), w(i,j,k)
!      write(*,*) '!@#$%^&*',den_zn*w(i,j,k)*(wt-wb)*vvdz(k)
      endif
       endif
       ENDIF            
       ENDIF            
       

      RK3ZO(I,J,K)=RK3Z

 1000  CONTINUE

      RETURN
      END


!C******************************************************************
      SUBROUTINE LHSW(RHS1,U,V,W,DF_Z,VF_Z,DENF_ZI)
!C******************************************************************
      USE PARAM_VAR
      
      USE FLOW_VAR
      USE FLOW_GEOM_VAR

      IMPLICIT NONE

      REAL*8 RHS1(M1M,M2M,M3M,3)
      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL*8 DF_Z(0:M1,0:M2,0:M3,3),VF_Z(0:M1,0:M2,0:M3,3)
      REAL*8 DENF_ZI(0:M1,0:M2,0:M3)

      real, dimension (:,:), allocatable :: AI,BI,CI,GI
      real, dimension (:,:), allocatable :: AJ,BJ,CJ,GJ
      real, dimension (:,:), allocatable :: AK,BK,CK,GK
      
      INTEGER*8 I,J,K
      REAL*8 DI1,DI2,DI3,DJ1,DJ2,DJ3,DK1,DK2,DK3
 
      I=1
      J=1
!C-----Z-DIRECTION
      allocate(AK(M1,M3),BK(M1,M3),CK(M1,M3),GK(M1,M3))
      DO 141 K=KBG,N3M

      DK1=2.*ACOEF*DENF_ZI(I,J,K)*VF_Z(I,J,KMV(K),3)*AKW(K)
      DK3=2.*ACOEF*DENF_ZI(I,J,K)*VF_Z(I,J,K,3)*CKW(K)
      DK2=-DK1-DK3

      AK(I,K)=DK1
      CK(I,K)=DK3
      BK(I,K)=1.+DK2
      GK(I,K)=RHS1(I,J,K,3)
  141 CONTINUE
  
      IF (IPZ .EQ. 1) THEN
       CALL TRDIAG3P(AK,BK,CK,GK,KBG,N3M,1,1)
      ELSE
       CALL TRDIAG3(AK,BK,CK,GK,GK,KBG,N3M,1,1)
      ENDIF

      DO 151 K=KBG,N3M
      W(I,J,K)=GK(I,K)+W(I,J,K)
  151 CONTINUE

       OPEN(121,FILE='0dw_lhs.DAT') 
       WRITE(121,*) 'VARIABLES= "Z","lhs_dw"'
       WRITE(121,*)'ZONE  K=',N3M,',F=POINT'
      DO K=1,N3M
      WRITE(121,117) ZP(K),GK(1,k)
      ENDDO
       CLOSE(121)
 117  FORMAT(1F16.8,1D20.8)  

      deallocate(AK,BK,CK,GK)
      !stop
      RETURN
      END

!C******************************************************************
      SUBROUTINE BC(RHS1,U,V,W)
!C******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL*8 RHS1(M1M,M2M,M3M,3)
      
      INTEGER*8 I,J,K
      REAL*8 UBC_AVG,Q_NP,AA
      

!$OMP PARALLEL DO private(J)
          DO K=1,N3M
          DO J=1,N2M
          U(1 ,J,K)=0d0
          U(N1,J,K)=0d0
          V(0 ,J,K)=V(1  ,J,K) 
          V(N1,J,K)=V(N1M,J,K)                     
          W(0 ,J,K)=W(1  ,J,K) 
          W(N1,J,K)=W(N1M,J,K)          
         ENDDO
         ENDDO      
!$OMP PARALLEL DO private(I)
          DO K=1,N3M
          DO I=1,N1M
          u(I,0 ,K)=u(I,1 ,K)
          u(I,N2,K)=u(I,N2M,K)              
          V(I,1,K)=0d0
          V(I,N2,K)=0d0
          W(I,0 ,K)=W(I,1 ,K)
          W(I,N2,K)=W(I,N2M,K)          
         ENDDO
         ENDDO      
!$OMP PARALLEL DO private(I)
          DO J=1,N2M
          DO I=1,N1M
          W(I,J,1)=0D0
          W(I,J,N3)=W(I,J,N3M)
         ENDDO
         ENDDO

        
       CALL BC_FLUX(Q_NP,AA,U,V,W) !CALCULATE MASS FLUX AT THE BOUNDARY AT (N+1)-STEP.
       IF (AA .NE. 0.) UBC_AVG=-Q_NP/AA

      !VELOCITY CORRECTION FOR MASS CONSERVATION.
        DO J=1,N2M
        DO I=1,N1M
         W(I,J,1)=0d0
         W(I,J,N3)=W(I,J,N3)+UBC_AVG
        ENDDO
        ENDDO

       
      !BOUNDARY CONDITION IS TREATED EXPLICITLY.
      ! DIRICHLET AT INLET, NEUMANN AT OUTLET
!$OMP PARALLEL DO private(I)
      DO J=1,N2M
      DO I=1,N1M
      RHS1(I,J,N3M,3)=RHS1(I,J,N3M,3)
     &           -ACOEF*CKW(N3M)*(W(I,J,N3M)-W(I,J,N3))  !CKW HAS '-'SIGN.     
      ENDDO
      ENDDO

      
      RETURN
      END

!******************************************************************
      SUBROUTINE BC_FLUX(QQ,AA,U,V,W)
!******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      
      IMPLICIT NONE

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL*8 QQ,AA
      
      INTEGER*8 I,J,K
      REAL*8 QQ1,QQ2,QQ3,Q_IBM

      !OUTWARD IS POSITIVE.
        QQ1=0.
        QQ2=0.
        QQ3=0.
        AA=0.
        AA=AA+2.*XL*YL
!$OMP PARALLEL DO private(I)
!$OMP&reduction(+:QQ3)
        DO J=1,N2M
        DO I=1,N1M
         QQ3=QQ3+(W(I,J,N3)-W(I,J,1))*SDX(I)*SDY(J)
        ENDDO
        ENDDO


      QQ=QQ3

      RETURN
      END


!C*******************************************************************
      SUBROUTINE LHSINIT
!C*******************************************************************
!C     CALCULATE COEFFICIENTS FOR LHS
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      USE FLOW_VAR
      
      IMPLICIT NONE
      
      INTEGER*8 IC,IM,IP,JC,JM,JP,KC,KM,KP

!$OMP PARALLEL DO private(IM)
      DO 10 IC=IBG,N1M
      IM=IMV(IC)
      AIU(IC)=-VVDX(IC)*SSDX(IM)                 !   :i-1
      CIU(IC)=-VVDX(IC)*SSDX(IC)                 !   :i+1
   10 CONTINUE

!$OMP PARALLEL DO private(IP)
      DO 20 IC=1,N1M
      IP=IPV(IC)
      AIVW(IC)=-VVDX(IC)*SSDX(IC)
      CIVW(IC)=-VVDX(IP)*SSDX(IC)
   20 CONTINUE

!$OMP PARALLEL DO private(JM)
      DO 30 JC=JBG,N2M
      JM=JMV(JC)
      AJV(JC)=-VVDY(JC)*SSDY(JM)
      CJV(JC)=-VVDY(JC)*SSDY(JC)
   30 CONTINUE

!$OMP PARALLEL DO private(JP)
      DO 40 JC=1,N2M
      JP=JPV(JC)
      AJUW(JC)=-VVDY(JC)*SSDY(JC)
      CJUW(JC)=-VVDY(JP)*SSDY(JC)
   40 CONTINUE

!$OMP PARALLEL DO private(KM)
      DO 50 KC=KBG,N3M
      KM=KMV(KC)
      AKW(KC)=-VVDZ(KC)*SSDZ(KM)
      CKW(KC)=-VVDZ(KC)*SSDZ(KC)
   50 CONTINUE

!$OMP PARALLEL DO private(KP)
      DO 60 KC=1,N3M
      KP=KPV(KC)
      AKUV(KC)=-VVDZ(KC)*SSDZ(KC)
      CKUV(KC)=-VVDZ(KP)*SSDZ(KC)
   60 CONTINUE

      RETURN
      END
      
!C=======================CONTINUITY EQUATION============================C
!C      CAL. TWO-PHASE CONTINUITY EQUATION                              C
!C                                                                      C
!C      FLUX IS CALCULATED AT REFINED LEVEL-SET GRID                    C
!C      MASS-REDISTRIBUTION ALGORITHM IS IMPLEMETAION                   C
!C      (D.KIM STANFORD THESIS 2011)                                    C
!C                                                                      C
!C                                            KIYOUNG KIM 2015.3.18     C
!C=======================CONTINUITY EQUATION============================C
!C =====================================================================
      SUBROUTINE CAL_CONTINUITY(DF_X,DF_Y,DF_Z,DF_C
     &,VF_X,VF_Y,VF_Z,U,V,W,PSI_XN,PSI_YN,PSI_ZN,PSI_CN
     &,PSI_X,PSI_Y,PSI_Z,PSI_C)
!C =====================================================================
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      ! USE LES_VAR
      
      IMPLICIT NONE
      
      real*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      real*8 PSI_XN(0:M1,0:M2,0:M3)
      real*8 PSI_YN(0:M1,0:M2,0:M3)
      real*8 PSI_ZN(0:M1,0:M2,0:M3)
      real*8 PSI_CN(0:M1,0:M2,0:M3)
      
      real*8 DF_X(0:M1,0:M2,0:M3,3),VF_X(0:M1,0:M2,0:M3,3)
      real*8 DF_Y(0:M1,0:M2,0:M3,3),VF_Y(0:M1,0:M2,0:M3,3)
      real*8 DF_Z(0:M1,0:M2,0:M3,3),VF_Z(0:M1,0:M2,0:M3,3)
      real*8 DF_C(0:M1,0:M2,0:M3,3)
      integer*8 REGION(-1:M1+1,-1:M2+1,-1:M3+1)

      real*8 PSI_X(0:M1,0:M2,0:M3),PSI_Y(0:M1,0:M2,0:M3),
     &       PSI_Z(0:M1,0:M2,0:M3),PSI_C(0:M1,0:M2,0:M3)
     
      integer*8 ITRACKING
      
      integer*8 I,J,K,II,JJ,KK,I2,J2,K2
      real*8 DT_CON
       
!$OMP PARALLEL DO private(I,J)
      DO K=-1,N3+1
      DO J=-1,N2+1
      DO I=-1,N1+1
      REGION(I,J,K)=0
      ENDDO
      ENDDO
      ENDDO

!C-----DEVIDE REGION WITH GRADIENT
      i=1
      j=1
!$OMP PARALLEL DO private(II,JJ,KK,I2,J2,K2)
       DO K=1,N3M
        IF (   ( PSI_CN(I,J,KPV(K))-PSI_CN(I,J,K) .NE. 0.)
     &    .OR. ( PSI_CN(I,J,K)-PSI_CN(I,J,KMV(K)) .NE. 0.) ) THEN
         DO K2=K-2,K+2
            KK=K2
             IF (IPZ .EQ. 1) THEN
              IF (KK .LE. 0) KK=KK+N3M
              IF (KK .GE. N3) KK=KK-N3M
              ENDIF
          REGION(I,J,KK)=2
         ENDDO
        ELSE IF (PSI_CN(I,J,K) .EQ. 1.) THEN 
          IF (REGION(I,J,K) .NE. 2)   REGION(I,J,K)=1
       ENDIF
      ENDDO   

!C-----CONTINUITY SOLVER

      DT_CON=0.05*DTCONST
      CALL CONTI_Z(REGION,DF_Z,VF_Z,PSI_Z,DT_CON,U,V,W,PSI_ZN)
      CALL CONTI_C(REGION,DF_C,PSI_C,DT_CON,U,V,W
     &,PSI_XN,PSI_YN,PSI_ZN,PSI_CN) !FOR MULTIGRID AT POISSON


         
      RETURN
      END

!C =====================================================================
      SUBROUTINE CONTI_Z(REGION,DF_Z,VF_Z,PSI_Z,DT_CON,U,V,W
     &,PSI_ZN)
!C =====================================================================
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      real*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      real*8 PSI_ZN(0:M1,0:M2,0:M3)

      real*8 DFX(0:M1,0:M2,0:M3),DFY(0:M1,0:M2,0:M3),DFZ(0:M1,0:M2,0:M3)
      real*8 EPS(0:M1,0:M2,0:M3),EPS_OLD(0:M1,0:M2,0:M3)
      integer*8 PSI_ITERATION,STEADY_EPS
      
      integer*8 REGION(-1:M1+1,-1:M2+1,-1:M3+1)
      real*8 DF_Z(0:M1,0:M2,0:M3,3),VF_Z(0:M1,0:M2,0:M3,3) 
      real*8 PSI_Z(0:M1,0:M2,0:M3)
      real*8 TAUC(3,2)

      integer*8 ITRACKING
      real*8 DT_CON
      
      integer*8 I,J,K,IP1,IM1,JP1,JM1,KP1,KM1
      real*8 UE,UW,VN,VS,WT,WB
      real*8 PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_ZN,DEN_Z
      
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
        DF_Z(I,J,K,1)=0.
        DF_Z(I,J,K,2)=0.
        DF_Z(I,J,K,3)=0.
       ENDDO
       ENDDO
       ENDDO

      !DEFINE PSI_Z
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1,UE,UW,VN,VS,WT,WB)
!$OMP&private(PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_ZN,DEN_Z)
       DO K=KBG,N3M
        KP1=KPV(K)
        KM1=KMV(K)
       J=1
       I=1

      IF ( REGION(I,J,K) .EQ. 0 ) THEN
       VF_Z(I,J,K,3)=VISM
       VF_Z(I,J,KM1,3)=VISM
       DF_Z(I,J,K,3)=DENM
       DF_Z(I,J,KM1,3)=DENM

       PSI_Z(I,J,K)=PSI_ZN(I,J,K)

      ELSE IF (REGION(I,J,K) .EQ. 1 ) THEN
       VF_Z(I,J,K,3)=VISP
       VF_Z(I,J,KM1,3)=VISP
       DF_Z(I,J,K,3)=DENP
       DF_Z(I,J,KM1,3)=DENP

       PSI_Z(I,J,K)=PSI_ZN(I,J,K)

      ELSE

      WT=0.5*(W(I,J,KP1)+W(I,J,K))
      WB=0.5*(W(I,J,K)+W(I,J,KM1))


      IF ( DF_Z(I,J,K,3) .EQ. 0. ) THEN
       IF (DABS(WT) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_ZN(I,J,K)+PSI_ZN(I,J,KP1))
       ELSE
        CALL PSIFACE_Z(I,J,K,WT,PSI_SUM,3)
       ENDIF
        VF_Z(I,J,K,3)=VISM+(VISP-VISM)*PSI_SUM
        DF_Z(I,J,K,3)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_Z(I,J,KM1,3) .EQ. 0. ) THEN
       IF (DABS(WB) .LE. 1.E-8) THEN
         PSI_SUM=0.5*(PSI_ZN(I,J,K)+PSI_ZN(I,J,KM1))
       ELSE
        CALL PSIFACE_Z(I,J,KM1,WB,PSI_SUM,3)        !K-1 NOT KMV(K)
       ENDIF
        VF_Z(I,J,KM1,3)=VISM+(VISP-VISM)*PSI_SUM
        DF_Z(I,J,KM1,3)=DENM+DEN_DIFF*PSI_SUM
      ENDIF

      FLUX_Z=( DF_Z(I,J,K,3)*WT-DF_Z(I,J,KM1,3)*WB )*VVDZ(K)

       DEN_ZN=DENM+DEN_DIFF*PSI_ZN(I,J,K)
       DEN_Z=DEN_ZN-DTCONST*FLUX_Z
       PSI_Z(I,J,K)=(DEN_Z-DENM)/DEN_DIFF

        ENDIF
       ENDDO


!$OMP PARALLEL DO private(I,J)
       DO K=KBG,N3M
       DO J=1,N2M
       DO I=1,N1M
         IF (PSI_Z(I,J,K) .GT. 1.) THEN
           PSI_Z(I,J,K)=1.
         ELSE IF (PSI_Z(I,J,K) .LT. 0.) THEN
           PSI_Z(I,J,K)=0.
         ENDIF
       ENDDO
       ENDDO
       ENDDO

C-----CONTINUITY BOUNDARY CONDITION : NEUMANN+PERIODIC CONDITION.
!$OMP PARALLEL DO private(J)
      DO K=KBG,N3M
      DO J=1,N2M
       PSI_Z(0,J,K)=PSI_Z(1,J,K)
       PSI_Z(N1,J,K)=PSI_Z(N1M,J,K)
      ENDDO
      ENDDO

!$OMP PARALLEL DO private(I)
      DO K=KBG,N3M
      DO I=1,N1M
       PSI_Z(I,0,K)=PSI_Z(I,1,K)
       PSI_Z(I,N2,K)=PSI_Z(I,N2M,K)
      ENDDO
      ENDDO

!$OMP PARALLEL DO private(I)
       DO J=1,N2M
       DO I=1,N1M
       PSI_Z(I,J,1)=PSI_Z(I,J,2)
       PSI_Z(I,J,N3)=PSI_Z(I,J,N3M)
      ENDDO
      ENDDO
C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.


       RETURN
       END

!C =====================================================================
      SUBROUTINE CONTI_C(REGION,DF_C,PSI_C,DT_CON,U,V,W
     &,PSI_XN,PSI_YN,PSI_ZN,PSI_CN)
!C =====================================================================
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      real*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      real*8 PSI_XN(0:M1,0:M2,0:M3)
      real*8 PSI_YN(0:M1,0:M2,0:M3)
      real*8 PSI_ZN(0:M1,0:M2,0:M3)
      real*8 PSI_CN(0:M1,0:M2,0:M3)

      integer*8 REGION(-1:M1+1,-1:M2+1,-1:M3+1)
      real*8 DF_C(0:M1,0:M2,0:M3,3)
      real*8 PSI_C(0:M1,0:M2,0:M3)

      integer*8 ITRACKING
      real*8 DT_CON
      
      integer*8 I,J,K,IP1,IM1,JP1,JM1,KP1,KM1
      real*8 UE,UW,VN,VS,WT,WB
      real*8 PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_CN,DEN_C
    
!$OMP PARALLEL DO private(I,J)
       DO K=0,N3
       DO J=0,N2
       DO I=0,N1
        DF_C(I,J,K,1)=0D0
        DF_C(I,J,K,2)=0D0
        DF_C(I,J,K,3)=0D0
       ENDDO
       ENDDO
       ENDDO

      !DEFINE PSI_C
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1,UE,UW,VN,VS,WT,WB)
!$OMP&private(PSI_SUM,FLUX_X,FLUX_Y,FLUX_Z,DEN_CN,DEN_C)
       DO K=1,N3M
        KP1=KPV(K)
        KM1=KMV(K)
       J=1
       I=1

      IF ( REGION(I,J,K) .EQ. 0 ) THEN     
       PSI_C(I,J,K)=PSI_CN(I,J,K)
       DF_C(I,J,K,3)=DENM
       DF_C(I,J,KM1,3)=DENM

      ELSE IF (REGION(I,J,K) .EQ. 1 ) THEN
       PSI_C(I,J,K)=PSI_CN(I,J,K)
       DF_C(I,J,K,3)=DENP
       DF_C(I,J,KM1,3)=DENP

      ELSE

      WT=W(I,J,KP1)
      WB=W(I,J,K)

      IF ( DF_C(I,J,K,3) .EQ. 0. ) THEN
       IF (DABS(WT) .LE. 1.E-8) THEN
         PSI_SUM=PSI_CN(I,J,KP1)
       ELSE
        CALL PSIFACE_C(I,J,K,WT,PSI_SUM,3)
       ENDIF
        DF_C(I,J,K,3)=DENM+DEN_DIFF*PSI_SUM
      ENDIF
      IF ( DF_C(I,J,KM1,3) .EQ. 0. ) THEN
       IF (DABS(WB) .LE. 1.E-8) THEN
         PSI_SUM=PSI_CN(I,J,K)
       ELSE
        CALL PSIFACE_C(I,J,KM1,WB,PSI_SUM,3)        !K-1 NOT KMV(K)
       ENDIF
        DF_C(I,J,KM1,3)=DENM+DEN_DIFF*PSI_SUM
      ENDIF

       FLUX_Z=( DF_C(I,J,K,3)*WT-DF_C(I,J,KM1,3)*WB )*SSDZ(K)

       DEN_CN=DENM+DEN_DIFF*PSI_CN(I,J,K)
       DEN_C=DEN_CN-DTCONST*( FLUX_Z )
       PSI_C(I,J,K)=(DEN_C-DENM)/DEN_DIFF

        ENDIF
       ENDDO


!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
         IF (PSI_C(I,J,K) .GT. 1D0) THEN
           PSI_C(I,J,K)=1D0
         ELSE IF (PSI_C(I,J,K) .LT. 0D0) THEN
           PSI_C(I,J,K)=0D0
         ENDIF
       ENDDO
       ENDDO
       ENDDO

!C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
!$OMP PARALLEL DO private(J)
      DO K=1,N3M
      DO J=1,N2M
       PSI_C(0,J,K)=PSI_C(1,J,K)
       PSI_C(N1,J,K)=PSI_C(N1M,J,K)
      ENDDO
      ENDDO
!$OMP PARALLEL DO private(I)
      DO K=1,N3M
      DO I=1,N1M
       PSI_C(I,0,K)=PSI_C(I,1,K)
       PSI_C(I,N2,K)=PSI_C(I,N2M,K)
      ENDDO
      ENDDO
!$OMP PARALLEL DO private(I)
       DO J=1,N2M
       DO I=1,N1M
       PSI_C(I,J,1)=PSI_C(I,J,2)
       PSI_C(I,J,N3)=PSI_C(I,J,N3M)
      ENDDO
      ENDDO

       RETURN
       END


!C=======================================================================
      SUBROUTINE PSIFACE_Z(I,J,K,UU,PSI_SUM,IDIR)
!C=======================================================================  
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE FLOW_VAR

      USE LVS_COUPLING
      
      IMPLICIT NONE
      
      INTEGER*8 I,J,K,IDIR
      REAL*8 UU,PSI_SUM
      
      INTEGER*8 I1,I2,J1,J2,K1,K2
      REAL*8 X1,X2,Y1,Y2,Z1,Z2

      !X-DIR
       IF ( IDIR .EQ. 1) THEN  
        IF ( UU .GE. 0. ) THEN
       X1=X(I+1)-UU*DTCONST
       X2=X(I+1)
!            IF ( X1 .LT. X(1) ) X1=X(1)
        I1=ICOU2(I)
        I2=ICOU1(I+1)
        ELSE
       X1=X(I+1)
       X2=X(I+1)-UU*DTCONST
!            IF ( X2 .GT. X(N1) ) X2=X(N1)
        I1=ICOU2(I+1)
        I2=ICOU1(I+2)
        ENDIF
        Y1=Y(J)
        Y2=Y(J+1)
        Z1=ZP(K-1)
        Z2=ZP(K)

        J1=JCOU2(J)
        J2=JCOU1(J+1)
        K1=KCOUMP2(K-1)
        K2=KCOUMP1(K)

       !Y-DIR
       ELSE IF ( IDIR .EQ. 2) THEN

        X1=X(I)
        X2=X(I+1)
        IF ( UU .GE. 0. ) THEN
       Y1=Y(J+1)-UU*DTCONST
       Y2=Y(J+1)
!            IF ( Y1 .LT. Y(1) ) Y1=Y(1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)
        ELSE
       Y1=Y(J+1)
       Y2=Y(J+1)-UU*DTCONST
!            IF ( Y2 .GT. Y(N2) ) Y2=Y(N2)
        J1=JCOU2(J+1)
        J2=JCOU1(J+2)
        ENDIF
        Z1=ZP(K-1)
        Z2=ZP(K)

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        K1=KCOUMP2(K-1)
        K2=KCOUMP1(K)

      !Z-DIR
      ELSE IF ( IDIR .EQ. 3) THEN

        X1=X(I)
        X2=X(I+1)
        Y1=Y(J)
        Y2=Y(J+1)
         IF ( UU .GE. 0. ) THEN
       Z1=ZP(K)-UU*DTCONST
       Z2=ZP(K)
!            IF ( Z1 .LT. Z(1) ) Z1=Z(1)
        K1=KCOUMP2(K-1)
        K2=KCOUMP1(K)
         ELSE
       Z1=ZP(K)
       Z2=ZP(K)-UU*DTCONST
!            IF ( Z2 .GT. Z(N3) ) Z2=Z(N3)
        K1=KCOUMP2(K)
        K2=KCOUMP1(K+1)
         ENDIF

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)

      ENDIF

       CALL INTER_PSI_FACE(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)

        RETURN
        END

!C=======================================================================
      SUBROUTINE PSIFACE_C(I,J,K,UU,PSI_SUM,IDIR)
!C=======================================================================  
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE FLOW_VAR

      USE LVS_COUPLING
      
      IMPLICIT NONE
      
      INTEGER*8 I,J,K,IDIR
      REAL*8 UU,PSI_SUM
      
      INTEGER*8 I1,I2,J1,J2,K1,K2
      REAL*8 X1,X2,Y1,Y2,Z1,Z2
      
       !X-DIR
       IF ( IDIR .EQ. 1) THEN
       IF ( UU .GE. 0.) THEN
        X1=X(I+1)-UU*DTCONST
        X2=X(I+1)
!              IF ( X1 .LT. X(1) ) X1=X(1)
        I1=ICOU2(I)      !IF CFL IS GREATER THAN 1, THIS MAY HAVE PROBLEM.
        I2=ICOU1(I+1)
       ELSE
        X1=X(I+1)
        X2=X(I+1)-UU*DTCONST
!              IF ( X2 .GT. X(N1) ) X2=X(N1)
        I1=ICOU2(I+1)      !IF CFL IS GREATER THAN 1, THIS MAY HAVE PROBLEM.
        I2=ICOU1(I+2)
       ENDIF
       Y1=Y(J)
       Y2=Y(J+1)
       Z1=Z(K)
       Z2=Z(K+1)

        J1=JCOU2(J)
        J2=JCOU1(J+1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

      !Y-DIR
      ELSE IF ( IDIR .EQ. 2) THEN

       X1=X(I)
       X2=X(I+1)
       IF ( UU .GE. 0.) THEN
       Y1=Y(J+1)-UU*DTCONST    
       Y2=Y(J+1)
!            IF ( Y1 .LT. Y(1) ) Y1=Y(1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)
       ELSE
       Y1=Y(J+1)
       Y2=Y(J+1)-UU*DTCONST
!            IF ( Y2 .GT. Y(N2) ) Y2=Y(N2)
        J1=JCOU2(J+1)
        J2=JCOU1(J+2)
       ENDIF
       Z1=Z(K)
       Z2=Z(K+1)

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)

       ELSE IF ( IDIR .EQ. 3) THEN

       X1=X(I)
       X2=X(I+1)
       Y1=Y(J)
       Y2=Y(J+1)
         IF ( UU .GE. 0. ) THEN
       Z1=Z(K+1)-UU*DTCONST
       Z2=Z(K+1)
!            IF ( Z1 .LT. Z(1) ) Z1=Z(1)
        K1=KCOU2(K)
        K2=KCOU1(K+1)
         ELSE
       Z1=Z(K+1)
       Z2=Z(K+1)-UU*DTCONST
!            IF ( Z2 .GT. Z(N3) ) Z2=Z(N3)
        K1=KCOU2(K+1)
        K2=KCOU1(K+2)
         ENDIF

        I1=ICOU2(I)
        I2=ICOU1(I+1)
        J1=JCOU2(J)
        J2=JCOU1(J+1)

      ENDIF

       CALL INTER_PSI_FACE(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,PSI_SUM)

       RETURN
       END

!C=======================================================================
      SUBROUTINE INTER_PSI_FACE(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2
     &,PSI_SUM)
!C=======================================================================
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE
      
      INTEGER*8 I1,I2,J1,J2,K1,K2
      REAL*8 X1,X2,Y1,Y2,Z1,Z2
      REAL*8 PSI_SUM
      
      INTEGER*8 II,JJ,KK
      REAL*8 VOLF
      REAL*8 AI_CUT,AJ_CUT

      !SHOULD BE DONE GEOMETRICALLY.

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
!C =====================================================================
      SUBROUTINE CAL_CONT_BC(PSI_X,PSI_Y,PSI_Z,PSI_C,
     & DENF_X,DENF_XI,DENF_Y,DENF_YI,DENF_Z,DENF_ZI,DENF_C)
!C =====================================================================
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      USE PHASE_BOILING
      
      IMPLICIT NONE
      INTEGER*8 I,J,K
      
      REAL*8    PSI_X(0:M1,0:M2,0:M3),PSI_Y(0:M1,0:M2,0:M3),
     &          PSI_Z(0:M1,0:M2,0:M3),PSI_C(0:M1,0:M2,0:M3)
     
      REAL*8    DENF_X(0:M1,0:M2,0:M3),DENF_XI(0:M1,0:M2,0:M3),
     &          DENF_Y(0:M1,0:M2,0:M3),DENF_YI(0:M1,0:M2,0:M3),
     &          DENF_Z(0:M1,0:M2,0:M3),DENF_ZI(0:M1,0:M2,0:M3),
     &          DENF_C(0:M1,0:M2,0:M3)

!$OMP PARALLEL DO
         DO K=1,N3M
         J=1
         I=1
         DENF_Z(I,J,K)=DENM+DEN_DIFF*PSI_Z(I,J,K)
         DENF_ZI(I,J,K)=1./DENF_Z(I,J,K)
         ENDDO

         DENF_Z(1,1,N3)= DENF_Z(1,1,N3M)
        DENF_ZI(1,1,N3)=DENF_ZI(1,1,N3M) !FOR VARIABLE POISSON EQN.
         DENF_Z(1,1,0)= DENF_Z(1,1,1)
        DENF_ZI(1,1,0)=DENF_ZI(1,1,1) !FOR VARIABLE POISSON EQN.
        
!$OMP PARALLEL DO
       DO K=1,N3M
        DENF_C(1,1,K)=DENM+DEN_DIFF*PSI_C(1,1,K)
       ENDDO      
 
      
      RETURN
      END

!C =====================================================================
      SUBROUTINE CAL_CONTINUITY2(PSI_CN,PSI_X,PSI_Y,PSI_Z,PSI_C)
!C =====================================================================
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      USE LVS_COUPLING
      USE PHASE_BOILING
      
      ! USE LES_VAR
      
      IMPLICIT NONE

      INTEGER*8 I,J,K,II,JJ,KK,I1,J1,K1,I2,J2,K2      
      INTEGER*8 REGION(-1:M1+1,-1:M2+1,-1:M3+1)

      REAL*8    PSI_X(0:M1,0:M2,0:M3),PSI_Y(0:M1,0:M2,0:M3),
     &          PSI_Z(0:M1,0:M2,0:M3),PSI_C(0:M1,0:M2,0:M3)
      REAL*8    DVOL
      REAL*8    PSI_CN(0:M1,0:M2,0:M3)
      REAL*8    DEN_X,DEN_Y,DEN_Z,DEN_C,X1,X2,Y1,Y2,Z1,Z2
      REAL*8    S_SUM
      

      
      IF ((DENM .NE. DENP) .OR. (VISM .NE. VISP) ) THEN
      
!$OMP PARALLEL DO private(I,J)
      DO K=-1,N3+1
      DO J=-1,N2+1
      DO I=-1,N1+1
      REGION(I,J,K)=0
      ENDDO
      ENDDO
      ENDDO

            J=1
            I=1
!C-----DEVIDE REGION WITH GRADIENT
!$OMP PARALLEL DO private(KK,K2)
       DO K=1,N3M
        IF (   ( PSI_CN(1,1,KPV(K))-PSI_CN(1,1,K) .NE. 0.)
     &    .OR. ( PSI_CN(1,1,K)-PSI_CN(1,1,KMV(K)) .NE. 0.) ) THEN
         DO K2=K-2,K+2
            KK=K2
             IF (IPZ .EQ. 1) THEN
              IF (KK .LE. 0) KK=KK+N3M
              IF (KK .GE. N3) KK=KK-N3M
             ENDIF         
          REGION(1,1,KK)=2
         ENDDO
        ELSE IF (PSI_CN(1,1,K) .EQ. 1.) THEN 
          IF (REGION(1,1,K) .NE. 2)   REGION(1,1,K)=1
       ENDIF
      ENDDO
      ENDIF !      IF ((DENM .NE. DENP) .OR. (VISM .NE. VISP) ) THEN

            J=1
            I=1      
     
       S_SUM=0D0
!$OMP PARALLEL DO private(I,J)
!$OMP&private(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,S_SUM)
       DO 203 K=KBG,N3M
       DO 203 J=1,N2M
       DO 203 I=1,N1M
      IF ( REGION(I,J,K) .EQ. 2 ) THEN           
         X1=X(I)
         X2=X(I+1)
         I1=ICOU2(I)
         I2=ICOU1(I+1)
         Y1=Y(J)
         Y2=Y(J+1)
         J1=JCOU2(J)
         J2=JCOU1(J+1)         
         Z1=ZP(K-1)
         Z2=ZP(K)
         K1=KCOU2(K)
         K2=KCOU1(K+1)
       CALL INTER_PSI_MTP(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,S_SUM
     & ,SFGAS2)
       PSI_Z(I,J,K)=PSI_Z(I,J,K)-S_SUM
       ENDIF
  203 CONTINUE

!$OMP PARALLEL DO private(I,J)
       DO K=KBG,N3M
       DO J=1,N2M
       DO I=1,N1M
         IF (PSI_Z(I,J,K) .GT. 1.) THEN
           PSI_Z(I,J,K)=1.
         ELSE IF (PSI_Z(I,J,K) .LT. 0.) THEN
           PSI_Z(I,J,K)=0.
         ENDIF
       ENDDO
       ENDDO
       ENDDO  
            J=1
            I=1       

      S_SUM=0D0       
!$OMP PARALLEL DO private(I,J)
!$OMP&private(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,S_SUM)
       DO 204 K=1,N3M
       DO 204 J=1,N2M
       DO 204 I=1,N1M
      IF ( REGION(I,J,K) .EQ. 2 ) THEN
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
       CALL INTER_PSI_MTP(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,S_SUM
     & ,SFGAS2)
       PSI_C(I,J,K)=PSI_C(I,J,K)-S_SUM   
        ENDIF
  204 CONTINUE

            J=1
            I=1
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3M
       DO J=1,N2M
       DO I=1,N1M
         IF (PSI_C(I,J,K) .GT. 1.) THEN
           PSI_C(I,J,K)=1.
         ELSE IF (PSI_C(I,J,K) .LT. 0.) THEN
           PSI_C(I,J,K)=0.
         ENDIF
       ENDDO
       ENDDO
       ENDDO     
     
            J=1
            I=1     

!C-----CONTINUITY BOUNDARY CONDITION : NEUMAN+PERIODIC CONDITION.
!$OMP PARALLEL DO private(I)
       DO J=1,N2M
       DO I=1,N1M
       PSI_Z(I,J,1)=PSI_Z(I,J,2)
       PSI_Z(I,J,N3)=PSI_Z(I,J,N3M)
      ENDDO
      ENDDO
!$OMP PARALLEL DO private(I)
       DO J=1,N2M
       DO I=1,N1M
       PSI_C(I,J,1)=PSI_C(I,J,2)
       PSI_C(I,J,N3)=PSI_C(I,J,N3M)
      ENDDO
      ENDDO
      
         
      RETURN
      END
      
!C===========================CORRECTION_STEP============================C
!C      TWO-PHASE PRESSURE-CORRECTION ALGORITHM                          C
!C                                                                      C
!C      CALCULATE CONTANT POISSON EQUATION ITERATIVELY INSTEAD OF        C
!C     SOLVING VARIABLE POISSON EQUATION ( D.KIM STANFORD THESIS 2011)  C
!C                                                                      C
!C                                            KIYOUNG KIM 2012.8.10     C
!C===========================CORRECTION_STEP============================C
!C*******************************************************************
      SUBROUTINE DIVGS(U,V,W,DIVGSUM)
!C*******************************************************************
!C     POISSON EQUATION PRE-PROCESSING
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      
      USE PHASE_BOILING
      USE LVS_COUPLING
      USE TWO_PHASE_PROPERTY
      
      IMPLICIT NONE

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 DIVGSUM(M1M,M2M,M3M)

      INTEGER*8 I,J,K,II,JJ,KK,IMM,JMM,KMM,N
      REAL*8 DIVG1,DIVG2,DIVG3
      REAL*8 S_SUM,X1,X2,Y1,Y2,Z1,Z2,dvol
      INTEGER*8 I1,I2,J1,J2,K1,K2
      
      divgsum=0d0
!$OMP PARALLEL DO private(DIVG3)
      DO K=1,N3M-1
       DIVG3=(W(1,1,KPV(K))-W(1,1,K))*SSDZ(K)
       DIVGSUM(1,1,K)=(DIVG3)
      ENDDO
      divgsum(1,1,n3m)=0d0
!!$OMP PARALLEL DO private(I,J,DIVG3)
!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!       DIVG3=(W(I,J,KPV(K))-W(I,J,K))*SSDZ(K)
!       DIVGSUM(I,J,K)=(DIVG3)
!      ENDDO
!      ENDDO
!      ENDDO
      
        OPEN(121,FILE='0VEL_DIVGS1.DAT') 
       !WRITE(121,*) 'VARIABLES="X","Y","Z","dvg0"'
       WRITE(121,*) 'VARIABLES= "Z","dvg0"'
       !WRITE(121,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
       WRITE(121,*)'ZONE  K=',N3M,',F=POINT'
      DO K=1,N3M
      !DO J=1,N2M
      !DO I=1,N1M
      !WRITE(121,117)XP(I),YP(J),ZP(K),DIVGSUM(i,j,k)
      WRITE(121,117) ZP(K),DIVGSUM(1,1,k)
      !ENDDO
      !ENDDO
      ENDDO
       CLOSE(121)
       

!$OMP PARALLEL DO private(I,J)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
       QMASS(I,J,K)=0.
      ENDDO
      ENDDO
      ENDDO

      
      IF (IPHS .EQ. 1) THEN 
!$OMP PARALLEL DO private(I,J,X1,X2,I1,I2,Y1,Y2,J1,J2,Z1,Z2,K1,K2,S_SUM)
!$OMP&private(DVOL)
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
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
       CALL INTER_PSI_MTP(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,S_SUM
     & ,sfgas)
       DVOL=SDX(I)*SDY(J)*SDZ(K)
       DIVGSUM(I,J,K)=DIVGSUM(I,J,K)-S_SUM*(DEN_DIFF/(DENP*denm))!*dvol
      ENDDO
      ENDDO
      ENDDO      
      ENDIF
      
!!=====SAVE==============================================================    
        OPEN(120,FILE='0VEL_DIVGS2.DAT') 
       !WRITE(120,*) 'VARIABLES="X","Y","Z","dvg-s"'
       WRITE(120,*) 'VARIABLES= "Z","dvg-s"'
       !WRITE(120,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
       WRITE(120,*)'ZONE  K=',N3M,',F=POINT'
      DO K=1,N3M
      !DO J=1,N2M
      !DO I=1,N1M
      !WRITE(120,117)XP(I),YP(J),ZP(K),DIVGSUM(i,j,k)
      WRITE(120,117) ZP(K),DIVGSUM(1,1,k)
      !ENDDO
      !ENDDO
      ENDDO
       CLOSE(120)
! 117  FORMAT(3F16.8,1D20.8)
 117  FORMAT(1F16.8,1D20.8)
            ! if (m .eq. 2) stop
!!=======================================================================           
      !stop
!$OMP PARALLEL DO private(I,J)
      DO 80 K=2,N3M
      DO 80 J=1,N2M
      DO 80 I=1,N1M
      DIVGSUM(I,J,K)=DIVGSUM(I,J,K)*DTCONSTI
   80 CONTINUE
 

      RETURN
      END

C*******************************************************************
      SUBROUTINE CONVRGE1(DVMAX,U,V,W)
C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      USE IBM_VAR
      
      use flow_var
      USE PHASE_BOILING
      USE LVS_COUPLING
      USE TWO_PHASE_PROPERTY      
      IMPLICIT NONE
      
      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL*8 DVMAX
      
      INTEGER*8 I,J,K
      REAL*8 DVG11,DVG12,DVG13,DIVGMAX_TMP
      
      REAL*8      S_SUM,X1,X2,Y1,Y2,Z1,Z2,dvol
      INTEGER*8   I1,I2,J1,J2,K1,K2
      
      DVMAX=0.
!$OMP PARALLEL DO private(I,J,DVG11,DVG12,DVG13,DIVGMAX_TMP)
!$OMP&private(X1,X2,I1,I2,Y1,Y2,J1,J2,Z1,Z2,K1,K2,S_SUM)
!$OMP&reduction(MAX:DVMAX)
      DO 13 K=1,N3M
      DO 13 J=1,N2M
      DO 13 I=1,N1M
       DVG13=(W(I,J,KPV(K))-W(I,J,K))*SSDZ(K)
       DIVGMAX_TMP=DVG13  !check qmass..
       if (iphs .eq. 1) then 
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
       CALL INTER_PSI_MTP(X1,X2,Y1,Y2,Z1,Z2,I1,I2,J1,J2,K1,K2,S_SUM
     & ,sfgas)
       DIVGMAX_TMP=DIVGMAX_TMP -S_SUM*(DEN_DIFF/DENP)      
       endif
       !WRITE(*,*) DIVGMAX_TMP       
       DVMAX=DMAX1(DABS(DIVGMAX_TMP),DVMAX)

   13 CONTINUE
   
       
       RETURN
       END
       
!C =====================================================================      
      SUBROUTINE FSM_CHOI(DENF_XI,DENF_YI,DENF_ZI,U,V,W,P,PSI_ZN,DENF_Z
     &,QVOL_ORI,VOL1)
!C ===================================================================== 
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 P(0:M1,0:M2,0:M3)

      REAL*8 PSI_ZN(0:M1,0:M2,0:M3)
      REAL*8 DENF_Z(0:M1,0:M2,0:M3)
 
      REAL*8 DENF_XI(0:M1,0:M2,0:M3),DENF_YI(0:M1,0:M2,0:M3)
     &                                         ,DENF_ZI(0:M1,0:M2,0:M3) 
      REAL*8 VOL1
      REAL*8 QVOL_ORI
      
      INTEGER*8 I,J,K,JM1,KM1
      REAL*8 QVOLH,PHCAP
      REAL*8 FUNCBODY

!$OMP PARALLEL DO private(I,J,KM1)
      DO 42 K=KBG,N3M
         KM1=KMV(K)
      DO 42 J=1,N2M
      DO 42 I=1,N1M
         W(I,J,K)=W(I,J,K)
     &   +DTCONST*( (P(I,J,K)-P(I,J,KM1))*VVDZ(K))*DENF_ZI(I,J,K)
   42  CONTINUE

       RETURN
       END

!C =====================================================================      
      SUBROUTINE PRESSURE_CORRECTION(U,V,W,P,DENF_Y,DENF_C
     &,DENF_XI,DENF_YI,DENF_ZI)
!C =====================================================================        
      USE PARAM_VAR
      
      USE FLOW_GEOM_VAR
      USE FLOW_VAR
      
      IMPLICIT NONE
      
      REAL*8 DENF_XI(0:M1,0:M2,0:M3),DENF_YI(0:M1,0:M2,0:M3)
     &                                         ,DENF_ZI(0:M1,0:M2,0:M3)

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 P(0:M1,0:M2,0:M3)
      REAL*8 DIVGSUM(M1M,M2M,M3M)
      REAL*8 DENF_Y(0:M1,0:M2,0:M3)
      REAL*8 DENF_C(0:M1,0:M2,0:M3)      
      REAL*8 FUNCBODY
      INTEGER*8 I,J,K
      REAL*8 PHIREF,DVOL,VOL
 
 
       CALL DIVGS(U,V,W,DIVGSUM)
 
       CALL POISSON(P,DIVGSUM,DENF_XI,DENF_YI,DENF_ZI,DENF_Y,DENF_C)

       CALL UCALC(DENF_XI,DENF_YI,DENF_ZI,U,V,W,P)

!C     SET THE AVERAGE PHI AT THE UPPER WALL TO BE ZERO.
      PHIREF=0.
      VOL=0.
!$OMP PARALLEL DO private(I,J,DVOL)
!$OMP&reduction(+:PHIREF,VOL)
      DO 80 K=1,N3M
      DO 80 J=1,N2M
      DO 80 I=1,N1M
        DVOL=SDX(I)*SDY(J)*SDZ(K)
        PHIREF=PHIREF+P(I,J,K)*DVOL
        VOL=VOL+DVOL
   80 CONTINUE
      PHIREF=P(1,1,1)!PHIREF/VOL

!$OMP PARALLEL DO private(I,J)
      DO 93 K=1,N3M
      DO 93 J=1,N2M
      DO 93 I=1,N1M
        P(I,J,K)=P(I,J,K)-PHIREF
   93 CONTINUE

!$OMP PARALLEL DO private(J)
      DO K=1,N3M
      DO J=1,N2M
         P(0,J,K)=P(1,J,K)
      ENDDO
      ENDDO
!$OMP PARALLEL DO private(I)
      DO K=1,N3M
      DO I=1,N1M
         P(I,0,K)=P(I,1,K)
      ENDDO
      ENDDO
!$OMP PARALLEL DO private(I)
      DO J=1,N2M
      DO I=1,N1M
         P(I,J,0) =P(I,J,1)
      ENDDO
      ENDDO

!!=====SAVE==============================================================    

        OPEN(118,FILE='0TDMA_P.DAT') 

       WRITE(118,*) 'VARIABLES="Z","P","w"'
       WRITE(118,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
      DO K=1,N3M
      J=1
      I=1
      WRITE(118,117) Z(K),P(I,J,K),w(i,j,k)
      ENDDO
       CLOSE(118)
 117  FORMAT(1F16.8,2D20.8)
       ! stop
!!=======================================================================            
       RETURN
       END
       
       
!C=======================================================================
      SUBROUTINE POISSON(P,DIVGSUM,DENF_XI,DENF_YI,DENF_ZI
     &,DENF_Y,DENF_C)
!C=======================================================================
      USE FLOW_VAR
      USE MG_OPTION

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE       
      INTEGER ICG_LIMIT,I,J,K,IP1,JP1,KP1,ITER
      REAL TEST1,ZRESID,ZRESID_OLD,ALPHA,RESMAX,BETA
      REAL DENF_Y(0:M1,0:M2,0:M3)
      REAL DENF_C(0:M1,0:M2,0:M3)
      
      REAL AW(M1MD,M2M,M3MD),AE(M1MD,M2M,M3MD)
     &    ,AS(M1MD,M2M,M3MD),AN(M1MD,M2M,M3MD)
     &    ,AB(M1MD,M2M,M3MD),AF(M1MD,M2M,M3MD),AC(M1MD,M2M,M3MD)

      REAL GAM(M3),BET
      REAL UU(M3,M1)

      REAL DIVGSUM(M1M,M2M,M3M)
      REAL P(0:M1,0:M2,0:M3)

      REAL DENF_XI(0:M1,0:M2,0:M3),DENF_YI(0:M1,0:M2,0:M3)
     &                                         ,DENF_ZI(0:M1,0:M2,0:M3)  

      real, dimension (:,:), allocatable :: AK,BK,CK,GK
      
      REAL*8 DI1,DI2,DI3,DJ1,DJ2,DJ3,DK1,DK2,DK3
      integer*8 kb
      real*8 phiref
      

       TEST1=RESID_POI/DTCONST
        RESMAX=0d0          
       J=1
       I=1

!C-----Z-DIRECTION
      allocate(AK(M3,M1),BK(M3,M1),CK(M3,M1),GK(M3,M1))

      DO 141 K=1,N3M
        KP1=KPV(K)
      AK(K,1)=VVDZ(K)*SSDZ(K)*DENF_ZI(1,1,K)*(1.-FIXKL(K))
      CK(K,1)=VVDZ(KP1)*SSDZ(K)*DENF_ZI(1,1,KP1)!*(1.-FIXKU(K))
      BK(K,1)=-( AK(K,1)+ CK(K,1) )
      GK(K,1)= divgsum(1,1,K)
  141 CONTINUE

      UU=0d0 

        kb=1!kbg
      !call CTRDIAG2(AK,BK,CK,GK,1,n3m,UU,1)
       GAM=0d0
       BET=1./BK(KB,1)
       UU(KB,1)=gk(KB,1)*BET
      
      I=1
      DO 20 k=kb+1,n3m
      GAM(k)=CK(k-1,1)*BET
      BET=1./(BK(k,1)-AK(k,1)*GAM(k))
      UU(k,1)=(gk(k,1)-AK(k,1)*UU(k-1,1))*BET
  20  CONTINUE
  
      DO 30 k=n3m-1,kb,-1
      UU(k,1)=UU(k,1)-GAM(k+1)*UU(k+1,1)
  30  CONTINUE
  
      DO 151 K=1,N3M
      P(1,1,K)=UU(K,1)
 151  CONTINUE
 

  
      deallocate(AK,BK,CK,GK)

!$OMP PARALLEL DO private(I)
      DO J=1,N2M
      DO I=1,N1M
         P(I,J,0) =P(I,J,1)
         P(I,J,n3) =P(I,J,n3m)
      ENDDO
      ENDDO
      return
      end
 
      
!C*******************************************************************
      SUBROUTINE UCALC(DENF_XI,DENF_YI,DENF_ZI,U,V,W,P)
!C*******************************************************************
!C     CALCULATING U FROM UHAT
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      IMPLICIT NONE

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 P(0:M1,0:M2,0:M3)

      REAL*8 DENF_XI(0:M1,0:M2,0:M3),DENF_YI(0:M1,0:M2,0:M3)
     &                                         ,DENF_ZI(0:M1,0:M2,0:M3)
     
      INTEGER*8 I,J,K,JM1,KM1
      REAL*8 FUNCBODY

!$OMP PARALLEL DO private(I,J,KM1)
      DO 42 K=KBG,N3M
         KM1=KMV(K)
      DO 42 J=1,N2M
      DO 42 I=1,N1M
       W(I,J,K)=W(I,J,K)
     &        -DTCONST*(P(I,J,K)-P(I,J,KM1))*VVDZ(K)*DENF_ZI(I,J,K)
   42  CONTINUE


      !When Neumann condition is used, p(1,j,k)-p(0,j,k)=0
      !boundary condition like U(1,J,K)=U(2,J,K) is not needed!!
      !U(1,J,K)=U(2,J,K) IS UADATED WITH INTERMEDIATE VELOCITY.
!$OMP PARALLEL DO private(I)
          DO J=1,N2M
          DO I=1,N1M
          U(I,J,0)=0D0
          V(I,J,0)=0D0
          w(i,j,1)=0d0
          U(I,J,N3)=U(I,J,N3M)
          V(I,J,N3)=V(I,J,N3M)
          w(i,j,n3)=w(i,j,n3m)
         ENDDO
         ENDDO

      RETURN
      END
      
!C*******************************************************************
      SUBROUTINE PASSIVE_SCALAR(U,V,W,PSI_CN,DF_C)
!C*******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE HEAT_VAR
      
      IMPLICIT NONE
      INTEGER*8   I,J,K
      REAL*8      U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)    
      REAL*8      PSI_CN(0:M1,0:M2,0:M3)  
      REAL*8      ANPSO(M1M,M2M,M3M),DENSC,DENSCF  
      REAL*8      DH(0:M1,0:M2,0:M3,3)
      !CALCULATE HX,HY,HZ FOR SHARP T CACULATION
      
      IF ((IPHS .EQ. 1) .AND. (M .NE. 1)) THEN 
         
         CALL SHARP_DIS(DH)
         CALL CAL_RHSPS2(U,V,W,ANPSO,PSI_CN,DF_C)
         
 !$OMP PARALLEL DO private(I,J,DENSC,DENSCF)
         DO K=1,N3M
         DO J=1,N2M
         DO I=1,N1M
              IF (PSI_CN(I,J,K) .GT. 0.5D0) THEN 
              DENSC=DENSCM+DENSC_DIFF
              DENSCF=DENSC
              ELSE
              DENSC=DENSCM
              DENSCF=DENSC
              ENDIF
           RHSPS(I,J,K)=( (DENSC-DENSCF)*T(I,J,K)+RHSPS(I,J,K))/DENSCF
         ENDDO
         ENDDO
         ENDDO      
         
        CALL LHSPS2(PSI_CN,PSI_C,DF_C)       

        CALL TEMP_GRAD(PSI_CN)                                     !PURPOSE: OBTAIN M_DOTS
        CALL CAL_CONTINUITY2(PSI_CN,PSI_X,PSI_Y,PSI_Z,PSI_C)       !PURPOSE: ADD M_DOT TO HATTED DEN.
        CALL CAL_CONT_BC(PSI_X,PSI_Y,PSI_Z,PSI_C,       
     &       DENF_X,DENF_XI,DENF_Y,DENF_YI,DENF_Z,DENF_ZI,DENF_C)        
       
      ELSE
         
         CALL CAL_RHSPS(U,V,W,ANPSO,PSI_CN,DF_C)        
 !$OMP PARALLEL DO private(I,J,DENSC,DENSCF)
         DO K=1,N3M
         DO J=1,N2M
         DO I=1,N1M
           DENSC=DENSCM+DENSC_DIFF*PSI_CN(I,J,K)
           DENSCF=DENSCM+DENSC_DIFF*PSI_C(I,J,K)               
           RHSPS(I,J,K)=( (DENSC-DENSCF)*T(I,J,K)+RHSPS(I,J,K))/DENSCF
         ENDDO
         ENDDO
         ENDDO                
         CALL LHSPS(PSI_CN,PSI_C,DF_C)         
         
      ENDIF

         
!=====SAVE==============================================================        
!         OPEN(138,FILE='0TEST_t.DAT')
!       WRITE(138,*) 'VARIABLES="X","Y","Z","t0"'
!       WRITE(138,*)'ZONE I=',N1M,',J=',N2M,',K=',N3M,',F=POINT'
!      DO K=1,N3M
!      DO J=1,N2M
!      DO I=1,N1M
!      WRITE(138,137)XP(I),YP(J),ZP(K),t(i,j,k)
!      ENDDO
!      ENDDO
!      ENDDO
!       CLOSE(138)
! 137  FORMAT(3F16.8,1D20.8)        
!======================================================================= 
            

      RETURN
      END
!C******************************************************************
!C     SECTION FOR SOLVING UNCOUPLED MOMENTUM EQS. TO GET UHAT
!C     - LHS-EI.F
!C     - FOR SPECIFIC FEATURES, SEE 'rhsnlhs-ei.f'
!C     - THIS IS FASTER THAN 'LHS FOR LES' BECAUSE THIS DOESN'T
!C       CALCULATE COEFS. EVERYTIME
!C     - CAN HANDLE OPEN( 2 TYPES ) OR CLOSED TOP ACCORDING TO ITOPEN.
!C
!C                                  LAST REVISED BY SEONGWON KANG
!C                                            1998.07.06.   13:00
!C******************************************************************

!C*******************************************************************
      SUBROUTINE TRDIAG1(A,B,C,R,UU,L1,L2,LL1,LL2)
!C*******************************************************************
!C     SOLVE THE TRIDIAGONAL MATRIX (X-1 TRIDIAGONAL) WITH
!C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      DIMENSION GAM(M2,M1),A(M2,M1),B(M2,M1),C(M2,M1),
     &          R(M2,M1),UU(M2,M1),BET(M2)

      DO 10 I=LL1,LL2
      BET(I)=1./B(I,L1)
      UU(I,L1)=R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)=C(I,J-1)*BET(I)
      BET(I)=1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J)=(R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE
      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J)=UU(I,J)-GAM(I,J+1)*UU(I,J+1)
   30 CONTINUE

      RETURN
      END

!C*******************************************************************
      SUBROUTINE TRDIAG2(A,B,C,R,UU,L1,L2,LL1,LL2)
!C*******************************************************************
!C     SOLVE THE TRIDIAGONAL MATRIX (X-2 TRIDIAGONAL) WITH
!C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      DIMENSION GAM(M1,M2),A(M1,M2),B(M1,M2),C(M1,M2),
     &          R(M1,M2),UU(M1,M2),BET(M1)

      DO 10 I=LL1,LL2
      BET(I)=1./B(I,L1)
      UU(I,L1)=R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)=C(I,J-1)*BET(I)
      BET(I)=1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J)=(R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE
      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J)=UU(I,J)-GAM(I,J+1)*UU(I,J+1)
   30 CONTINUE

      RETURN
      END

!C*******************************************************************
      SUBROUTINE TRDIAG3(A,B,C,R,UU,L1,L2,LL1,LL2)
!C*******************************************************************
!C     SOLVE THE TRIDIAGONAL MATRIX (X-3 TRIDIAGONAL) WITH
!C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      DIMENSION GAM(M1,M3),A(M1,M3),B(M1,M3),C(M1,M3),
     &          R(M1,M3),UU(M1,M3),BET(M1)

      DO 10 I=LL1,LL2
      BET(I)=1./B(I,L1)
      UU(I,L1)=R(I,L1)*BET(I)
   10 CONTINUE

      DO 20 J=L1+1,L2
      DO 20 I=LL1,LL2
      GAM(I,J)=C(I,J-1)*BET(I)
      BET(I)=1./(B(I,J)-A(I,J)*GAM(I,J))
      UU(I,J)=(R(I,J)-A(I,J)*UU(I,J-1))*BET(I)
   20 CONTINUE
      DO 30 J=L2-1,L1,-1
      DO 30 I=LL1,LL2
      UU(I,J)=UU(I,J)-GAM(I,J+1)*UU(I,J+1)
      R(I,J)=UU(I,J)
   30 CONTINUE

      RETURN
      END

!C*******************************************************************
      SUBROUTINE TRDIAG1P(A,B,C,F,J1,J2,L1,L2)
!C*******************************************************************
!C     INVERT A PERIODIC MATRIX (X-3 TRIDIAGONAL) WITH
!C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      REAL*8 A(M2,M1),B(M2,M1),C(M2,M1),F(M2,M1),Q(M2,M1),
     %     S(M2,M1),QE(M2,M1),FN(M2),PN(M2)

      JA=J1+1                                 
      JJ=J1+J2                                
      DO 10 K=L1,L2       !for L1 to L2 set (j1,2 is responsible for discretized dir)                   
      BINV=1./B(K,J1)     !calculate the inverse of the first line j1                    
      Q(K,J1)=-C(K,J1)*BINV  ! the first line C                 
      S(K,J1)=-A(K,J1)*BINV  ! the first line A                 
      FN(K)=F(K,J2)          ! the last line RHS1 related (not devided)                  
      F(K,J1)=F(K,J1)*BINV   ! the first line RHS1-related                  
   10 CONTINUE

!C     FORWARD ELIMINATION SWEEP
      DO 20 J=JA,J2          ! J1 is excluded: no worry about the division above
      DO 20 K=L1,L2
      PN(K)=1./(B(K,J)+A(K,J)*Q(K,J-1))
      Q(K,J)=-C(K,J)*PN(K)
      S(K,J)=-A(K,J)*S(K,J-1)*PN(K)
      F(K,J)=(F(K,J)-A(K,J)*F(K,J-1))*PN(K)
   20 CONTINUE

!C     BACKWARD PASS
      DO 30 K=L1,L2
      S(K,J2)=1.
      QE(K,J2)=0.
   30 CONTINUE
      DO 40 I=JA,J2
      J=JJ-I
      DO 40 K=L1,L2
      S(K,J)=S(K,J)+Q(K,J)*S(K,J+1)
      QE(K,J)=F(K,J)+Q(K,J)*QE(K,J+1)
   40 CONTINUE
      DO 50 K=L1,L2
      F(K,J2)=(FN(K)-C(K,J2)*QE(K,J1)-A(K,J2)*QE(K,J2-1))
     %       /(C(K,J2)*S(K,J1)+A(K,J2)*S(K,J2-1)+B(K,J2))
   50 CONTINUE

!C     BACKWARD ELIMINATION PASS
      DO 60 I=JA,J2
      J=JJ-I
      DO 60 K=L1,L2
      F(K,J)=F(K,J2)*S(K,J)+QE(K,J)
   60 CONTINUE

      RETURN
      END

!C*******************************************************************
      SUBROUTINE TRDIAG2P(A,B,C,F,J1,J2,L1,L2)
!C*******************************************************************
!C     INVERT A PERIODIC MATRIX (X-2 TRIDIAGONAL) WITH
!C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      REAL*8 A(M1,M2),B(M1,M2),C(M1,M2),F(M1,M2),Q(M1,M2),
     %     S(M1,M2),QE(M1,M2),FN(M1),PN(M1)

      JA=J1+1
      JJ=J1+J2
      DO 10 K=L1,L2
      BINV=1./B(K,J1)
      Q(K,J1)=-C(K,J1)*BINV
      S(K,J1)=-A(K,J1)*BINV
      FN(K)=F(K,J2)
      F(K,J1)=F(K,J1)*BINV
   10 CONTINUE

C     FORWARD ELIMINATION SWEEP
      DO 20 J=JA,J2
      DO 20 K=L1,L2
      PN(K)=1./(B(K,J)+A(K,J)*Q(K,J-1))
      Q(K,J)=-C(K,J)*PN(K)
      S(K,J)=-A(K,J)*S(K,J-1)*PN(K)
      F(K,J)=(F(K,J)-A(K,J)*F(K,J-1))*PN(K)
   20 CONTINUE

C     BACKWARD PASS
      DO 30 K=L1,L2
      S(K,J2)=1.
      QE(K,J2)=0.
   30 CONTINUE
      DO 40 I=JA,J2
      J=JJ-I
      DO 40 K=L1,L2
      S(K,J)=S(K,J)+Q(K,J)*S(K,J+1)
      QE(K,J)=F(K,J)+Q(K,J)*QE(K,J+1)
   40 CONTINUE
      DO 50 K=L1,L2
      F(K,J2)=(FN(K)-C(K,J2)*QE(K,J1)-A(K,J2)*QE(K,J2-1))
     %       /(C(K,J2)*S(K,J1)+A(K,J2)*S(K,J2-1)+B(K,J2))
   50 CONTINUE

C     BACKWARD ELIMINATION PASS
      DO 60 I=JA,J2
      J=JJ-I
      DO 60 K=L1,L2
      F(K,J)=F(K,J2)*S(K,J)+QE(K,J)
   60 CONTINUE

      RETURN
      END

!C*******************************************************************
      SUBROUTINE TRDIAG3P(A,B,C,F,J1,J2,L1,L2)
!C*******************************************************************
!C     INVERT A PERIODIC MATRIX (X-3 TRIDIAGONAL) WITH
!C     DIFFERENT COEFFICIENTS.

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      REAL*8 A(M1,M3),B(M1,M3),C(M1,M3),F(M1,M3),Q(M1,M3),
     %     S(M1,M3),QE(M1,M3),FN(M1),PN(M1)

      JA=J1+1
      JJ=J1+J2
      DO 10 K=L1,L2
      BINV=1./B(K,J1)         
      Q(K,J1)=-C(K,J1)*BINV   
      S(K,J1)=-A(K,J1)*BINV
      FN(K)=F(K,J2)
      F(K,J1)=F(K,J1)*BINV    
   10 CONTINUE

C     FORWARD ELIMINATION SWEEP
      DO 20 J=JA,J2
      DO 20 K=L1,L2
      PN(K)=1./(B(K,J)+A(K,J)*Q(K,J-1))
      Q(K,J)=-C(K,J)*PN(K)             
      S(K,J)=-A(K,J)*S(K,J-1)*PN(K)
      F(K,J)=(F(K,J)-A(K,J)*F(K,J-1))*PN(K)
   20 CONTINUE

C     BACKWARD PASS
      DO 30 K=L1,L2
      S(K,J2)=1.
      QE(K,J2)=0.
   30 CONTINUE
      DO 40 I=JA,J2
      J=JJ-I
      DO 40 K=L1,L2
      S(K,J)=S(K,J)+Q(K,J)*S(K,J+1)
      QE(K,J)=F(K,J)+Q(K,J)*QE(K,J+1)
   40 CONTINUE
      DO 50 K=L1,L2
      F(K,J2)=(FN(K)-C(K,J2)*QE(K,J1)-A(K,J2)*QE(K,J2-1))
     %       /(C(K,J2)*S(K,J1)+A(K,J2)*S(K,J2-1)+B(K,J2))
   50 CONTINUE

C     BACKWARD ELIMINATION PASS
      DO 60 I=JA,J2
      J=JJ-I
      DO 60 K=L1,L2
      F(K,J)=F(K,J2)*S(K,J)+QE(K,J)
   60 CONTINUE

      RETURN
      END

!C----------------------------- PTDIAG1A -------------------------------
!C     SOLVE THE PENTADIAGONAL MATRIX (X-1 PENTADIAGONAL)

      SUBROUTINE PTDIAG1A(A,B,C,D,E,F,UU,I1,IN,J1,JN,O,Q,R)

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      DIMENSION A(M2,M1),B(M2,M1),C(M2,M1),D(M2,M1),E(M2,M1)
      DIMENSION O(M2,M1),Q(M2,M1),R(M2,M1),F(M2,M1),UU(M2,M1)
      DIMENSION PDENO2(M2),PDENOJ(M2)

      J2=J1+1
      J3=J2+1
      JNM=JN-1
      DO 1 I=I1,IN
      O(I,J1)=-B(I,J1)/A(I,J1)
      Q(I,J1)=-C(I,J1)/A(I,J1)
      R(I,J1)=F(I,J1)/A(I,J1)
    1 CONTINUE
      DO 2 I=I1,IN
      PDENO2(I)=1./(A(I,J2)+D(I,J2)*O(I,J1))
      O(I,J2)=-(B(I,J2)+D(I,J2)*Q(I,J1))*PDENO2(I)
      Q(I,J2)=-C(I,J2)*PDENO2(I)
      R(I,J2)=(F(I,J2)-D(I,J2)*R(I,J1))*PDENO2(I)
    2 CONTINUE
      DO 10 J=J3,JN
      JM=J-1
      JMM=JM-1
      DO 10 I=I1,IN
      PDENOJ(I)=1./(A(I,J)+E(I,J)*Q(I,JMM)
     1             +(D(I,J)+E(I,J)*O(I,JMM))*O(I,JM))
      O(I,J)=-(B(I,J)+(D(I,J)+E(I,J)*O(I,JMM))*Q(I,JM))*PDENOJ(I)
      Q(I,J)=-C(I,J)*PDENOJ(I)
      R(I,J)=(F(I,J)-E(I,J)*R(I,JMM)-(D(I,J)+E(I,J)*O(I,JMM))*R(I,JM))
     1      *PDENOJ(I)
   10 CONTINUE
      DO 11 I=I1,IN
      UU(I,JN)=R(I,JN)
   11 CONTINUE
      DO 12 I=I1,IN
      UU(I,JNM)=O(I,JNM)*UU(I,JN)+R(I,JNM)
   12 CONTINUE
      DO 20 J=JNM-1,J1,-1
      DO 20 I=I1,IN
      UU(I,J)=O(I,J)*UU(I,J+1)+Q(I,J)*UU(I,J+2)+R(I,J)
   20 CONTINUE

      RETURN
      END


!C******************************************************************
      SUBROUTINE CAL_RHSPS(U,V,W,ANPSO,PSI_CN,DF_C)
!C******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE IBM_VAR
      USE HEAT_VAR
      
      USE TWO_PHASE_PROPERTY
      
      REAL*8 PSI_CN(0:M1,0:M2,0:M3)
      REAL*8 DF_C(0:M1,0:M2,0:M3,3)

      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL*8 ANPSO(M1M,M2M,M3M)

!!!!! FOR MEAN TEMPERATURE!!!!!!!
      TGR=0D0
      if (ibmon .eq. 1) then 
       TGR_TMP=0.
!$OMP PARALLEL DO private(I,J,KP1,DENSC,UU)
!$OMP&reduction(+:TGR_TMP)
       DO K=1,N3M
        KP1=KPV(K)
       DO J=1,N2M
       DO I=1,N1M
         DENSC=DENSCM+DENSC_DIFF*PSI_CN(I,J,K)
         UU=0.5*(W(I,J,KP1)+W(I,J,K))

         TGR_TMP=TGR_TMP+DENSC*UU*(SDX(I)*SDY(J)*SDZ(K))
       ENDDO
       ENDDO
       ENDDO

       TGR=-HEAT_SUM/TGR_TMP !0.5*YL IS R (CHARACTERISTIC LENGTH)
      endif
       
       
      REPRI=1./(RE*PRM)
      
       IF (DEN_DIFF .EQ. 0.) THEN
        DEN_DIFFI=1.
       ELSE
        DEN_DIFFI=1./DEN_DIFF
       ENDIF
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1,PSI_TMP)
!$OMP&private(DENSCE,DENSCW,DENSCN,DENSCS,DENSCT,DENSCB)
!$OMP&private(TCE,TCW,TCN,TCS,TCT,TCB)
!$OMP&private(PSE,PSW,PSC,PSF,PSN,PSS,ANPSX,ANPSY,ANPSZ,ANPS)
!$OMP&private(ALPSE,ALPSW,ALPSN,ALPSS,ALPSC,ALPSF,ALPSX,ALPSY,ALPSZ)
!$OMP&private(ALPS,DENSC,UU,AMEANT1,AMEANT2)
      DO K=1,N3M
        KP1=KPV(K)
        KM1=KMV(K)
      DO J=1,N2M
        JP1=JPV(J)
        JM1=JMV(J)
      DO I=1,N1M
        IP1=IPV(I)
        IM1=IMV(I)

       PSI_TMP=(DF_C(I,J,K,1)-DENM)*DEN_DIFFI
       DENSCE=DENSCM+DENSC_DIFF*PSI_TMP
       TCE=TCM+TC_DIFF*PSI_TMP
         
       PSI_TMP=(DF_C(IM1,J,K,1)-DENM)*DEN_DIFFI
       DENSCW=DENSCM+DENSC_DIFF*PSI_TMP
       TCW=TCM+TC_DIFF*PSI_TMP
       
       PSI_TMP=(DF_C(I,J,K,2)-DENM)*DEN_DIFFI
       DENSCN=DENSCM+DENSC_DIFF*PSI_TMP
       TCN=TCM+TC_DIFF*PSI_TMP

       PSI_TMP=(DF_C(I,JM1,K,2)-DENM)*DEN_DIFFI
       DENSCS=DENSCM+DENSC_DIFF*PSI_TMP
       TCS=TCM+TC_DIFF*PSI_TMP
       
       PSI_TMP=(DF_C(I,J,K,3)-DENM)*DEN_DIFFI
       DENSCT=DENSCM+DENSC_DIFF*PSI_TMP
         TCT=TCM+TC_DIFF*PSI_TMP
       
       PSI_TMP=(DF_C(I,J,KM1,3)-DENM)*DEN_DIFFI
       DENSCB=DENSCM+DENSC_DIFF*PSI_TMP
         TCB=TCM+TC_DIFF*PSI_TMP

      !CONVECTION
      PSE=0.5*(SDX(IP1)*T(I,J,K)+SDX(I)*T(IP1,J,K))*VVDX(IP1) ! interpolation
      PSW=0.5*(SDX(IM1)*T(I,J,K)+SDX(I)*T(IM1,J,K))*VVDX(I)
      PSC=0.5*(SDZ(KP1)*T(I,J,K)+SDZ(K)*T(I,J,KP1))*VVDZ(KP1)
      PSF=0.5*(SDZ(KM1)*T(I,J,K)+SDZ(K)*T(I,J,KM1))*VVDZ(K)
      PSN=0.5*(SDY(JP1)*T(I,J,K)+SDY(J)*T(I,JP1,K))*VVDY(JP1)
     &   *(1.-FIXJU(J))+T(I,N2,K)*FIXJU(J)
      PSS=0.5*(SDY(JM1)*T(I,J,K)+SDY(J)*T(I,JM1,K))*VVDY(J)
     &   *(1.-FIXJL(J))+T(I,0,K)*FIXJL(J)
      ANPSX=(DENSCE*U(IP1,J,K)*PSE-DENSCW*U(I,J,K)*PSW)*SSDX(I)
      ANPSY=(DENSCN*V(I,JP1,K)*PSN-DENSCS*V(I,J,K)*PSS)*SSDY(J)
      ANPSZ=(DENSCT*W(I,J,KP1)*PSC-DENSCB*W(I,J,K)*PSF)*SSDZ(K)
      ANPS=0d0!-(ANPSX+ANPSY+ANPSZ)

      !DIFFUSION
      ALPSE=(T(IP1,J,K)-T(I,J,K))*VVDX(IP1)
      ALPSW=(T(I,J,K)-T(IM1,J,K))*VVDX(I)
      ALPSN=(T(I,JP1,K)-T(I,J,K))*VVDY(JP1)
      ALPSS=(T(I,J,K)-T(I,JM1,K))*VVDY(J)
      ALPSC=(T(I,J,KP1)-T(I,J,K))*VVDZ(KP1)
      ALPSF=(T(I,J,K)-T(I,J,KM1))*VVDZ(K)
      ALPSX=(TCE*ALPSE-TCW*ALPSW)*SSDX(I)
      ALPSY=(TCN*ALPSN-TCS*ALPSS)*SSDY(J)
      ALPSZ=(TCT*ALPSC-TCB*ALPSF)*SSDZ(K)
      ALPS=0d0!REPRI*(ALPSX+ALPSY+ALPSZ)


       DENSC=DENSCM+DENSC_DIFF*PSI_CN(I,J,K)
       UU=0.5*(W(I,J,KP1)+W(I,J,K))
       AMEANT1=DENSC*UU*TGR

       AMEANT2=REPRI*(TCt-TCb)*SSDz(k)*TGR

      RHSPS(I,J,K)=( ( GAMMA(MSUB)*ANPS+RO(MSUB)*ANPSO(I,J,K) )
     &              +2.*ALPHA_RK3*ALPS
     &              +2.*ALPHA_RK3*(AMEANT1-AMEANT2) )*DT


      ANPSO(I,J,K)=ANPS

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END


!C******************************************************************
      SUBROUTINE CAL_RHSPS2(U,V,W,ANPSO,PSI_CN,DF_C)
!C******************************************************************
      USE FLOW_VAR
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE IBM_VAR
      USE HEAT_VAR     
      USE TWO_PHASE_PROPERTY
      USE PHASE_BOILING
      
      REAL*8 PSI_CN(0:M1,0:M2,0:M3)
      REAL*8 DF_C(0:M1,0:M2,0:M3,3)
      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8 ANPSO(M1M,M2M,M3M)
      REAL*8 PSI_TMP6,PSI_TMP5,PSI_TMP4,PSI_TMP3,PSI_TMP2,PSI_TMP1
      INTEGER*8 IM
!!!!! FOR MEAN TEMPERATURE!!!!!!!
      TGR=0D0
      if (ibmon .eq. 1) then 
       TGR_TMP=0.
!$OMP PARALLEL DO private(I,J,KP1,DENSC,UU)
!$OMP&reduction(+:TGR_TMP)
       DO K=1,N3M
        KP1=KPV(K)
       DO J=1,N2M
       DO I=1,N1M
         DENSC=DENSCM+DENSC_DIFF*PSI_CN(I,J,K)
         UU=0.5*(W(I,J,KP1)+W(I,J,K))
         TGR_TMP=TGR_TMP+DENSC*UU*(SDX(I)*SDY(J)*SDZ(K))
       ENDDO
       ENDDO
       ENDDO
       TGR=-HEAT_SUM/TGR_TMP !0.5*YL IS R (CHARACTERISTIC LENGTH)
      endif
       
       
      REPRI=1./(RE*PRM)
      
       IF (DEN_DIFF .EQ. 0.) THEN
        DEN_DIFFI=1.
       ELSE
        DEN_DIFFI=1./DEN_DIFF
       ENDIF
!$OMP PARALLEL DO private(I,J,IP1,IM1,JP1,JM1,KP1,KM1,PSI_TMP)
!$OMP&private(DENSCE,DENSCW,DENSCN,DENSCS,DENSCT,DENSCB)
!$OMP&private(TCE,TCW,TCN,TCS,TCT,TCB)
!$OMP&private(PSE,PSW,PSC,PSF,PSN,PSS,ANPSX,ANPSY,ANPSZ,ANPS)
!$OMP&private(ALPSE,ALPSW,ALPSN,ALPSS,ALPSC,ALPSF,ALPSX,ALPSY,ALPSZ)
!$OMP&private(ALPS,DENSC,UU,AMEANT1,AMEANT2)
!$OMP&private(IM,PSI_TMP6,PSI_TMP5,PSI_TMP4,PSI_TMP3,PSI_TMP2,PSI_TMP1)
      DO K=1,N3M
        KP1=KPV(K)
        KM1=KMV(K)
      DO J=1,N2M
        JP1=JPV(J)
        JM1=JMV(J)
      DO I=1,N1M
        IP1=IPV(I)
        IM1=IMV(I)
        
       PSI_TMP=PSI_CN(I,J,K)
!------1CELL AROUND THE INTERFACE CELL-----------------------------------      
      PSI_TMP6=(DF_C(I,J,KM1,3)-DENM)*DEN_DIFFI
      PSI_TMP5=(DF_C(I,J,K  ,3)-DENM)*DEN_DIFFI 
      PSI_TMP4=(DF_C(I,JM1,K,2)-DENM)*DEN_DIFFI 
      PSI_TMP3=(DF_C(I,J  ,K,2)-DENM)*DEN_DIFFI 
      PSI_TMP2=(DF_C(IM1,J,K,1)-DENM)*DEN_DIFFI 
      PSI_TMP1=(DF_C(I  ,J,K,1)-DENM)*DEN_DIFFI 
      IF ((DABS(PSI_TMP -0.5D0) .NE. 0.5D0) .OR.  ! INTERFACE 자체거나     
     &    (DABS(PSI_TMP1-0.5D0) .NE. 0.5D0) .OR.  ! 차분에 사용하는 CELL에 INTERFACE가 있거나
     &    (DABS(PSI_TMP2-0.5D0) .NE. 0.5D0) .OR.
     &    (DABS(PSI_TMP3-0.5D0) .NE. 0.5D0) .OR.
     &    (DABS(PSI_TMP4-0.5D0) .NE. 0.5D0) .OR.
     &    (DABS(PSI_TMP5-0.5D0) .NE. 0.5D0) .OR.
     &    (DABS(PSI_TMP6-0.5D0) .NE. 0.5D0)) THEN
      
       IF (PSI_TMP .GT. 0.5D0) THEN ! LIQUID
       TCE=TCP;       TCW=TCP;
       TCN=TCP;       TCS=TCP;
       TCT=TCP;       TCB=TCP
       DENSCE=DENSCP;DENSCW=DENSCP;
       DENSCN=DENSCP;DENSCS=DENSCP;
       DENSCT=DENSCP;DENSCB=DENSCP       
       IM=1
       ELSEIF (PSI_TMP .LE. 0.5D0) THEN ! GAS
       TCE=TCM;       TCW=TCM;
       TCN=TCM;       TCS=TCM;
       TCT=TCM;       TCB=TCM      
       DENSCE=DENSCM;DENSCW=DENSCM;
       DENSCN=DENSCM;DENSCS=DENSCM;
       DENSCT=DENSCM;DENSCB=DENSCM 
       IM=0
       ENDIF
      !CONVECTION
      PSE=0.5D0*(SDX(IP1)*TEMP(I,J,K,IM)+SDX(I)*TEMP(IP1,J,K,IM))
     &                                                        *VVDX(IP1) ! interpolation
      PSW=0.5D0*(SDX(IM1)*TEMP(I,J,K,IM)+SDX(I)*TEMP(IM1,J,K,IM))
     &                                                        *VVDX(I)
      PSC=0.5D0*(SDZ(KP1)*TEMP(I,J,K,IM)+SDZ(K)*TEMP(I,J,KP1,IM))
     &                                                        *VVDZ(KP1)
      PSF=0.5D0*(SDZ(KM1)*TEMP(I,J,K,IM)+SDZ(K)*TEMP(I,J,KM1,IM))
     &                                                        *VVDZ(K)
      PSN=0.5D0*(SDY(JP1)*TEMP(I,J,K,IM)+SDY(J)*TEMP(I,JP1,K,IM))
     &    *VVDY(JP1) *(1.-FIXJU(J))+TEMP(I,N2,K,IM)*FIXJU(J)
      PSS=0.5D0*(SDY(JM1)*TEMP(I,J,K,IM)+SDY(J)*TEMP(I,JM1,K,IM))
     &    *VVDY(J)   *(1.-FIXJL(J))+TEMP(I,0 ,K,IM)*FIXJL(J)
      ANPSX=(DENSCE*U(IP1,J,K)*PSE-DENSCW*U(I,J,K)*PSW)*SSDX(I)
      ANPSY=(DENSCN*V(I,JP1,K)*PSN-DENSCS*V(I,J,K)*PSS)*SSDY(J)
      ANPSZ=(DENSCT*W(I,J,KP1)*PSC-DENSCB*W(I,J,K)*PSF)*SSDZ(K)
      ANPS=-(ANPSX+ANPSY+ANPSZ)   
      !DIFFUSION
      ALPSE=(TEMP(IP1,J,K,IM)-TEMP(I  ,J,K,IM))*VVDX(IP1)
      ALPSW=(TEMP(I  ,J,K,IM)-TEMP(IM1,J,K,IM))*VVDX(I)
      ALPSN=(TEMP(I,JP1,K,IM)-TEMP(I,J  ,K,IM))*VVDY(JP1)
      ALPSS=(TEMP(I,J  ,K,IM)-TEMP(I,JM1,K,IM))*VVDY(J)
      ALPSC=(TEMP(I,J,KP1,IM)-TEMP(I,J,K  ,IM))*VVDZ(KP1)
      ALPSF=(TEMP(I,J,K  ,IM)-TEMP(I,J,KM1,IM))*VVDZ(K)
      ALPSX=(TCE*ALPSE-TCW*ALPSW)*SSDX(I)
      ALPSY=(TCN*ALPSN-TCS*ALPSS)*SSDY(J)
      ALPSZ=(TCT*ALPSC-TCB*ALPSF)*SSDZ(K)
      ALPS=REPRI*(ALPSX+ALPSY+ALPSZ)      
!-----NORMAL CELL-------------------------------------------------------       
      ELSE
       IF (PSI_TMP .EQ. 1D0) THEN ! LIQUID
        TCE=TCP;       TCW=TCP;
        TCN=TCP;       TCS=TCP;
        TCT=TCP;       TCB=TCP  
        DENSCE=DENSCP;DENSCW=DENSCP;
        DENSCN=DENSCP;DENSCS=DENSCP;
        DENSCT=DENSCP;DENSCB=DENSCP         
       ELSEIF (PSI_TMP .EQ. 0D0) THEN ! GAS
        TCE=TCM;       TCW=TCM;
        TCN=TCM;       TCS=TCM;
        TCT=TCM;       TCB=TCM
        DENSCE=DENSCM;DENSCW=DENSCM;
        DENSCN=DENSCM;DENSCS=DENSCM;
        DENSCT=DENSCM;DENSCB=DENSCM 
       ENDIF
      !CONVECTION
      PSE=0.5*(SDX(IP1)*T(I,J,K)+SDX(I)*T(IP1,J,K))*VVDX(IP1) ! interpolation
      PSW=0.5*(SDX(IM1)*T(I,J,K)+SDX(I)*T(IM1,J,K))*VVDX(I)
      PSC=0.5*(SDZ(KP1)*T(I,J,K)+SDZ(K)*T(I,J,KP1))*VVDZ(KP1)
      PSF=0.5*(SDZ(KM1)*T(I,J,K)+SDZ(K)*T(I,J,KM1))*VVDZ(K)
      PSN=0.5*(SDY(JP1)*T(I,J,K)+SDY(J)*T(I,JP1,K))*VVDY(JP1)
     &   *(1.-FIXJU(J))+T(I,N2,K)*FIXJU(J)
      PSS=0.5*(SDY(JM1)*T(I,J,K)+SDY(J)*T(I,JM1,K))*VVDY(J)
     &   *(1.-FIXJL(J))+T(I,0,K)*FIXJL(J)
      ANPSX=(DENSCE*U(IP1,J,K)*PSE-DENSCW*U(I,J,K)*PSW)*SSDX(I)
      ANPSY=(DENSCN*V(I,JP1,K)*PSN-DENSCS*V(I,J,K)*PSS)*SSDY(J)
      ANPSZ=(DENSCT*W(I,J,KP1)*PSC-DENSCB*W(I,J,K)*PSF)*SSDZ(K)
      ANPS=-(ANPSX+ANPSY+ANPSZ)   
      !DIFFUSION
      ALPSE=(T(IP1,J,K)-T(I,J,K))*VVDX(IP1)
      ALPSW=(T(I,J,K)-T(IM1,J,K))*VVDX(I)
      ALPSN=(T(I,JP1,K)-T(I,J,K))*VVDY(JP1)
      ALPSS=(T(I,J,K)-T(I,JM1,K))*VVDY(J)
      ALPSC=(T(I,J,KP1)-T(I,J,K))*VVDZ(KP1)
      ALPSF=(T(I,J,K)-T(I,J,KM1))*VVDZ(K)
      ALPSX=(TCE*ALPSE-TCW*ALPSW)*SSDX(I)
      ALPSY=(TCN*ALPSN-TCS*ALPSS)*SSDY(J)
      ALPSZ=(TCT*ALPSC-TCB*ALPSF)*SSDZ(K)
      ALPS=REPRI*(ALPSX+ALPSY+ALPSZ)      
      ENDIF
!!DIFFUSION
!       ALPSE=(T(IP1,J,K)-T(I,J,K))*VVDX(IP1)
!       ALPSW=(T(I,J,K)-T(IM1,J,K))*VVDX(I)
!       ALPSN=(T(I,JP1,K)-T(I,J,K))*VVDY(JP1)
!       ALPSS=(T(I,J,K)-T(I,JM1,K))*VVDY(J)
!       ALPSC=(T(I,J,KP1)-T(I,J,K))*VVDZ(KP1)
!       ALPSF=(T(I,J,K)-T(I,J,KM1))*VVDZ(K)
!       ALPSX=(TCE*ALPSE-TCW*ALPSW)*SSDX(I)
!       ALPSY=(TCN*ALPSN-TCS*ALPSS)*SSDY(J)
!       ALPSZ=(TCT*ALPSC-TCB*ALPSF)*SSDZ(K)
!       ALPS=REPRI*(ALPSX+ALPSY+ALPSZ)        



       DENSC=DENSCM+DENSC_DIFF*PSI_CN(I,J,K)
       UU=0.5*(W(I,J,KP1)+W(I,J,K))
       AMEANT1=DENSC*UU*TGR

       AMEANT2=REPRI*(TCt-TCb)*SSDz(k)*TGR

      RHSPS(I,J,K)=( ( GAMMA(MSUB)*ANPS+RO(MSUB)*ANPSO(I,J,K) )
     &              +2.*ALPHA_RK3*ALPS
     &              +2.*ALPHA_RK3*(AMEANT1-AMEANT2) )*DT
      
      if (iphs.eq.1) then 
      
      RHSPS(I,J,K)=( GAMMA(MSUB)*T(I,J,K)*
     & (W(I,J,KP1)-W(I,J,K))*VVDZ(K) )*DT
!     & ((UE-UW)*SSDX(I)+(VN-VS)*SSDY(J)+(WT-WB)*VVDZ(K))  )*DT
      endif
      
      ANPSO(I,J,K)=ANPS

      ENDDO
      ENDDO
      ENDDO

      RETURN
      END
      
!------------------------------- LHSPS --------------------------------
!******************************************************************
      SUBROUTINE LHSPS(PSI_CN,PSI_C,DF_C)
!******************************************************************
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE HEAT_VAR
      
      USE TWO_PHASE_PROPERTY

      ! REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

      REAL*8 PSI_CN(0:M1,0:M2,0:M3)
      REAL*8 DF_C(0:M1,0:M2,0:M3,3)
      REAL*8 PSI_C(0:M1,0:M2,0:M3)

      real, dimension (:,:), allocatable :: AI,BI,CI,GI
      real, dimension (:,:), allocatable :: AJ,BJ,CJ,GJ
      real, dimension (:,:), allocatable :: AK,BK,CK,GK

       IF (DEN_DIFF .EQ. 0.) THEN
        DEN_DIFFI=1.
       ELSE
        DEN_DIFFI=1./DEN_DIFF
       ENDIF

      CREPR=0.!FLOAT(ILES)*RE*PRA
      ACOPS=ALPHA_RK3*DT/(RE*PRM)!0.5*DT*REI*PRAI

!$OMP PARALLEL private(I,K,DK1,DK2,DK3)
!$OMP&private(AK,CK,BK,GK,DENSCF,COEF,PSI_TMP,TCT,TCB,TALPHC,TALPHF)
      allocate(AK(M1,M3),BK(M1,M3),CK(M1,M3),GK(M1,M3))
!$OMP DO
      DO J=1,N2M

      DO K=1,N3M
      DO I=1,N1M
        DENSCF=DENSCM+DENSC_DIFF*PSI_C(I,J,K)
        COEF=ACOPS/DENSCF
        PSI_TMP=(DF_C(I,J,KMV(K),3)-DENM)*DEN_DIFFI
        TCB=TCM+TC_DIFF*PSI_TMP 
        PSI_TMP=(DF_C(I,J,K,3)-DENM)*DEN_DIFFI
        TCT=TCM+TC_DIFF*PSI_TMP

        !if (iphs .eq. 1) then 
        !if (psi_c(i,j,k) .lt.0.5d0) then 
        !coef=acops/denscm
        !TCT=TCM
        !TCB=TCM            
        !else
        !coef=acops/denscp
        !TCT=TCP
        !TCB=TCP            
        !endif
        !endif

!!!!$OMP PARALLEL DO private(KP)
!!!      DO 60 KC=1,N3M
!!!      KP=KPV(KC)
!!!      AKUV(KC)=-VVDZ(KC)*SSDZ(KC)
!!!      CKUV(KC)=-VVDZ(KP)*SSDZ(KC)
!!!   60 CONTINUE
      AK(I,K)=COEF*TCB*AKUV(K)!*(1.-FIXKL(K))
      CK(I,K)=COEF*TCT*CKUV(K)
      BK(I,K)=1.-(AK(I,K)+CK(I,K))
      GK(I,K)=RHSPS(I,J,K)
      
      if (dabs(psi_cn(i,j,k)-0.5d0) .ne. 0.5) then 
      if (psi_cn(i,j,k) .le. 0.5d0) then 
      ck(i,k)=0d0
      elseif (psi_cn(i,j,k) .gt. 0.5d0) then 
      ak(i,k)=0d0
      endif
      endif
      
      ENDDO
      ENDDO

      IF (IPZ .EQ. 1) THEN
       CALL TRDIAG3P(AK,BK,CK,GK,1,N3M,1,N1M)
      ELSE
       CALL TRDIAG3(AK,BK,CK,GK,GK,1,N3M,1,N1M)
      ENDIF

      DO K=1,N3M
      DO I=1,N1M
      RHSPS(I,J,K)=GK(I,K)
      ENDDO
      ENDDO
      
!     UPDATE T
      DO K=1,N3M
      DO I=1,N1M
      T(I,J,K)=GK(i,k)+T(I,J,K)
      ENDDO
      ENDDO         
      ENDDO

!$OMP END DO
      deallocate(AK,BK,CK,GK)
!$OMP END PARALLEL

 1000 CONTINUE

!     PERIODICITY

!$OMP PARALLEL DO private(I)
      DO J=1,N2M
      DO I=1,N1M
         T(I,J,0) =0d0
         T(I,J,N3M)=1d0 
         T(I,J,N3)=T(I,J,N3M)        
      ENDDO
      ENDDO  

      RETURN
      END


!*****************************************************************
!
!! REAL_TIME returns the REAL*8 time in seconds.

      subroutine real_time(seconds)

      implicit none
      INTEGER*8 clock_count,clock_max,clock_rate
      REAL*8 seconds

      call system_clock ( clock_count, clock_rate, clock_max )

      seconds=real(clock_count)/real(clock_rate)

      return
      end


!*******************************************************************************
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Modified:
!
!    06 August 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
      subroutine timestamp
        implicit none

        character ( len = 8 ) ampm
        INTEGER*8 d
        INTEGER*8 h
        INTEGER*8 m
        INTEGER*8 mm
        character ( len = 9 ), parameter, dimension(12) :: month = (/
     &    'January  ', 'February ', 'March    ', 'April    ',
     &    'May      ', 'June     ', 'July     ', 'August   ',
     &    'September', 'October  ', 'November ', 'December ' /)
        INTEGER*8 n
        INTEGER*8 s
        INTEGER*8 values(8)
        INTEGER*8 y

        call date_and_time ( values = values )

        y = values(1)
        m = values(2)
        d = values(3)
        h = values(5)
        n = values(6)
        s = values(7)
        mm = values(8)

        if ( h < 12 ) then
          ampm = 'AM'
        else if ( h == 12 ) then
          if ( n == 0 .and. s == 0 ) then
            ampm = 'Noon'
          else
            ampm = 'PM'
          end if
        else
          h = h - 12
          if ( h < 12 ) then
            ampm = 'PM'
          else if ( h == 12 ) then
            if ( n == 0 .and. s == 0 ) then
              ampm = 'Midnight'
            else
              ampm = 'AM'
            end if
          end if
        end if

      write ( *, '(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' )
     &d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )
      write(*,*)' '

        return
      end



!******************************************************************
!       SECTION FOR GEOMETRY
!      - READ GRID FILE & COMPUTE VARIABLES ABOUT GEOMETRY
!      - GEOMETRY VARIABLES ARE LISTED IN 'GEOM.H'
!      - DATA READ FROM FILE ARE BLOCK 'DIM','COORD','SCALES',
!        'GEOMINPUT'
!      - COMPUTED VARIABLES ARE BLOCK 'GEOMETRY','INDEX','FIX',
!        'POSITION','VAR','VVAR'
!
!                               Revised by DONGJOO KIM
!                               Turbulence & Flow Control Lab.
!                               School of Mech. & Aerospace Eng.
!                               Seoul National Univ.
!                               9/27/1999
!
!
!      -REVISED FOR TWO-PHASE FLOW
!
!                               09/03/2015
!******************************************************************

      SUBROUTINE GEOM(gridfile)

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING

      CHARACTER*8 gridfile

      !FLOW_GEOM_VAR
      ALLOCATE(IPV(0:M1+1),JPV(0:M2+1),KPV(0:M3+1)
     &                ,IMV(-1:M1),JMV(-1:M2),KMV(-1:M3))
      ALLOCATE(FIXIL(M1M),FIXIU(M1M),FIXJL(M2M),FIXJU(M2M)
     &            ,FIXKL(M3M),FIXKU(M3M))
      ALLOCATE( X(0:M1+1),Y(0:M2+1),Z(0:M3+1))
      ALLOCATE( XP(0:M1+1),YP(0:M2+1),ZP(0:M3+1))
      ALLOCATE( SDX(0:M1),SDY(0:M2),SDZ(0:M3),VDX(M1),VDY(M2),VDZ(M3))
      ALLOCATE( SSDX(0:M1),SSDY(0:M2),SSDZ(0:M3)
     &             ,VVDX(M1),VVDY(M2),VVDZ(M3))

      !LVS_GEOM
      ALLOCATE( XF(M1L:M1U),YF(M2L:M2U),ZF(M3L:M3U))
      ALLOCATE( XPF(M1L:M1U),YPF(M2L:M2U),ZPF(M3L:M3U))
      ALLOCATE( IPF(0:M1F+1),IMF(-1:M1F),JPF(0:M2F+1),JMF(-1:M2F)
     &,KPF(0:M3F+1),KMF(-1:M3F))

      !LVS_COUPLING
      ALLOCATE( ICOU_VEL(M1F),JCOU_VEL(M2F),KCOU_VEL(M3F))
      ALLOCATE( ICOUMP_VEL(M1F),JCOUMP_VEL(M2F),KCOUMP_VEL(M3F))
      ALLOCATE( ICOU1(-1:M1+2),JCOU1(-1:M2+2),KCOU1(-1:M3+2),
     &          ICOUMP1(-1:M1+2),JCOUMP1(-1:M2+2),KCOUMP1(-1:M3+2))
      ALLOCATE( ICOU2(-1:M1+2),JCOU2(-1:M2+2),KCOU2(-1:M3+2),
     &          ICOUMP2(-1:M1+2),JCOUMP2(-1:M2+2),KCOUMP2(-1:M3+2))

      OPEN(11,FILE=gridfile)
      CALL READGRID
      CLOSE(11)

       CALL GEOMETRIC
       CALL GEOMETRIC_LVS
       CALL INDICES
       CALL INDICES_LVS
!      CALL GEOM_COUPLING   !LAST LINE IN THIS SUBROUTINE.

      IF (IPX .EQ. 1) THEN
      FIXIL(1)=0.
      FIXIU(N1M)=0.
      SDX(0)=SDX(N1M)
      SDX(N1)=SDX(1)
      SSDX(0)=SSDX(N1M)
      SSDX(N1)=SSDX(1)
      VDX(1)=0.5*(SDX(1)+SDX(N1M))
      VDX(N1)=0.5*(SDX(1)+SDX(N1M))
      X(0)=X(1)-SDX(1)
      XP(0)=XP(1)-VDX(1)
      XP(N1)=XP(N1M)+VDX(N1)
      DO 48 I=1,N1
   48 VVDX(I)=1./VDX(I)
      IMV(1)=N1M
      IPV(N1M)=1
      !FOR LEVEL-SET
        IPF(N1FM)=1
        IMF(1)=N1FM
      ENDIF
      IF (IPX.EQ.1) THEN
         IBG=1
      ELSE
         IBG=2
      ENDIF
      
      IF (IPY .EQ. 1) THEN
      FIXJL(1)=0.
      FIXJU(N2M)=0.
      SDY(0)=SDY(N2M)
      SDY(N2)=SDY(1)
      SSDY(0)=SSDY(N2M)
      SSDY(N2)=SSDY(1)
      VDY(1)=0.5*(SDY(1)+SDY(N2M))
      VDY(N2)=0.5*(SDY(1)+SDY(N2M))
      Y(0)=Y(1)-SDY(1)
      YP(0)=YP(1)-VDY(1)
      YP(N2)=YP(N2M)+VDY(N2)
      DO 49 J=1,N2
   49 VVDY(J)=1./VDY(J)
      JMV(1)=N2M
      JPV(N2M)=1
      !FOR LEVEL-SET
        JPF(N2FM)=1
        JMF(1)=N2FM
      ENDIF
      IF (IPY.EQ.1) THEN
         JBG=1
      ELSE
         JBG=2
      ENDIF
      
      IF (IPZ .EQ. 1) THEN
      FIXKL(1)=0.
      FIXKU(N3M)=0.
      SDZ(0)=SDZ(N3M)
      SDZ(N3)=SDZ(1)
      SSDZ(0)=SSDZ(N3M)
      SSDZ(N3)=SSDZ(1)
      VDZ(1)=0.5*(SDZ(1)+SDZ(N3M))
      VDZ(N3)=0.5*(SDZ(1)+SDZ(N3M))
      Z(0)=Z(1)-SDZ(1)
      ZP(0)=ZP(1)-VDZ(1)
      ZP(N3)=ZP(N3M)+VDZ(N3)
      DO 47 K=1,N3
   47 VVDZ(K)=1./VDZ(K)
      KMV(1)=N3M
      KPV(N3M)=1
      !FOR LEVEL-SET
        KPF(N3FM)=1       
        KMF(1)=N3FM
      ENDIF
      IF (IPZ.EQ.1) THEN
         KBG=1
      ELSE
         KBG=2
      ENDIF

      CALL GEOM_COUPLING !need ipf,imf,kpf,kmf

      RETURN
      END

!******************************************************************
       SUBROUTINE READGRID
!******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_GEOM_VAR

      !FOR NAVIER-STOKES
      READ(11,*) N1,N2,N3
      READ(11,*) XL,YL,ZL

      READ(11,*) (X(I),I=1,N1)
      READ(11,*) (Y(J),J=1,N2)
      READ(11,*) (Z(K),K=1,N3)

      IF ((N1.GT.M1) .OR. (N2.GT.M2) .OR. (N3.GT.M3)) THEN
      write(*,*) n1,m1,n2,m2,n3,m3
         PRINT*, 'NSE-ARRAY SIZE CAN NOT HANDLE THIS GRID.'
         STOP
      END IF
      
      !FOR LEVEL-SET(G-GRID)
      READ(11,*) N1F,N2F,N3F
      READ(11,*) XLG,YLG,ZLG

      READ(11,*) (XF(I),I=1,N1F)
      READ(11,*) (YF(J),J=1,N2F)
      READ(11,*) (ZF(K),K=1,N3F)

      IF ((N1F.GT.M1F) .OR. (N2F.GT.M2F) .OR. (N3F.GT.M3F)) THEN
         PRINT*, 'LVS-ARRAY SIZE CAN NOT HANDLE THIS GRID.'
         STOP
      END IF

      N1M=N1-1
      N2M=N2-1
      N3M=N3-1
     
      N1FM=N1F-1
      N2FM=N2F-1
      N3FM=N3F-1

      N1F_BD=DFLOAT(N1FM)/DFLOAT(N1M)
      N2F_BD=DFLOAT(N2FM)/DFLOAT(N2M)
      N3F_BD=DFLOAT(N3FM)/DFLOAT(N3M)

      X(0)=X(1) !X(0),Y(0),Z(0) IS ADDED FOR GRID_COUPLING
      Y(0)=Y(1) !WHEN USING PERIODIC CONDITION
      Z(0)=Z(1)      

      !GRID_INFORMATION_WRITE
      WRITE(*,101) N1,N2,N3,N1M*N2*N3M       
      WRITE(*,102) N1FM,N2FM,N3FM,N1FM*N2FM*N3FM
      WRITE(*,103) XL,YL,ZL      
 101  FORMAT('N1 = ',I4,' N2 = ',I4,' N3 = ',I4,' NSE_TOT = ',I10)
 102  FORMAT('N1F =',I4,' N2F =',I4,' N3F =',I4,' REFINE_TOT = ',I10) 
 103  FORMAT('XL = ',F10.5,' YL = ',F10.5,' ZL = ',F10.5) 

      RETURN
      END
      
!******************************************************************
       SUBROUTINE GEOMETRIC
!******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      DO 10 I=1,N1M
   10 SDX(I)=X(I+1)-X(I)
      DO 20 I=2,N1M
   20 VDX(I)=0.5*(SDX(I)+SDX(I-1))
      VDX(1)=0.5*SDX(1)             ! LEFT BOUNDED COND.
      VDX(N1)=0.5*SDX(N1M)          ! RIGHT BOUNDED COND.

      DO 30 I=1,N1M
   30 SSDX(I)=1./SDX(I)
      DO 40 I=1,N1
   40 VVDX(I)=1./VDX(I)

      DO 15 J=1,N2M
   15 SDY(J)=Y(J+1)-Y(J)
      DO 25 J=2,N2M
   25 VDY(J)=0.5*(SDY(J)+SDY(J-1))
      VDY(1)=0.5*SDY(1)             ! LOWER WALL BOUNDED COND.
      VDY(N2)=0.5*SDY(N2M)          ! UPPER WALL BOUNDED COND.

      DO 35 J=1,N2M
   35 SSDY(J)=1./SDY(J)
      DO 45 J=1,N2
   45 VVDY(J)=1./VDY(J)

      IF (N3M.EQ.1) THEN
      ZL=1.
      SDZ(1)=ZL/N3M
      SSDZ(1)=1./SDZ(1)
      VDZ(1)=SDZ(1)
      VVDZ(1)=SSDZ(1)

      ELSE
      DO 17 K=1,N3M
   17 SDZ(K)=Z(K+1)-Z(K)
      DO 27 K=2,N3M
   27 VDZ(K)=0.5*(SDZ(K)+SDZ(K-1))
      VDZ(1)=0.5*SDZ(1)             ! LOWER WALL BOUNDED COND.
      VDZ(N3)=0.5*SDZ(N3M)          ! UPPER WALL BOUNDED COND.

      DO 37 K=1,N3M
   37 SSDZ(K)=1./SDZ(K)
      DO 47 K=1,N3
   47 VVDZ(K)=1./VDZ(K)
      ENDIF

C-----set by zero
      SDX(0)=0.
      SDX(M1)=0.
      SDY(0)=0.
      SDY(M2)=0.
      SDZ(0)=0.
      SDZ(M3)=0.

      DO 5 I=1,N1M
    5 XP(I)=X(I)+0.5*SDX(I)
      XP(N1)=X(N1)
      XP(0)=X(1)

      DO 6 J=1,N2M
    6 YP(J)=Y(J)+0.5*SDY(J)
      YP(N2)=Y(N2)
      YP(0)=Y(1)

      DO 7 K=1,N3M
    7 ZP(K)=Z(K)+0.5*SDZ(K)
      ZP(N3)=Z(N3)
      ZP(0)=Z(1)

      RETURN
      END
      
!******************************************************************
       SUBROUTINE GEOMETRIC_LVS
!******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_GEOM_VAR



        SDXF=XF(2)-XF(1)
        SDYF=YF(2)-YF(1)   
        SDZF=ZF(2)-ZF(1)

        SSDXF=1./SDXF
        SSDYF=1./SDYF
        SSDZF=1./SDZF
        
        !X-DIR
        DO I=1,N1FM
        	XPF(I)=XF(I)+0.5*SDXF
        ENDDO
        DO I=0,2*N1F_BD
          XF(-I)=XF(1)-SDXF*(I+1)
          XF(N1F+1+I)=XF(N1F)+SDXF*(I+1)
        ENDDO
        DO I=0,2*N1F_BD
          XPF(-I)=XF(-I)+0.5*SDXF
          XPF(N1F+I)=XF(N1F+I)+0.5*SDXF
        ENDDO

        !Y-DIR
        DO J=1,N2FM
        	YPF(J)=YF(J)+0.5*SDYF
        ENDDO
        DO J=0,2*N2F_BD
          YF(-J)=YF(1)-SDYF*(J+1)
          YF(N2F+1+J)=YF(N2F)+SDYF*(J+1)
        ENDDO
        DO J=0,2*N2F_BD
          YPF(-J)=YF(-J)+0.5*SDYF
          YPF(N2F+J)=YF(N2F+J)+0.5*SDYF
        ENDDO

        !Z-DIR
        DO K=1,N3FM
        	ZPF(K)=ZF(K)+0.5*SDZF
        ENDDO
        DO K=0,2*N3F_BD
          ZF(-K)=ZF(1)-SDZF*(K+1)
          ZF(N3F+1+K)=ZF(N3F)+SDZF*(K+1)
        ENDDO
        DO K=0,2*N3F_BD
          ZPF(-K)=ZF(-K)+0.5*SDZF
          ZPF(N3F+K)=ZF(N3F+K)+0.5*SDZF
        ENDDO

      RETURN
      END
      
!******************************************************************
       SUBROUTINE INDICES
!******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      !NAVIER-STOKE
!-----STREAMWISE DIRECTION
      DO 10 I=0,N1
      IPV(I)=I+1
      IMV(I)=I-1
   10 CONTINUE
      IPV(N1+1)=N1+2
      IMV(-1)=-2
      
!-----NORMAL DIRECTION
      DO 20 J=0,N2
      JPV(J)=J+1
      JMV(J)=J-1
   20 CONTINUE
      JPV(N2+1)=N2+2
      JMV(-1)=-2
      
!-----SPANWISE DIRECTION
      DO 30 K=0,N3
      KPV(K)=K+1
      KMV(K)=K-1
   30 CONTINUE
      KPV(N3+1)=N3+2
      KMV(-1)=-2
      
      !FOR BOUNDARY CONDITION.
      DO 40 I=1,N1M
      FIXIL(I)=0.
      FIXIU(I)=0.
   40 CONTINUE
      FIXIL(1)=1.
      FIXIU(N1M)=1.

      DO 50 J=1,N2M
      FIXJL(J)=0.
      FIXJU(J)=0.
   50 CONTINUE
      FIXJL(1)=1.
      FIXJU(N2M)=1.

      DO 60 K=1,N3M
      FIXKL(K)=0.
      FIXKU(K)=0.
   60 CONTINUE
      FIXKL(1)=1.
      FIXKU(N3M)=1.


      RETURN
      END

!******************************************************************
       SUBROUTINE INDICES_LVS
!******************************************************************
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_GEOM_VAR
      
!------ Axial direction (X)
        DO 40 I=0,N1F
        IPF(I)=I+1
        IMF(I)=I-1
  40    CONTINUE
        IPF(N1F+1)=N1F+2 !FOR DIRICHLET AND NEUMANN CONDTIION IN TRANSPORT, REINITIALIZATION.
        IMF(-1)=-2

!------ Radial direction (R)
        DO 50 J=0,N2F
        JPF(J)=J+1
        JMF(J)=J-1
  50    CONTINUE
        JPF(N2F+1)=N2F+2
        JMF(-1)=-2
        
!------ Azimuthal direction (T)
        DO 60 K=0,N3F
        KPF(K)=K+1
        KMF(K)=K-1
  60    CONTINUE
        KPF(N3F+1)=N3F+2
        KMF(-1)=-2

      RETURN
      END

! =====================================================================
       SUBROUTINE GEOM_COUPLING
! =====================================================================
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_GEOM_VAR
      USE LVS_COUPLING
      
      !FOR FIND FLOW-SOLVER GRID USING REFINED GRID
      !FOR VELOCITY COUPLING

      !X-DIR
      IF ( N1M .EQ. 1 ) THEN
      DO I=1,N1FM
        ICOU_VEL(I)=1
        ICOUMP_VEL(I)=1
       ENDDO
      ELSE

      DO I=1,N1F
         XX=XF(I)
       DO II=1,N1M
        IF ( ABS(XX-XP(II)) .LE. 1.E-10 ) THEN
         ICOU_VEL(I)=II
         GOTO 11
        ELSE IF ( XP(1) .GT. XX ) THEN
         ICOU_VEL(I)=1
         GOTO 11
        ELSE IF ( XP(II) .GT. XX )  THEN
          ICOU_VEL(I)=II-1
         GOTO 11
        ELSE IF ( XP(N1M) .LT. XX ) THEN
         ICOU_VEL(I)=N1M
         GOTO 11
         ENDIF
       ENDDO

          WRITE(*,*) 'X-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 11   ENDDO
 
      DO I=1,N1F
         XX=XPF(I)
       DO II=1,N1M
        IF ( ABS(XX-XP(II)) .LE. 1.E-10 ) THEN
         ICOUMP_VEL(I)=II
         GOTO 10
        ELSE IF ( XP(1) .GT. XX ) THEN
         ICOUMP_VEL(I)=1
         GOTO 10
        ELSE IF ( XP(II) .GT. XX )  THEN
          ICOUMP_VEL(I)=II-1
         GOTO 10
        ELSE IF ( XP(N1M) .LT. XX ) THEN
         ICOUMP_VEL(I)=N1M
         GOTO 10
         ENDIF
       ENDDO

          WRITE(*,*) 'X-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 10   ENDDO
      ENDIF
      
      !Y-DIR
      DO J=1,N2F
         YY=YF(J)
       DO JJ=1,N2M
        IF ( ABS(YY-YP(JJ)) .LE. 1.E-10 ) THEN
         JCOU_VEL(J)=JJ
         GOTO 21
        ELSE IF ( YP(1) .GT. YY ) THEN
         JCOU_VEL(J)=1
         GOTO 21
        ELSE IF ( YP(JJ) .GT. YY ) THEN
          JCOU_VEL(J)=JJ-1
         GOTO 21
        ELSE IF ( YP(N2M) .LT. YY ) THEN
          JCOU_VEL(J)=N2M
           GOTO 21
         ENDIF
       ENDDO

          WRITE(*,*) 'Y-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 21   ENDDO

      DO J=1,N2F
         YY=YPF(J)
       DO JJ=1,N2M
        IF ( ABS(YY-YP(JJ)) .LE. 1.E-10 ) THEN
         JCOUMP_VEL(J)=JJ
         GOTO 20
        ELSE IF ( YP(1) .GT. YY ) THEN
         JCOUMP_VEL(J)=1
         GOTO 20
        ELSE IF ( YP(JJ) .GT. YY ) THEN
          JCOUMP_VEL(J)=JJ-1
         GOTO 20
        ELSE IF ( YP(N2M) .LT. YY ) THEN
           JCOUMP_VEL(J)=N2M
           GOTO 20
         ENDIF
       ENDDO

          WRITE(*,*) 'Y-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 20   ENDDO
 
      !Z-DIR
      IF ( N3M .EQ. 1 ) THEN
      DO K=1,N3F
        KCOU_VEL(K)=1
        KCOUMP_VEL(K)=1
       ENDDO
      ELSE
      DO K=1,N3F
         ZZ=ZF(K)
       DO KK=1,N3M
        IF ( ABS(ZZ-ZP(KK)) .LE. 1.E-10 ) THEN
         KCOU_VEL(K)=KK
         GOTO 31
        ELSE IF ( ZP(1) .GT. ZZ ) THEN
         KCOU_VEL(K)=1
         GOTO 31
        ELSE IF ( ZP(KK) .GT. ZZ ) THEN
          KCOU_VEL(K)=KK-1
         GOTO 31
        ELSE IF ( ZP(N3M) .LT. ZZ ) THEN
         KCOU_VEL(K)=N3M
         GOTO 31
         ENDIF
       ENDDO

          WRITE(*,*) 'Z-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 31   ENDDO
 
      DO K=1,N3F
         ZZ=ZPF(K)
       DO KK=1,N3M
        IF ( ABS(ZZ-ZP(KK)) .LE. 1.E-10 ) THEN
         KCOUMP_VEL(K)=KK
         GOTO 30
        ELSE IF ( ZP(1) .GT. ZZ ) THEN
         KCOUMP_VEL(K)=1
         GOTO 30
        ELSE IF ( ZP(KK) .GT. ZZ ) THEN
          KCOUMP_VEL(K)=KK-1
         GOTO 30
        ELSE IF ( ZP(N3M) .LT. ZZ ) THEN
         KCOUMP_VEL(K)=N3M
         GOTO 30
         ENDIF
       ENDDO

          WRITE(*,*) 'Z-DIR FLOW_SOLVER->REFINED COUPLING PROBLEM!'
 30   ENDDO
      ENDIF
 
      
      !X-DIR
      DO I=1,N1
       X_CRI=X(I)
        DO IF=1,N1F
         IF ( ABS(XF(IF)-X_CRI) .LE. 1.E-10 ) THEN
          ICOU1(I)=IF-1
          ICOU2(I)=IF
           GOTO 100
         ELSE IF( XF(IF) .GT. X_CRI ) THEN
          ICOU1(I)=IF-1
          ICOU2(I)=IF-1
          GOTO 100
         ENDIF
        ENDDO
           WRITE(*,*) 'CAL ICOU HAVE PROBLEM!!'
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 100  ENDDO
         ICOU1(0)=ICOU1(1)-(ICOU1(N1M)-ICOU1(N1M-1))
         ICOU1(-1)=ICOU1(0)-(ICOU1(N1M-1)-ICOU1(N1M-2))
         ICOU1(N1+1)=ICOU1(N1)+(ICOU1(2)-ICOU1(1))
         ICOU1(N1+2)=ICOU1(N1+1)+(ICOU1(3)-ICOU1(2))
         ICOU2(0)=ICOU2(1)-(ICOU2(N1M)-ICOU2(N1M-1))
         ICOU2(-1)=ICOU2(0)-(ICOU2(N1M-1)-ICOU2(N1M-2))
         ICOU2(N1+1)=ICOU2(N1)+(ICOU2(2)-ICOU2(1))
         ICOU2(N1+2)=ICOU2(N1+1)+(ICOU2(3)-ICOU2(2))

      DO I=1,N1M
       X_CRI=XP(I)
        DO IF=1,N1F
        IF ( ABS(XF(IF)-X_CRI) .LE. 1.E-10 ) THEN
          ICOUMP1(I)=IF-1
          ICOUMP2(I)=IF
           GOTO 101
         ELSE IF( XF(IF) .GT. X_CRI ) THEN
          ICOUMP1(I)=IF-1
          ICOUMP2(I)=IF-1
          GOTO 101
         ENDIF
        ENDDO
           WRITE(*,*) 'CAL ICOUMP HAVE PROBLEM!!'
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 101  ENDDO
         ICOUMP1(0)=ICOUMP1(1)-(ICOUMP1(N1M)-ICOUMP1(N1M-1))
         ICOUMP1(-1)=ICOUMP1(0)-(ICOUMP1(N1M-1)-ICOUMP1(N1M-2))
         ICOUMP1(N1)=ICOUMP1(N1M)+(ICOUMP1(N1M)-ICOUMP1(N1M-1)) !APPROXIMATION
         ICOUMP1(N1+1)=ICOUMP1(N1)+(ICOUMP1(2)-ICOUMP1(1))
         ICOUMP2(0)=ICOUMP2(1)-(ICOUMP2(N1M)-ICOUMP2(N1M-1))
         ICOUMP2(-1)=ICOUMP2(0)-(ICOUMP2(N1M-1)-ICOUMP2(N1M-2))
         ICOUMP2(N1)=ICOUMP2(N1M)+(ICOUMP2(N1M)-ICOUMP2(N1M-1)) !APPROXIMATION
         ICOUMP2(N1+1)=ICOUMP2(N1)+(ICOUMP2(2)-ICOUMP2(1))

      !Y-DIR
      DO J=1,N2
       Y_CRI=Y(J)
        DO JF=1,N2F
        IF ( ABS(YF(JF)-Y_CRI) .LE. 1.E-10 ) THEN
          JCOU1(J)=JF-1
          JCOU2(J)=JF
           GOTO 200
         ELSE IF( YF(JF) .GT. Y_CRI ) THEN
          JCOU1(J)=JF-1
          JCOU2(J)=JF-1
          GOTO 200
         ENDIF
        ENDDO
           WRITE(*,*) 'CAL JCOU HAVE PROBLEM!!'
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 200  ENDDO
         JCOU1(0)=JCOU1(1)-(JCOU1(N2M)-JCOU1(N2M-1))
         JCOU1(-1)=JCOU1(0)-(JCOU1(N2M-1)-JCOU1(N2M-2))
         JCOU1(N2+1)=JCOU1(N2)+(JCOU1(2)-JCOU1(1))
         JCOU1(N2+2)=JCOU1(N2+1)+(JCOU1(3)-JCOU1(2))
         JCOU2(0)=JCOU2(1)-(JCOU2(N2M)-JCOU2(N2M-1))
         JCOU2(-1)=JCOU2(0)-(JCOU2(N2M-1)-JCOU2(N2M-2))
         JCOU2(N2+1)=JCOU2(N2)+(JCOU2(2)-JCOU2(1))
         JCOU2(N2+2)=JCOU2(N2+1)+(JCOU2(3)-JCOU2(2))

      DO J=1,N2M
       Y_CRI=YP(J)
        DO JF=1,N2F
        IF ( ABS(YF(JF)-Y_CRI) .LE. 1.E-10 ) THEN
          JCOUMP1(J)=JF-1
          JCOUMP2(J)=JF
           GOTO 201
         ELSE IF( YF(JF) .GT. Y_CRI ) THEN
          JCOUMP1(J)=JF-1
          JCOUMP2(J)=JF-1
          GOTO 201

         ENDIF
        ENDDO
          WRITE(*,*) 'CAL JCOUMP HAVE PROBLEM!!'
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 201  ENDDO
         JCOUMP1(0)=JCOUMP1(1)-(JCOUMP1(N2M)-JCOUMP1(N2M-1))
         JCOUMP1(-1)=JCOUMP1(0)-(JCOUMP1(N2M-1)-JCOUMP1(N2M-2))
         JCOUMP1(N2)=JCOUMP1(N2M)+(JCOUMP1(N2M)-JCOUMP1(N2M-1)) !APPROXIMATION
         JCOUMP1(N2+1)=JCOUMP1(N2)+(JCOUMP1(2)-JCOUMP1(1))
         JCOUMP2(0)=JCOUMP2(1)-(JCOUMP2(N2M)-JCOUMP2(N2M-1))
         JCOUMP2(-1)=JCOUMP2(0)-(JCOUMP2(N2M-1)-JCOUMP2(N2M-2))
         JCOUMP2(N2)=JCOUMP2(N2M)+(JCOUMP2(N2M)-JCOUMP2(N2M-1)) !APPROXIMATION
         JCOUMP2(N2+1)=JCOUMP2(N2)+(JCOUMP2(2)-JCOUMP2(1))
         
       !Z-DIR
      DO 300 K=1,N3
       Z_CRI=Z(K)
        DO KF=1,N3F
        IF ( ABS(ZF(KF)-Z_CRI) .LE. 1.E-10 ) THEN
          KCOU1(K)=KF-1
          KCOU2(K)=KF
           GOTO 300
         ELSE IF( ZF(KF) .GT. Z_CRI) THEN
          KCOU1(K)=KF-1
          KCOU2(K)=KF-1
          GOTO 300
         ENDIF
        ENDDO
           WRITE(*,*) 'CAL KCOU HAVE PROBLEM!!', K
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 300  ENDDO
         KCOU1(0)=KCOU1(1)-(KCOU1(N3M)-KCOU1(N3M-1))
         KCOU1(-1)=KCOU1(0)-(KCOU1(N3M-1)-KCOU1(N3M-2))
         KCOU1(N3+1)=KCOU1(N3)+(KCOU1(2)-KCOU1(1))
         KCOU1(N3+2)=KCOU1(N3+1)+(KCOU1(3)-KCOU1(2))
         KCOU2(0)=KCOU2(1)-(KCOU2(N3M)-KCOU2(N3M-1))
         KCOU2(-1)=KCOU2(0)-(KCOU2(N3M-1)-KCOU2(N3M-2))
         KCOU2(N3+1)=KCOU2(N3)+(KCOU2(2)-KCOU2(1))
         KCOU2(N3+2)=KCOU2(N3+1)+(KCOU2(3)-KCOU2(2))
         
      DO 301 K=1,N3M
       Z_CRI=ZP(K)
        DO KF=1,N3F
        IF ( ABS(ZF(KF)-Z_CRI) .LE. 1.E-10 ) THEN
          KCOUMP1(K)=KF-1
          KCOUMP2(K)=KF
           GOTO 301
         ELSE IF( ZF(KF) .GT. Z_CRI) THEN
          KCOUMP1(K)=KF-1
          KCOUMP2(K)=KF-1
          GOTO 301
         ENDIF
        ENDDO
           WRITE(*,*) 'CAL KCOUMP HAVE PROBLEM!!'
      !IF YOU RUN FLOWFIELD.F, THIS IS CAUSED BY THE 
      !DIFFERENCE OF "N1F_BD, N2F_BD, N3F_BD" OF TWO FIELDS
      !%IF YOU USE SPARE GRID FOR SIMULATE PIPE USING IBM,
      !SIMULATION COMES HERE DEPENDING ON THE NUMBER OF REFINEMENTS.
 301  ENDDO
         KCOUMP1(0)=KCOUMP1(1)-(KCOUMP1(N3M)-KCOUMP1(N3M-1))
         KCOUMP1(-1)=KCOUMP1(0)-(KCOUMP1(N3M-1)-KCOUMP1(N3M-2))
         KCOUMP1(N3)=KCOUMP1(N3M)+(KCOUMP1(N3M)-KCOUMP1(N3M-1)) !APPROXIMATION
         KCOUMP1(N3+1)=KCOUMP1(N3)+(KCOUMP1(2)-KCOUMP1(1))
         KCOUMP2(0)=KCOUMP2(1)-(KCOUMP2(N3M)-KCOUMP2(N3M-1))
         KCOUMP2(-1)=KCOUMP2(0)-(KCOUMP2(N3M-1)-KCOUMP2(N3M-2))
         KCOUMP2(N3)=KCOUMP2(N3M)+(KCOUMP2(N3M)-KCOUMP2(N3M-1)) !APPROXIMATION
         KCOUMP2(N3+1)=KCOUMP2(N3)+(KCOUMP2(2)-KCOUMP2(1))
      
!      DO I=1,N2
!       WRITE(*,678) I,ICOU1(I),JCOU1(I),KCOU1(I)
!      ENDDO
!      DO I=1,N2
!       WRITE(*,678) I,ICOUMP1(I),JCOUMP1(I),KCOUMP1(I)
!      ENDDO
!      DO I=1,N2
!       WRITE(*,678) I,ICOUMP2(I),JCOUMP2(I),KCOUMP2(I)
!      ENDDO
! 678   FORMAT(4I5)  
 
        RETURN
        END   