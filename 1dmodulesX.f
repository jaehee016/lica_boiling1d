!-----------------------------------------------------------------------
       MODULE PARAM_VAR
!-----------------------------------------------------------------------             
 !-----GEOMETRY
       PARAMETER (M1 =2 ,M2 = 2,M3 =257)
       PARAMETER (M1M=M1-1,M2M=M2-1,M3M=M3-1)
 !-----MULTIGRID
       PARAMETER (M1MD=2*M1M)
       PARAMETER (M3MD=2*M3M)
       PARAMETER (M13MD=2*M1M*M3M)
       PARAMETER (MLEV=20)
       PARAMETER (M1MH=M1M/2+1)
       PARAMETER (M3MH=M3M/2+1)
 !-----HISTORY
       PARAMETER (MHIST=1000000)     !# of time steps
 !-----IMMERSED BOUNDARY
       PARAMETER (MBODY0=   273408)
       PARAMETER (MBODY =   6295296)

      PARAMETER (MLVS=1)  !THE NUMBER OF LVS
      PARAMETER (M_MAX=11) !MAXIMUM BAND NUMBER
      PARAMETER (M1F=3,M2F=3,M3F=513)
      PARAMETER (M1FM=M1F-1,M2FM=M2F-1,M3FM=M3F-1)
      PARAMETER(M1F_BD=FLOAT(M1FM)/FLOAT(M1M)
     &,M2F_BD=FLOAT(M2FM)/FLOAT(M2M),M3F_BD=FLOAT(M3FM)/FLOAT(M3M))
      PARAMETER (MF=M1FM*M2FM*M3FM)
      PARAMETER (MF_BAND=0.1*MF)

      PARAMETER(M1L=-2*M1F_BD,M1U=M1F+2*M1F_BD+1)
      PARAMETER(M2L=-2*M2F_BD,M2U=M2F+2*M2F_BD+1)
      PARAMETER(M3L=-2*M3F_BD,M3U=M3F+2*M3F_BD+1)
      END MODULE PARAM_VAR

      
!-----------------------------------------------------------------------      
      MODULE FLOW_VAR
!-----------------------------------------------------------------------      
        REAL*8    ::    RE,REI,REM,REP
        INTEGER*8 ::    M,MSUB
        REAL*8    ::    DT,TIME
        REAL*8    ::    DTCONST,DTCONSTI
        REAL*8    ::    ACOEF,ACOEFI
        REAL*8    ::    CFLMAX
        INTEGER*8 ::    NTST,NPIN,NPRINT,NTIME
        INTEGER*8 ::    IREAD,IWRT,IRESET,IPZERO
        INTEGER*8 ::    ITRACE,NTEST,NTEST2,NTPRINT
        INTEGER*8 ::    IDTOPT,IDTOLD,IMPL
        INTEGER*8 ::    ITOPEN
        REAL*8    ::    ALPHA_RK3,GAMMA(3),RO(3),UFC
        INTEGER*8 ::    ICH
        INTEGER*8 ::    IPOISS
        REAL*8    ::    PMI,PMI_CONST,PMI_DUDY,PMI_GRAV
        
        REAL*8, ALLOCATABLE :: SUR_X(:,:,:),SUR_Y(:,:,:),SUR_Z(:,:,:)
        REAL*8, ALLOCATABLE :: UO(:,:,:),VO(:,:,:),WO(:,:,:)
        REAL*8, ALLOCATABLE :: AIU(:),CIU(:),AIVW(:),CIVW(:)
     &                      ,AJV(:),CJV(:),AJUW(:),CJUW(:)
     &                      ,AKW(:),CKW(:),AKUV(:),CKUV(:)
      END MODULE FLOW_VAR


!-----------------------------------------------------------------------
      MODULE FLOW_GEOM_VAR
!-----------------------------------------------------------------------      
      INTEGER*8 ::      N1,N1M,N2,N2M,N3,N3M
      INTEGER*8 ::      IPX,IPY,IPZ,IBG,JBG,KBG
      INTEGER*8, ALLOCATABLE ::IPV(:),JPV(:),KPV(:),IMV(:),JMV(:),KMV(:)
      
      REAL*8 :: XL,YL,ZL      
      REAL*8, ALLOCATABLE ::  FIXIL(:),FIXIU(:),
     &                        FIXJL(:),FIXJU(:),
     &                        FIXKL(:),FIXKU(:)
                      
      REAL*8, ALLOCATABLE ::  X(:),Y(:),Z(:)
      REAL*8, ALLOCATABLE ::  XP(:),YP(:),ZP(:)
      REAL*8, ALLOCATABLE ::  SDX(:),SDY(:),SDZ(:)
      REAL*8, ALLOCATABLE ::  VDX(:),VDY(:),VDZ(:)
      REAL*8, ALLOCATABLE ::  SSDX(:),SSDY(:),SSDZ(:)
      REAL*8, ALLOCATABLE ::  VVDX(:),VVDY(:),VVDZ(:)
      END MODULE FLOW_GEOM_VAR

     
!-----------------------------------------------------------------------      
      MODULE LVS_VAR
!-----------------------------------------------------------------------      
      REAL*8,     ALLOCATABLE :: ALPHI(:,:,:)
      REAL*8,     ALLOCATABLE :: ALPHI_EXT(:,:,:,:)
      INTEGER*8,  ALLOCATABLE :: NUMA(:,:)
      INTEGER*8,  ALLOCATABLE :: I_B(:,:,:),J_B(:,:,:),K_B(:,:,:)
      REAL*8,     ALLOCATABLE :: PSIF(:,:,:)
      REAL*8 :: PSIF_BACKGROUND
      
      INTEGER*8 :: ITRACKING,MCLS,MGLOBAL    
      INTEGER*8 :: N_MAX                     
      INTEGER*8 :: NLVS
      INTEGER*8 :: NUMG,NUMA_MAX
     
      INTEGER*8, ALLOCATABLE :: MASK_BUB(:,:,:)
      INTEGER*8, ALLOCATABLE :: MASK_GLOBAL(:,:,:)
      
      INTEGER*8, ALLOCATABLE :: CONTACT(:,:),CONTACT_OLD(:,:)        !2017-08-16
      REAL*8   , ALLOCATABLE :: CONTACT_START(:,:),CONTACT_TIME(:,:) !2017-08-16
      INTEGER*8, ALLOCATABLE :: MASK_COL(:,:,:)                      !2017-08-20
      REAL*8   , ALLOCATABLE :: ALPHI_COL(:,:,:)                     !2017-08-20
      REAL*8   , ALLOCATABLE :: VELVS(:,:)                           !2017-09-22_T_CON                  
      REAL*8   , ALLOCATABLE :: T_CON(:,:)                           !2017-09-22_T_CON                  
      INTEGER*8, ALLOCATABLE :: LVSON(:)
      REAL*8   , ALLOCATABLE :: XLAGPRT(:,:),ULAGPRT(:,:),VLAGPRT(:,:) ! U: FLOW, V: PARTICLE
      REAL*8   , ALLOCATABLE :: DVDT(:,:),DUDX(:,:,:)
      INTEGER*8 :: NPRT,LAGHIST,ILAG
      END MODULE LVS_VAR

      
!-----------------------------------------------------------------------        
      MODULE LVS_GEOM_VAR
!-----------------------------------------------------------------------        
      INTEGER*8 :: N1F,N1FM,N2F,N2FM,N3F,N3FM,N3FH,N1F_BD,N2F_BD,N3F_BD
      REAL*8, ALLOCATABLE :: XF(:),YF(:),ZF(:)
      REAL*8, ALLOCATABLE :: XPF(:),YPF(:),ZPF(:)
      INTEGER*8, ALLOCATABLE :: IPF(:),IMF(:),
     &                          JPF(:),JMF(:),
     &                          KPF(:),KMF(:)
      REAL*8 :: SDXF,SDYF,SDZF,SSDXF,SSDYF,SSDZF
      END MODULE LVS_GEOM_VAR
       
!-----------------------------------------------------------------------        
      MODULE LVS_COUPLING
!-----------------------------------------------------------------------        
      INTEGER*8,ALLOCATABLE :: ICOU_VEL(:),JCOU_VEL(:),KCOU_VEL(:)
      INTEGER*8,ALLOCATABLE :: ICOUMP_VEL(:),JCOUMP_VEL(:),KCOUMP_VEL(:)
      INTEGER*8,ALLOCATABLE :: ICOU1(:),JCOU1(:),KCOU1(:),
     &                         ICOUMP1(:),JCOUMP1(:),KCOUMP1(:)
      INTEGER*8,ALLOCATABLE :: ICOU2(:),JCOU2(:),KCOU2(:),
     &                         ICOUMP2(:),JCOUMP2(:),KCOUMP2(:)
      END MODULE LVS_COUPLING

      
!-----------------------------------------------------------------------        
      MODULE TWO_PHASE_PROPERTY
!-----------------------------------------------------------------------        
        REAL*8 :: DENM,DENP,DENR,DEN_DIFF
        REAL*8 :: VISM,VISP,VISR,VIS_DIFF
        REAL*8 :: SCM,SCP,SCR,SC_DIFF,DENSCM,DENSCP,DENSC_DIFF !SPECIFIC CAPACITY
        REAL*8 :: TCM,TCP,TCR,TC_DIFF !THERMAL CONDUCTIVITY
        REAL*8 :: PRM,PRP
        REAL*8 :: SURF,GRAVITY
        REAL*8 :: FR,SURF_J,FK,RE_AIR
        REAL*8 :: TIME_ST,TIME_GV
        
      END MODULE TWO_PHASE_PROPERTY

!-----------------------------------------------------------------------        
      MODULE MG_OPTION
!-----------------------------------------------------------------------        
       INTEGER*8 :: NLEV,IWC,NBLI,MGITR,IOLDV,IMGSOR
       REAL*8 :: RESID_POI,WWSOR
       INTEGER*8 :: IPCG
      END MODULE MG_OPTION

      
!-----------------------------------------------------------------------  
      MODULE IBM_VAR
!-----------------------------------------------------------------------        
      INTEGER*8 IBMON
      INTEGER*8 MASSON
      REAL*8, ALLOCATABLE :: QMASS(:,:,:)
      
      INTEGER*8, ALLOCATABLE :: NFC_INTP(:),NFC_INNER(:)
      INTEGER*8, ALLOCATABLE :: IFC(:,:),JFC(:,:),KFC(:,:)
      REAL*8, ALLOCATABLE :: FCV(:,:),FCVAVG(:,:)
      
      INTEGER*8, ALLOCATABLE :: MPI(:,:,:)
      REAL*8, ALLOCATABLE :: GFI(:,:,:,:,:)
      REAL*8, ALLOCATABLE :: FORCESUM(:)
      
      
      INTEGER*8 :: NFC_INTP_PS,NFC_INNER_PS
      INTEGER*8, ALLOCATABLE :: IFC_PS(:),JFC_PS(:),KFC_PS(:)
      REAL*8, ALLOCATABLE :: FCV_PS(:),FCVAVG_PS(:)

      INTEGER*8, ALLOCATABLE :: MPI_PS(:,:)
      REAL*8, ALLOCATABLE :: GFI_PS(:,:,:,:)
      REAL*8, ALLOCATABLE :: RHSPS_IBM(:,:,:)
      END MODULE IBM_VAR

!-----------------------------------------------------------------------        
      MODULE HEAT_VAR
!-----------------------------------------------------------------------        
      INTEGER*8 :: IPS
      REAL*8 :: PRA,PRAI
      REAL*8 :: HEAT_SUM
      REAL*8, ALLOCATABLE :: T(:,:,:)
      REAL*8, ALLOCATABLE :: TALPH(:,:,:)
      REAL*8, ALLOCATABLE :: RHSPS(:,:,:)
      END MODULE HEAT_VAR

!-----------------------------------------------------------------------  
      MODULE FLD_AVG
!-----------------------------------------------------------------------       
      INTEGER*8 :: IAVG,NPRIAVG
      REAL*8, ALLOCATABLE :: URAVGZ(:,:),UTAVGZ(:,:) !RT
      REAL*8, ALLOCATABLE :: UAVGZ(:,:),VAVGZ(:,:)   !XYZ
      REAL*8, ALLOCATABLE :: WAVGZ(:,:),UIUJAVGZ(:,:,:),
     &                     PAVGZ(:,:),P2AVGZ(:,:),
     &                     VOLFZ(:,:) !TNUAVGZ(:,:),
!     &                     SIJAVGZ(:,:,:),VORAVGZ(:,:,:),
!     &                     SSAVGZ(:,:),SGSDISSZ(:,:)
      REAL*8 :: TIMEINIT,TIMEEND
      INTEGER*8 :: IHISTINIT,IHISTEND
      REAL*8 :: PMIAVG

      REAL*8, ALLOCATABLE :: URAVGZ0(:,:),UTAVGZ0(:,:) !RT
      REAL*8, ALLOCATABLE :: UAVGZ0(:,:),VAVGZ0(:,:)   !XYZ
      REAL*8, ALLOCATABLE :: WAVGZ0(:,:),UIUJAVGZ0(:,:,:),
     &                     PAVGZ0(:,:),P2AVGZ0(:,:),TIME0(:,:)

      REAL*8, ALLOCATABLE :: URAVGZ1(:,:),UTAVGZ1(:,:) !RT
      REAL*8, ALLOCATABLE :: UAVGZ1(:,:),VAVGZ1(:,:)   !XYZ
      REAL*8, ALLOCATABLE :: WAVGZ1(:,:),UIUJAVGZ1(:,:,:),
     &                     PAVGZ1(:,:),P2AVGZ1(:,:),TIME1(:,:)

      REAL*8, ALLOCATABLE :: PSAVGZ(:,:),PSRMSZ(:,:),
     &                    PSUIAVGZ(:,:,:),TALPHAVGZ(:,:)
      REAL*8, ALLOCATABLE :: PSAVGZ0(:,:),PSRMSZ0(:,:),PSUIAVGZ0(:,:,:)
      REAL*8, ALLOCATABLE :: PSAVGZ1(:,:),PSRMSZ1(:,:),PSUIAVGZ1(:,:,:)
      END MODULE FLD_AVG

!-----------------------------------------------------------------------       
      MODULE TIME_VAR
!-----------------------------------------------------------------------       
      REAL*8 TIME_BEGIN,TIME_END
      REAL*8 TOTAL_TIME_B,TOTAL_TIME_E      
      REAL*8, ALLOCATABLE :: SGSTIME_B(:),SGSTIME_E(:)
     &   ,ANSETIME_B(:),ANSETIME_E(:),POISSTIME_B(:),POISSTIME_E(:)
     &   ,CONTINUITY_B(:),CONTINUITY_E(:),ALVS_B(:),ALVS_E(:)

      END MODULE TIME_VAR

!-----------------------------------------------------------------------  
      MODULE LES_VAR
!-----------------------------------------------------------------------        
      INTEGER*8 ILES,ISGSMOD,IFILTER
      REAL*8 CSGS,CSGSPS
      LOGICAL DYNMON
      
      REAL*8, ALLOCATABLE :: CFX1 (:,:),CFX2 (:,:),CFZP1(:,:),CFZP2(:,:)
      REAL*8, ALLOCATABLE :: TNU(:,:,:)
      REAL*8, ALLOCATABLE :: TNUX(:,:,:),TNUY(:,:,:),TNUZ(:,:,:)
      END MODULE LES_VAR                              
!      module mon_var
!      real*8                       :: TIME_BEGIN,TIME_END,SGSTIME_B(3),SGSTIME_E(3)
!      real*8                       :: RHSTIME_B(3),RHSTIME_E(3),LHSTIME_B(3),LHSTIME_E(3)
!      real*8                       :: POISSTIME_B(3),POISSTIME_E(3)
!      end module mon_var
!-----------------------------------------------------------------------        
      MODULE PHASE_BOILING
!-----------------------------------------------------------------------        
      INTEGER*8  :: IPHS
      REAL*8 :: JA_AIR,T_SAT
      INTEGER*8, ALLOCATABLE :: IMASK_PH(:,:,:)
      REAL*8, ALLOCATABLE :: TEMP(:,:,:,:)
      !REAL*8, ALLOCATABLE :: MFDOT(:,:,:)
      REAL*8, ALLOCATABLE :: SFGAS(:,:,:)
      REAL*8, ALLOCATABLE :: SFGAS2(:,:,:)
      REAL*8, ALLOCATABLE :: MFLVS(:,:,:)
      REAL*8, ALLOCATABLE :: SGLVS(:,:,:)
      REAL*8, ALLOCATABLE :: SGAS_X(:,:,:),SGAS_Y(:,:,:)
      REAL*8, ALLOCATABLE :: SGAS_Z(:,:,:)
      REAL*8, ALLOCATABLE :: temp_B(:),temp_E(:)
      REAL*8, ALLOCATABLE :: VOL_TOT_IPHS(:)

      END MODULE PHASE_BOILING
