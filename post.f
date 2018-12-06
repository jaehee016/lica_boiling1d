CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                     C
C     POST-PROCESSING CODE FOR INSTANTANEOUS FLOW FIELD DATA          C
C                                                                     C
C     Referred S. Kang and D. Kim's postprocessing codes              C
C     Referred WOONGJAE CHANG`s postprocessing codes  2002.11.20.     C
C                                                                     C 
C       TWO-PHASE POST-PROCESSING CODE   KIYOUNG KIM  2015. 7.        C 
C                                                                     C 
C---------------------------------------------------------------------C

      PROGRAM MAIN

       USE PARAM_VAR
       USE FLOW_VAR
       USE IBM_VAR
        
       USE FLOW_GEOM_VAR
        
       USE LVS_VAR

       USE HEAT_VAR
      
       REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
       REAL P(0:M1,0:M2,0:M3)

      REAL PSI_XN(0:M1,0:M2,0:M3)
      REAL PSI_YN(0:M1,0:M2,0:M3)
      REAL PSI_ZN(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)

      INTEGER BUBCOL(MLVS)
      INTEGER BUBCOL2(MLVS)
      integer nn2
      
      
       CHARACTER*20 DUMMY
       CHARACTER*30 gridfile, inputfile
       CHARACTER*3 tfn1

       COMMON/CENTER/VPW(0:M1,0:M2,0:M3,0:9)

       COMMON/SKIP1/ISKIP,JSKIP,KSKIP
       COMMON/SKIP2/XBIG,YBIG,ZBIG,XEND,YEND,ZEND

      !LVS_VAR
      ALLOCATE( ALPHI(MF_BAND,M_MAX,MLVS))
      ALLOCATE( ALPHI_EXT(-2:M1F+3,-2:M2F+3,-2:M3F+3,1))
      ALLOCATE( NUMA(M_MAX,MLVS))
      ALLOCATE( I_B(MF_BAND,M_MAX,MLVS),J_B(MF_BAND,M_MAX,MLVS)
     &                           ,K_B(MF_BAND,M_MAX,MLVS))
      ALLOCATE( PSIF(M1L:M1U,M2L:M2U,M3L:M3U))

      ALLOCATE(MASK_BUB(M1M,M2M,M3M))
      ALLOCATE(MASK_GLOBAL(0:M1F,0:M2F,0:M3F))
      
      ALLOCATE( ALPHI_COL(M1L:M1U,M2L:M2U,M3L:M3U))                 ! 2017-08-20
      ALLOCATE( MASK_COL(0:M1F,0:M2F,0:M3F))                        ! 2017-08-20
      ALLOCATE( CONTACT(MLVS,MLVS),CONTACT_OLD(MLVS,MLVS))          ! 2017-08-16 
      ALLOCATE( CONTACT_START(MLVS,MLVS),CONTACT_TIME(MLVS,MLVS))   ! 2017-08-16 
      ALLOCATE( VELVS(MLVS,3))                                      ! 2017-09-24
      ALLOCATE( T_CON(MLVS,MLVS))                                      ! 2017-09-24
      ALLOCATE( LVSON(MLVS))
      ALLOCATE( XLAGPRT(10000,4), uLAGPRT(10000,6), vlagprt(10000,3)) !20171222, 2018-01-07
      ALLOCATE( DVDT(10000,6),DUDX(6,10000,3))
                  
      !HEAT_VAR
      ALLOCATE(T(0:M1,0:M2,0:M3))
      
        OPEN(2,FILE='post.in')
        WRITE(*,87)
        READ(2,*) DUMMY
        READ(2,*) DUMMY
        READ(2,*) IMOVIE
        WRITE(*,*) IMOVIE
          IF ( IMOVIE .EQ. 0 ) THEN
        READ(2,*) DUMMY
        READ(2,*) ITIMESTEP
        WRITE(*,*) ITIMESTEP
        READ(2,*) DUMMY
        READ(2,*) DUMMY
         IADD=0
         IHOWMANY=0
          ELSE IF ( IMOVIE .EQ. 1 ) THEN
        READ(2,*) DUMMY
        READ(2,*) DUMMY
        READ(2,*) DUMMY
        READ(2,*) ITIMESTEP,IADD,IHOWMANY
        WRITE(*,*) ITIMESTEP,IADD,IHOWMANY
         ENDIF
        READ(2,*) DUMMY
        READ(2,*) IBMON
        READ(2,*) DUMMY
        READ(2,*) XBIG,XEND,YBIG,YEND,ZBIG,ZEND
        READ(2,*) DUMMY
        READ(2,*) ISKIP,JSKIP,KSKIP
        READ(2,*) DUMMY
        READ(2,*) IVEL_CYLIND,IVOR_CYLIND
        READ(2,*) DUMMY
        READ(2,*) DUMMY
        READ(2,*) IDIMEN
        READ(2,*) DUMMY
        READ(2,*) I_POS,J_POS,K_POS
        READ(2,*) DUMMY
        READ(2,*) LVSGRID
        READ(2,*) DUMMY
        READ(2,*) DUMMY
        WRITE(*,*)
        READ(2,'(a)') gridfile
        WRITE(*,98) gridfile

  87  FORMAT('----------- BASIC OPTION -----------')
      WRITE(*,100) IPX,IPY,IPZ
      WRITE(*,101) XBIG,XEND
      WRITE(*,102) YBIG,YEND
      WRITE(*,103) ZBIG,ZEND
      WRITE(*,104) ISKIP,JSKIP,KSKIP
      WRITE(*,105) IDIMEN
      WRITE(*,106) I_POS,J_POS,K_POS
      WRITE(*,107) IVEL_CYLIND,IVOR_CYLIND
      WRITE(*,108) LVSGRID
 100  FORMAT('IPX=',I5,'IPY=',I5,'IPZ=',I5)
 101  FORMAT('XBIG=',ES13.5,'XEND=',E13.5)
 102  FORMAT('YBIG=',ES13.5,'YEND=',E13.5)
 103  FORMAT('ZBIG=',ES13.5,'ZEND=',E13.5)
 104  FORMAT('ISKIP=',I5,'JSKIP=',I5,'KSKIP=',I5)
 105  FORMAT('IDIMEN=',I5)
 106  FORMAT('I_POS=',I5,'J_POS=',I5,'K_POS=',I5)
 107  FORMAT('IVEL_CYLIND=',I5,'IVOR_CYLIND=',I5)
 108  FORMAT('LVSGRID=',I5)

  88  FORMAT('----------- GRID INFORMATION -----------')
  94    FORMAT('IMOVIE = ',I5)
  95    FORMAT('ITIMESTEP = ',I5)
  96    FORMAT('ITIMESTEP = ',I5,' IADD = ',I5,' IHOWMANY = ',I8)
  98    FORMAT('GRIDFILE :       ',A30)


C-----INPUT PARAMETER such as IPX,IPY,IPZ,IPS
      tfn1='fld' 
      idg1=ITIMESTEP/100000
      idg2=(ITIMESTEP-idg1*100000)/10000
      idg3=(ITIMESTEP-idg1*100000-idg2*10000)/1000
      idg4=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6= ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
      inputfile=tfn1//char(idg1+48)//char(idg2+48)//char(idg3+48)
     &              //char(idg4+48)//char(idg5+48)//char(idg6+48)

      OPEN(12,FILE=inputfile)
      READ(12,1001) N1,N2,N3
      READ(12,1002) IHIST,M,TIME,DT
      READ(12,1003) IPS,IPX,IPY,IPZ,PRAA
 1001   FORMAT(3I8)
 1002   FORMAT(2I8,2ES20.12)
 1003   FORMAT(4I8,1ES20.12)
      CLOSE(12)
C-----INPUT PARAMETER

C********* GEOMETRY
        CALL GEOM(gridfile)
C********************

!        OPEN(201,FILE='0POST.DAT',POSITION='APPEND')
!      WRITE(201,333) 'VARIABLES="TIME","VOL1","VOL2","Q1","Q2","Q_TOT"
!     &,"BULK1","BULK2"'
!  333   FORMAT(A100)
!        CLOSE(201)

        IF ( IMOVIE .EQ. 0 ) IHOWMANY=1    !SNAP SHOT

      DO 5000 IFLDNUM=1,IHOWMANY

c-----field read
      tfn1='fld' 
      idg1=ITIMESTEP/100000
      idg2=(ITIMESTEP-idg1*100000)/10000
      idg3=(ITIMESTEP-idg1*100000-idg2*10000)/1000
      idg4=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6= ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
      inputfile=tfn1//char(idg1+48)//char(idg2+48)//char(idg3+48)
     &              //char(idg4+48)//char(idg5+48)//char(idg6+48)

      WRITE(*,89)
  89  FORMAT('----------- FIELD INFORMATION -----------')
      WRITE(*,90) inputfile
   90 FORMAT('FIELD NAME:       ',A30)

      OPEN(12,FILE=inputfile)
      CALL READFLD(12,U,V,W,P)
      CLOSE(12)
c-----field read
      nn2=n_max-1
      DO LLVS=1,NLVS

       DO N=1,N_MAX
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         ALPHI_EXT(I,J,K,1)=ALPHI(NN,N,LLVS)
        ENDDO
       ENDDO

       CALL CAL_PSIF(1,NN2-1,LLVS,U,V,W,PSI_CN,BUBCOL,bubcol2
     &,vol_tot_ori )
       CALL GRID_COUPLING(1,PSI_XN,PSI_YN,PSI_ZN,PSI_CN)
      ENDDO
       IF (NLVS .NE. 0) DEALLOCATE(SUR_X,SUR_Y,SUR_Z)
        CALL P_FLUC(P,PSI_CN)
        CALL CENTERVEL(U,V,W,P)

       IF ( IDIMEN .EQ. 2 ) THEN

       IF (LVSGRID .EQ. 1) CALL VORNLAMBDA2
       CALL XCUT(I_POS,LVSGRID,ITIMESTEP,PSI_CN,IVEL_CYLIND,IVOR_CYLIND)
       CALL YCUT(J_POS,LVSGRID,ITIMESTEP,PSI_CN,IVEL_CYLIND,IVOR_CYLIND)
       CALL ZCUT(K_POS,LVSGRID,ITIMESTEP,PSI_CN,IVEL_CYLIND,IVOR_CYLIND)
       CALL CAL_VOL(W,PSI_CN)

       ELSE IF ( IDIMEN .EQ. 3 ) THEN
        CALL THREE_DIMENSION(LVSGRID,ITIMESTEP,PSI_CN,U,V,W
     &,IVEL_CYLIND,IVOR_CYLIND)
        CALL CAL_VOL(W,PSI_CN)
       ENDIF   !IDIMEN ENDIF

 4900  ITIMESTEP=ITIMESTEP+IADD
 5000  CONTINUE    !IMOVIE ENDDO 

       CLOSE(201)  
       CLOSE(146)  

      STOP
      END

C============================================================================
      SUBROUTINE READFLD(NV,U,V,W,P)          ! FIELD READING PART (BINARY FORMAT)
C============================================================================
      USE FLOW_VAR
      USE TWO_PHASE_PROPERTY
      
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR

      USE HEAT_VAR
      
      USE IBM_VAR

        REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
        REAL P(0:M1,0:M2,0:M3)
        real alphi_save(-2:m1f+3,-2:m2f+3,-2:m3f+3,mlvs)
       INTEGER LVS
      COMMON/GLOBAL_MASS_CORRECTION1/VOL_TOT_ORI(MLVS)

      I_B=0.
      J_B=0.
      Z_B=0.
      ALPHI_EXT=0.

!     dum for future use
      READ(12,1001) N1,N2,N3
      READ(12,1002) IHIST,M,TIME,DT
      READ(12,1003) IPS,IPX,IPY,IPZ,PRAA
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
        if (n .lt. 10) alphi_save(i,j,k,llvs)=alphi_tmp
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
      READ(12,1006) PRM,SCR,TCR
      READ(12,1005) ((( T(I,J,K) ,I=0,N1),J=0,N2),K=0,N3)
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

C============================FOR TWO-PHASE=============================C
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
C============================FOR TWO-PHASE=============================C

       WRITE(*,*) 'TIME = ',TIME

!BOUNDARY_CONDITION
        IF (IPX .EQ. 0) THEN
          DO K=1,N3M
          DO J=1,N2M
          U(1,J,K)=U(2,J,K)
          V(0,J,K)=V(1,J,K)
          W(0,J,K)=W(1,J,K)
          U(N1,J,K)=U(N1M,J,K)
          V(N1,J,K)=V(N1M,J,K)
          W(N1,J,K)=W(N1M,J,K)
         ENDDO
         ENDDO
        ENDIF
        IF (IPY .EQ. 0) THEN
          DO K=1,N3M
          DO I=1,N1M
          U(I,0,K)=U(I,1,K)
          V(I,1,K)=V(I,2,K)
          W(I,0,K)=W(I,1,K)
          U(I,N2,K)=U(I,N2M,K)
          V(I,N2,K)=V(I,N2M,K)
          W(I,N2,K)=W(I,N2M,K) 
         ENDDO
         ENDDO
        ENDIF
        IF (IPZ .EQ. 0) THEN
          DO J=1,N2M
          DO I=1,N1M
          U(I,J,0)=U(I,J,1)
          V(I,J,0)=V(I,J,1)
          W(I,J,1)=W(I,J,2)
          U(I,J,N3)=U(I,J,N3M)
          V(I,J,N3)=V(I,J,N3M)
          W(I,J,N3)=W(I,J,N3M) 
         ENDDO
         ENDDO
        ENDIF


!=====SAVE
       llvs=1
       OPEN(146,FILE='0alphi_savepost.DAT')
       WRITE(146,*) 'VARIABLES="X","Y","Z","alphi"'
      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
      DO K=1,N3FM
      DO J=1,N2FM
      DO I=1,N1FM
        WRITE(146,152) XPF(I),YPF(J),ZPF(K),alphi_save(i,j,k,llvs)
      ENDDO
      ENDDO
      ENDDO
       CLOSE(146)
 152  FORMAT(4F15.8) 


        
      RETURN
      END

C==============================================================
C    COMPUTE THE VELOCITIES OF CELL-CENTER
      SUBROUTINE CENTERVEL(U,V,W,P)
C==============================================================
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE HEAT_VAR
        
        REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
        REAL P(0:M1,0:M2,0:M3)
        
        COMMON/CENTER/VPW(0:M1,0:M2,0:M3,0:9)

        DO K=1,N3M
         KP1=KPV(K)
        DO J=1,N2M
         JP1=JPV(J)
        DO I=1,N1M
         IP1=IPV(I)
         VPW(I,J,K,1)=0.5*(U(I,J,K)+U(IP1,J,K))
         VPW(I,J,K,2)=0.5*(V(I,J,K)+V(I,JP1,K))
         VPW(I,J,K,3)=0.5*(W(I,J,K)+W(I,J,KP1))
         VPW(I,J,K,4)=P(I,J,K)
         VPW(I,J,K,0)=T(I,J,K)
        END DO
        END DO
        END DO
        

      
      RETURN
      END

C=======================================================================
      SUBROUTINE P_FLUC(P,PSI_CN)
C=======================================================================
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE TWO_PHASE_PROPERTY

        REAL P(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)

        IF (DENR .EQ. 1) THEN
         IPSI=0	
        ELSE
         IPSI=1
        ENDIF

        P_AVG=0.
        VOL=0.
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
        	XX=XP(I)
        	YY=YP(J)
        	ZZ=ZP(K)
         IF (FUNCBODY(XX,YY,ZZ) .GT. 0. .AND.
     &       FUNCBODY(XP(IPV(I)),YY,ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0. !NOT NEAR THE INTERFACE
     & .AND. FUNCBODY(XP(IMV(I)),YY,ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0.
     & .AND. FUNCBODY(XX,YP(JPV(J)),ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0.
     & .AND. FUNCBODY(XX,YP(JMV(J)),ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0.) THEN
         	IF (PSI_CN(I,J,K) .EQ. IPSI) THEN
         	 DVOL=SDX(I)*SDY(J)*SDZ(K)
         	 P_AVG=P_AVG+P(I,J,K)*DVOL
         	 VOL=VOL+DVOL
          ENDIF
         ENDIF	
        ENDDO
        ENDDO
        ENDDO
        
        IF (VOL .NE. 0.) THEN
         P_AVG=P_AVG/VOL
        ELSE
         WRITE(*,*) 'VOL=0. FOR CAL. P_AVG'
         P_AVG=0.
        ENDIF

        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
        	XX=XP(I)
        	YY=YP(J)
        	ZZ=ZP(K)
         IF (FUNCBODY(XX,YY,ZZ) .GT. 0. .AND.
     &       FUNCBODY(XP(IPV(I)),YY,ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0. !NOT NEAR THE INTERFACE
     & .AND. FUNCBODY(XP(IMV(I)),YY,ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0.
     & .AND. FUNCBODY(XX,YP(JPV(J)),ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0.
     & .AND. FUNCBODY(XX,YP(JMV(J)),ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0.) THEN
          P(I,J,K)=P(I,J,K)-P_AVG
         ELSE
          P(I,J,K)=0.
         ENDIF	
        ENDDO
        ENDDO
        ENDDO

      RETURN
      END
      
C=======================================================================
      SUBROUTINE XCUT(I_POS,LVSGRID,ITIMESTEP,PSI_CN
     &,IVEL_CYLIND,IVOR_CYLIND)
C=======================================================================
      USE FLOW_VAR
      
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING
      
      USE HEAT_VAR

      REAL PSI_CN(0:M1,0:M2,0:M3)
      
        COMMON/CENTER/VPW(0:M1,0:M2,0:M3,0:9)
        COMMON/LAMDA2F/VLAMBDA2(0:M1,0:M2,0:M3)

        CHARACTER*3 tfn1
        CHARACTER*4 tfn2
        CHARACTER*2 tfn3       
        CHARACTER*30 outfile
        
        PI=ACOS(-1.)
        
      tfn1='fld'
      tfn2='.dat'
      tfn3='VW'    
        
      idg1=ITIMESTEP/100000
      idg2=(ITIMESTEP-idg1*100000)/10000
      idg3=(ITIMESTEP-idg1*100000-idg2*10000)/1000
      idg4=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6= ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10

      outfile=tfn3//char(idg1+48)//char(idg2+48)//char(idg3+48)//
     &        char(idg4+48)//char(idg5+48)//char(idg6+48)//tfn2 
     
      I=I_POS
      IF (LVSGRID .EQ. 0) THEN
        OPEN(1000,FILE=outfile,STATUS='UNKNOWN')
        WRITE(1000,*) 'VARIABLES="YY","ZZ","VV","WW","PP","VOLF"'
        WRITE(1000,*) 'ZONE I=',N2M,',J=',N3M,',F=POINT'

        DO K=1,N3M
        DO J=1,N2M
         YY=YP(J)
         ZZ=ZP(K)
         VV=VPW(I,J,K,2)
         WW=VPW(I,J,K,3)
         PP=VPW(I,J,K,4)
         WRITE(1000,10) YY,ZZ,VV,WW,PP,PSI_CN(I,J,K)
        END DO
        END DO

   10  FORMAT(2F12.6,4E14.6)
       CLOSE(1000)
       
      ELSE IF (LVSGRID .EQ. 1) THEN

        OPEN(1000,FILE=outfile,STATUS='UNKNOWN')
       IF (IVOR_CYLIND .EQ. 0) THEN
       	 IF (IVEL_CYLIND .EQ. 0) THEN
        WRITE(1000,*) 'VARIABLES="YY","ZZ","VV","WW","PP","VOLF"
     &,"VORX","VORY","VORZ","LAMDA2"'
         ELSE IF (IVEL_CYLIND .EQ. 1) THEN
        WRITE(1000,*) 'VARIABLES="YY","ZZ","UT","WW","PP","VOLF"
     &,"VORX","VORY","VORZ","LAMDA2"' 
        ENDIF
       ELSE IF (IVOR_CYLIND .EQ. 1) THEN
       	 IF (IVEL_CYLIND .EQ. 0) THEN
        WRITE(1000,*) 'VARIABLES="YY","ZZ","VV","WW","PP","VOLF"
     &,"VORR","VORT","VORZ","LAMDA2"'
         ELSE IF (IVEL_CYLIND .EQ. 1) THEN
        WRITE(1000,*) 'VARIABLES="YY","ZZ","UT","WW","PP","VOLF"
     &,"VORR","VORT","VORZ","LAMDA2"' 
        ENDIF
        ENDIF
        WRITE(1000,*) 'ZONE I=',N2M,',J=',N3M,',F=POINT'

        DO K=1,N3M
        DO J=1,N2M
         YY=YP(J)
         ZZ=ZP(K)
         VV=VPW(I,J,K,2)
         WW=VPW(I,J,K,3)
         PP=VPW(I,J,K,4)
         VOR1=VPW(I,J,K,5)
         VOR2=VPW(I,J,K,6)
         VOR3=VPW(I,J,K,7)

         XX=X(I)
         TT_TMP=0.
         IF (IVEL_CYLIND .EQ. 1) THEN
         TT_TMP=ATAN(YY/XX)
         IF (XX .LE. 0. .AND. YY .GE. 0.) THEN !2ND QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .LE. 0. .AND. YY .LE. 0.) THEN !3RD QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .GE. 0. .AND. YY .LE. 0.) THEN !4TH QUADRANT
           TT_TMP=TT_TMP+2.*PI
         ENDIF

         UR=UU*COS(TT_TMP)+VV*SIN(TT_TMP)
         UT=-UU*SIN(TT_TMP)+VV*COS(TT_TMP)

         UU=UR
         VV=UT
         ENDIF   

         IF (IVOR_CYLIND .EQ. 1) THEN
          IF (TT_TMP .EQ. 0.) THEN
         TT_TMP=ATAN(YY/XX)
         IF (XX .LE. 0. .AND. YY .GE. 0.) THEN !2ND QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .LE. 0. .AND. YY .LE. 0.) THEN !3RD QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .GE. 0. .AND. YY .LE. 0.) THEN !4TH QUADRANT
           TT_TMP=TT_TMP+2.*PI
         ENDIF
         ENDIF

         WR=VOR1*COS(TT_TMP)+VOR2*SIN(TT_TMP)
         WT=-VOR1*SIN(TT_TMP)+VOR2*COS(TT_TMP)

         VOR1=WR
         VOR2=WT
         ENDIF      

         WRITE(1000,20) YY,ZZ,VV,WW,PP,PSI_CN(I,J,K)
     &,VOR1,VOR2,VOR3,VLAMBDA2(I,J,K)
        END DO
        END DO

   20  FORMAT(2F12.6,8E14.6)
       CLOSE(1000)
      	
      ELSE IF (LVSGRID .EQ. 2) THEN
        OPEN(1000,FILE=outfile,STATUS='UNKNOWN')
        WRITE(1000,*) 'VARIABLES="YY","ZZ","ALPHI"'
        WRITE(1000,*) 'ZONE I=',N3FM,',J=',N2FM,',F=POINT'

        I=ICOU2(I_POS)
        DO K=1,N3FM
        DO J=1,N2FM
         YY=YPF(J)
         ZZ=ZPF(K)
         
       	 ALPHI_MIN=1000000.
       	 ALPHI_TMP=0.
         DO LLVS=1,NLVS
       	  IF (ALPHI_EXT(I,J,K,LLVS) .NE. 0.) THEN
       	   ALPHI_MIN=AMIN1(ABS(ALPHI_EXT(I,J,K,LLVS)),ALPHI_MIN)
            IF (ALPHI_MIN .EQ. ABS(ALPHI_EXT(I,J,K,LLVS))) 
     &                              ALPHI_TMP=ALPHI_EXT(I,J,K,LLVS)
          ENDIF
         ENDDO

         WRITE(1000,30) YY,ZZ,ALPHI_TMP
        END DO
        END DO

   30  FORMAT(2F12.6,1000E12.4)
       CLOSE(1000)

      ENDIF

       RETURN
       END

C==================================================================
      SUBROUTINE YCUT(J_POS,LVSGRID,ITIMESTEP,PSI_CN
     &,IVEL_CYLIND,IVOR_CYLIND)
C==================================================================
      USE FLOW_VAR
      
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING
      
      USE HEAT_VAR

      REAL PSI_CN(0:M1,0:M2,0:M3)
      
        COMMON/CENTER/VPW(0:M1,0:M2,0:M3,0:9)
        COMMON/LAMDA2F/VLAMBDA2(0:M1,0:M2,0:M3)

        CHARACTER*3 tfn1
        CHARACTER*4 tfn2
        CHARACTER*2 tfn3       
        CHARACTER*30 outfile
        
        PI=ACOS(-1.)
        
      tfn1='fld'
      tfn2='.dat'
      tfn3='UW'   

      idg1=ITIMESTEP/100000
      idg2=(ITIMESTEP-idg1*100000)/10000
      idg3=(ITIMESTEP-idg1*100000-idg2*10000)/1000
      idg4=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6= ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10

      outfile=tfn3//char(idg1+48)//char(idg2+48)//char(idg3+48)//
     &        char(idg4+48)//char(idg5+48)//char(idg6+48)//tfn2 

      J=J_POS
      IF (LVSGRID .EQ. 0) THEN
        OPEN(1000,FILE=outfile,STATUS='UNKNOWN')
        WRITE(1000,*) 'VARIABLES="XX","ZZ","UU","WW","PP","VOLF"'
        WRITE(1000,*) 'ZONE I=',N1M,',J=',N3M,',F=POINT'

        DO K=1,N3M
        DO I=1,N1M
         XX=XP(I)
         ZZ=ZP(K)
         UU=VPW(I,J,K,1)
         WW=VPW(I,J,K,3)
         PP=VPW(I,J,K,4)
         WRITE(1000,10) XX,ZZ,UU,WW,PP,PSI_CN(I,J,K)
        END DO
        END DO

  10   FORMAT(2F12.6,4E14.6)
       CLOSE(1000)
       
      ELSE IF (LVSGRID .EQ. 1) THEN
        OPEN(1000,FILE=outfile,STATUS='UNKNOWN')
       IF (IVOR_CYLIND .EQ. 0) THEN
       	 IF (IVEL_CYLIND .EQ. 0) THEN
        WRITE(1000,*) 'VARIABLES="XX","ZZ","UU","WW","PP","VOLF"
     &,"VORX","VORY","VORZ","LAMDA2"'
         ELSE IF (IVEL_CYLIND .EQ. 1) THEN
        WRITE(1000,*) 'VARIABLES="XX","ZZ","UR","WW","PP","VOLF"
     &,"VORX","VORY","VORZ","LAMDA2"'
        ENDIF
       ELSE IF (IVOR_CYLIND .EQ. 1) THEN
       	 IF (IVEL_CYLIND .EQ. 0) THEN
        WRITE(1000,*) 'VARIABLES="XX","ZZ","UU","WW","PP","VOLF"
     &,"VORR","VORT","VORZ","LAMDA2"'
         ELSE IF (IVEL_CYLIND .EQ. 1) THEN
        WRITE(1000,*) 'VARIABLES="XX","ZZ","UR","WW","PP","VOLF"
     &,"VORR","VORT","VORZ","LAMDA2"'
        ENDIF
        ENDIF
        WRITE(1000,*) 'ZONE I=',N1M,',J=',N3M,',F=POINT'

        DO K=1,N3M
        DO I=1,N1M
         XX=XP(I)
         ZZ=ZP(K)
         UU=VPW(I,J,K,1)
         WW=VPW(I,J,K,3)
         PP=VPW(I,J,K,4)
         VOR1=VPW(I,J,K,5)
         VOR2=VPW(I,J,K,6)
         VOR3=VPW(I,J,K,7)
         
         YY=Y(J)
         TT_TMP=0.
         IF (IVEL_CYLIND .EQ. 1) THEN
         TT_TMP=ATAN(YY/XX)
         IF (XX .LE. 0. .AND. YY .GE. 0.) THEN !2ND QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .LE. 0. .AND. YY .LE. 0.) THEN !3RD QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .GE. 0. .AND. YY .LE. 0.) THEN !4TH QUADRANT
           TT_TMP=TT_TMP+2.*PI
         ENDIF

         UR=UU*COS(TT_TMP)+VV*SIN(TT_TMP)
         UT=-UU*SIN(TT_TMP)+VV*COS(TT_TMP)

         UU=UR
         VV=UT
         ENDIF   

         IF (IVOR_CYLIND .EQ. 1) THEN
          IF (TT_TMP .EQ. 0.) THEN
         TT_TMP=ATAN(YY/XX)
         IF (XX .LE. 0. .AND. YY .GE. 0.) THEN !2ND QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .LE. 0. .AND. YY .LE. 0.) THEN !3RD QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .GE. 0. .AND. YY .LE. 0.) THEN !4TH QUADRANT
           TT_TMP=TT_TMP+2.*PI
         ENDIF
         ENDIF

         WR=VOR1*COS(TT_TMP)+VOR2*SIN(TT_TMP)
         WT=-VOR1*SIN(TT_TMP)+VOR2*COS(TT_TMP)

         VOR1=WR
         VOR2=WT
         ENDIF      

         WRITE(1000,20) XX,ZZ,UU,WW,PP,PSI_CN(I,J,K)
     &,VOR1,VOR2,VOR3,VLAMBDA2(I,J,K)
        END DO
        END DO

  20   FORMAT(2F12.6,8E14.6)
       CLOSE(1000)
      	
      ELSE IF (LVSGRID .EQ. 2) THEN
        OPEN(1000,FILE=outfile,STATUS='UNKNOWN')
        WRITE(1000,*) 'VARIABLES="XX","ZZ","UU","ALPHI"'
        WRITE(1000,*) 'ZONE I=',N1FM,',J=',N3FM,',F=POINT'

        J=JCOU2(J_POS)
        DO K=1,N3FM
        DO I=1,N1FM
         XX=XPF(I)
         ZZ=ZPF(K)
         
       	 ALPHI_MIN=1000000.
       	 ALPHI_TMP=0.
         DO LLVS=1,NLVS
       	  IF (ALPHI_EXT(I,J,K,LLVS) .NE. 0.) THEN
       	   ALPHI_MIN=AMIN1(ABS(ALPHI_EXT(I,J,K,LLVS)),ALPHI_MIN)
            IF (ALPHI_MIN .EQ. ABS(ALPHI_EXT(I,J,K,LLVS))) 
     &                              ALPHI_TMP=ALPHI_EXT(I,J,K,LLVS)
          ENDIF
         ENDDO
         
         WRITE(1000,30) XX,ZZ,ALPHI_TMP
        END DO
        END DO

  30   FORMAT(2F12.6,1000E12.4)
       CLOSE(1000)
      	
      ENDIF



       RETURN
       END

C===================================================================
      SUBROUTINE ZCUT(K_POS,LVSGRID,ITIMESTEP,PSI_CN
     &,IVEL_CYLIND,IVOR_CYLIND)
C===================================================================
       USE FLOW_VAR
      
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING
      
      USE HEAT_VAR

      REAL PSI_CN(0:M1,0:M2,0:M3)
      
        COMMON/CENTER/VPW(0:M1,0:M2,0:M3,0:9)
        COMMON/LAMDA2F/VLAMBDA2(0:M1,0:M2,0:M3)

        CHARACTER*3 tfn1
        CHARACTER*4 tfn2
        CHARACTER*2 tfn3       
        CHARACTER*30 outfile
        
        PI=ACOS(-1.)
        
      tfn1='fld'
      tfn2='.dat'
      tfn3='UV'   

      idg1=ITIMESTEP/100000
      idg2=(ITIMESTEP-idg1*100000)/10000
      idg3=(ITIMESTEP-idg1*100000-idg2*10000)/1000
      idg4=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6= ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10

      outfile=tfn3//char(idg1+48)//char(idg2+48)//char(idg3+48)//
     &        char(idg4+48)//char(idg5+48)//char(idg6+48)//tfn2 

      K=K_POS
      IF (LVSGRID .EQ. 0) THEN
        OPEN(1000,FILE=outfile,STATUS='UNKNOWN')
        WRITE(1000,*) 'VARIABLES="XX","YY","UU","VV","PP","VOLF"'
        WRITE(1000,*) 'ZONE I=',N1M,',J=',N2M,',F=POINT'

        DO J=1,N2M
        DO I=1,N1M
         XX=XP(I)
         YY=YP(J)
         UU=VPW(I,J,K,1)
         VV=VPW(I,J,K,2)
         PP=VPW(I,J,K,4)
        WRITE(1000,10) XX,YY,UU,VV,PP,PSI_CN(I,J,K)
        END DO
        END DO

  10   FORMAT(2F12.6,4E14.6)
       CLOSE(1000)

      ELSE IF (LVSGRID .EQ. 1) THEN
        OPEN(1000,FILE=outfile,STATUS='UNKNOWN')
        OPEN(1000,FILE=outfile,STATUS='UNKNOWN')
       IF (IVOR_CYLIND .EQ. 0) THEN
       	 IF (IVEL_CYLIND .EQ. 0) THEN
        WRITE(1000,*) 'VARIABLES="XX","YY","UU","VV","PP","VOLF"
     &,"VORX","VORY","VORZ","LAMDA2"'
         ELSE IF (IVEL_CYLIND .EQ. 1) THEN
        WRITE(1000,*) 'VARIABLES="XX","YY","UR","UT","PP","VOLF"
     &,"VORX","VORY","VORZ","LAMDA2"'
        ENDIF
       ELSE IF (IVOR_CYLIND .EQ. 1) THEN
       	 IF (IVEL_CYLIND .EQ. 0) THEN
        WRITE(1000,*) 'VARIABLES="XX","YY","UU","VV","PP","VOLF"
     &,"VORR","VORT","VORZ","LAMDA2"'
         ELSE IF (IVEL_CYLIND .EQ. 1) THEN
        WRITE(1000,*) 'VARIABLES="XX","YY","UR","UT","PP","VOLF"
     &,"VORR","VORT","VORZ","LAMDA2"'
        ENDIF
        ENDIF
        WRITE(1000,*) 'ZONE I=',N1M,',J=',N2M,',F=POINT'

        DO J=1,N2M
        DO I=1,N1M
         XX=XP(I)
         YY=YP(J)
         UU=VPW(I,J,K,1)
         VV=VPW(I,J,K,2)
         PP=VPW(I,J,K,4)
         VOR1=VPW(I,J,K,5)
         VOR2=VPW(I,J,K,6)
         VOR3=VPW(I,J,K,7)
         
         TT_TMP=0.
         IF (IVEL_CYLIND .EQ. 1) THEN
         TT_TMP=ATAN(YY/XX)
         IF (XX .LE. 0. .AND. YY .GE. 0.) THEN !2ND QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .LE. 0. .AND. YY .LE. 0.) THEN !3RD QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .GE. 0. .AND. YY .LE. 0.) THEN !4TH QUADRANT
           TT_TMP=TT_TMP+2.*PI
         ENDIF

         UR=UU*COS(TT_TMP)+VV*SIN(TT_TMP)
         UT=-UU*SIN(TT_TMP)+VV*COS(TT_TMP)

         UU=UR
         VV=UT
         ENDIF   

         IF (IVOR_CYLIND .EQ. 1) THEN
          IF (TT_TMP .EQ. 0.) THEN
         TT_TMP=ATAN(YY/XX)
         IF (XX .LE. 0. .AND. YY .GE. 0.) THEN !2ND QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .LE. 0. .AND. YY .LE. 0.) THEN !3RD QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .GE. 0. .AND. YY .LE. 0.) THEN !4TH QUADRANT
           TT_TMP=TT_TMP+2.*PI
         ENDIF
         ENDIF

         WR=VOR1*COS(TT_TMP)+VOR2*SIN(TT_TMP)
         WT=-VOR1*SIN(TT_TMP)+VOR2*COS(TT_TMP)

         VOR1=WR
         VOR2=WT
         ENDIF      

        WRITE(1000,20) XX,YY,UU,VV,PP,PSI_CN(I,J,K)
     &,VOR1,VOR2,VOR3,VLAMBDA2(I,J,K)
        END DO
        END DO

  20   FORMAT(2F12.6,8E14.6)
       CLOSE(1000)
      	
      ELSE IF (LVSGRID .EQ. 2) THEN
        OPEN(1000,FILE=outfile,STATUS='UNKNOWN')
        WRITE(1000,*) 'VARIABLES="XX","YY","ALPHI"'
        WRITE(1000,*) 'ZONE I=',N1FM,',J=',N2FM,',F=POINT'

        K=KCOU2(K_POS)
        DO J=1,N2FM
        DO I=1,N1FM
         XX=XPF(I)
         YY=YPF(J)
         
       	 ALPHI_MIN=1000000.
       	 ALPHI_TMP=0.
         DO LLVS=1,NLVS
       	  IF (ALPHI_EXT(I,J,K,LLVS) .NE. 0.) THEN
       	   ALPHI_MIN=AMIN1(ABS(ALPHI_EXT(I,J,K,LLVS)),ALPHI_MIN)
            IF (ALPHI_MIN .EQ. ABS(ALPHI_EXT(I,J,K,LLVS))) 
     &                              ALPHI_TMP=ALPHI_EXT(I,J,K,LLVS)
          ENDIF
         ENDDO
         
        WRITE(1000,30) XX,YY,ALPHI_TMP
        END DO
        END DO

  30   FORMAT(2F12.6,1000E12.4)
       CLOSE(1000)
      ENDIF
       

       RETURN
       END
C*******************************************************************
      SUBROUTINE THREE_DIMENSION(LVSGRID,ITIMESTEP,PSI_CN,U,V,W
     &,IVEL_CYLIND,IVOR_CYLIND)
C*******************************************************************
C       get 3D data
      USE FLOW_VAR
      
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR

      USE HEAT_VAR
      
      USE TWO_PHASE_PROPERTY

      REAL PSI_CN(0:M1,0:M2,0:M3)

      REAL, ALLOCATABLE ::CONV(:,:,:,:),CONV_FILTER(:,:,:,:)
      REAL, ALLOCATABLE ::ST(:,:,:,:),ST_FILTER(:,:,:,:)
      REAL, ALLOCATABLE ::DIF(:,:,:,:),DIF_FILTER(:,:,:,:)
      REAL, ALLOCATABLE::U_FILTER(:,:,:),V_FILTER(:,:,:),W_FILTER(:,:,:)
      REAL, ALLOCATABLE ::PSI_FILTER(:,:,:)
      REAL, ALLOCATABLE ::TAU2(:,:,:,:),TAU3(:,:,:,:)
      
       COMMON/CENTER/VPW(0:M1,0:M2,0:M3,0:9)
       COMMON/LAMDA2F/VLAMBDA2(0:M1,0:M2,0:M3)

       COMMON/SKIP1/ISKIP,JSKIP,KSKIP
       COMMON/SKIP2/XBIG,YBIG,ZBIG,XEND,YEND,ZEND

        CHARACTER*3 tfn1
        CHARACTER*4 tfn2
        CHARACTER*2 tfn3       
        CHARACTER*30 outfile

       REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)

        PI=ACOS(-1.)
        
      IF ( XBIG .EQ. 0. .AND. XEND .EQ. 0.) THEN
      	IBIG=1
      	IEND=N1M
      ELSE
       DO I=1,N1M
        IF (XP(I+1) .GE. XBIG) THEN
          IBIG=I
          GOTO 10
        ENDIF
       ENDDO
  10   CONTINUE
       DO I=1,N1M
        IF (XP(I+1) .GE. XEND) THEN
          IEND=I
          GOTO 11
        ENDIF
       ENDDO
  11   CONTINUE
      ENDIF
  
      IF ( YBIG .EQ. 0. .AND. YEND .EQ. 0.) THEN
      	JBIG=1
      	JEND=N2M
      ELSE
       DO J=1,N2M
         IF (YP(J+1) .GE. YBIG) THEN
            JBIG=J
            GOTO 12
         ENDIF
       ENDDO
 12   CONTINUE
       DO J=1,N2M
         IF (YP(J+1) .GE. YEND) THEN
            JEND=J
            GOTO 13
         ENDIF
       ENDDO
 13   CONTINUE
      ENDIF
      
      IF ( ZBIG .EQ. 0. .AND. ZEND .EQ. 0.) THEN
      	KBIG=1
      	KEND=N3M
      ELSE
       DO K=1,N3M
         IF (ZP(K+1) .GE. ZBIG) THEN
            KBIG=K
            GOTO 14
         ENDIF
       ENDDO
 14   CONTINUE
       DO K=1,N3M
         IF (ZP(K+1) .GE. ZEND) THEN
            KEND=K
            GOTO 15
         ENDIF
       ENDDO
 15   CONTINUE
      ENDIF

       N1MS=(IEND-IBIG)/ISKIP+1
       N2MS=(JEND-JBIG)/JSKIP+1
       N3MS=(KEND-KBIG)/KSKIP+1
       
      tfn1='fld'
      tfn2='.dat'
      tfn3='3d'    
        
      idg1=ITIMESTEP/100000
      idg2=(ITIMESTEP-idg1*100000)/10000
      idg3=(ITIMESTEP-idg1*100000-idg2*10000)/1000
      idg4=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6= ITIMESTEP-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10

      outfile=tfn3//char(idg1+48)//char(idg2+48)//char(idg3+48)//
     &        char(idg4+48)//char(idg5+48)//char(idg6+48)//tfn2 

       IF ( LVSGRID .EQ. 0) THEN

        OPEN(42,FILE=outfile)
        IF (IPS .EQ. 0) THEN
       WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","U","V","W","P","VOLF"
     &,"BODY"'   
        ELSE IF (IPS .EQ. 1) THEN
       WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","U","V","W","P","VOLF"
     &,"BODY","TEMP"'
        ENDIF
              
        WRITE(42,*) 'ZONE I=',N1MS,',J=',N2MS,',K=',N3MS,',F=POINT'
        DO K=KBIG,KEND,KSKIP
        DO J=JBIG,JEND,JSKIP
        DO I=IBIG,IEND,ISKIP
         XX=XP(I)
         YY=YP(J)
         ZZ=ZP(K)
         UU=VPW(I,J,K,1)
         VV=VPW(I,J,K,2)
         WW=VPW(I,J,K,3)
         PP=VPW(I,J,K,4)
         VOLF=PSI_CN(I,J,K)
         BODY=FUNCBODY(XX,YY,ZZ)
         
        IF (IPS .EQ. 0) THEN
         WRITE(42,43)  XX,YY,ZZ,UU,VV,WW,PP,VOLF,BODY
        ELSE IF (IPS .EQ. 1) THEN
         TEMP=T(I,J,K)
         WRITE(42,43)  XX,YY,ZZ,UU,VV,WW,PP,VOLF,BODY,TEMP
        ENDIF

       ENDDO
       ENDDO
       ENDDO
        CLOSE(42)
 43     FORMAT(3F10.5,7E12.4)

       ELSE IF ( LVSGRID .EQ. 1 ) THEN

       CALL VORNLAMBDA2

        OPEN(42,FILE=outfile)
       IF (IVOR_CYLIND .EQ. 0) THEN
       	 IF (IVEL_CYLIND .EQ. 0) THEN
          IF (IPS .EQ. 0) THEN
        WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","U","V","W","P"
     &,"VOLF","BODY","wx","wy","wz","L2"' 
          ELSE IF (IPS .EQ. 1) THEN
        WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","U","V","W","P"
     &,"VOLF","BODY","wx","wy","wz","L2","TEMP"'
          ENDIF
         ELSE IF (IVEL_CYLIND .EQ. 1) THEN
          IF (IPS .EQ. 0) THEN
        WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","UR","UT","W","P"
     &,"VOLF","BODY","wx","wy","wz","L2"'    
          ELSE IF (IPS .EQ. 1) THEN
        WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","UR","UT","W","P"
     &,"VOLF","BODY","wx","wy","wz","L2","TEMP"'    
          ENDIF
        ENDIF
       ELSE IF (IVOR_CYLIND .EQ. 1) THEN
       	 IF (IVEL_CYLIND .EQ. 0) THEN
          IF (IPS .EQ. 0) THEN
        WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","U","V","W","P"
     &,"VOLF","BODY","wR","wT","wz","L2"'   
          ELSE IF (IPS .EQ. 1) THEN
        WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","U","V","W","P"
     &,"VOLF","BODY","wR","wT","wz","L2","TEMP"'   
          ENDIF
         ELSE IF (IVEL_CYLIND .EQ. 1) THEN
          IF (IPS .EQ. 0) THEN
        WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","UR","UT","W","P"
     &,"VOLF","BODY","wR","wT","wz","L2"'  
          ELSE IF (IPS .EQ. 1) THEN
        WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","UR","UT","W","P"
     &,"VOLF","BODY","wR","wT","wz","L2","TEMP"'   
          ENDIF

        ENDIF
        ENDIF
       
        WRITE(42,*) 'ZONE I=',N1MS,',J=',N2MS,',K=',N3MS,',F=POINT'
        DO K=KBIG,KEND,KSKIP
        DO J=JBIG,JEND,JSKIP
        DO I=IBIG,IEND,ISKIP
         XX=XP(I)
         YY=YP(J)
         ZZ=ZP(K)
         UU=VPW(I,J,K,1)
         VV=VPW(I,J,K,2)
         WW=VPW(I,J,K,3)
         PP=VPW(I,J,K,4)
         VOLF=PSI_CN(I,J,K)
         BODY=FUNCBODY(XX,YY,ZZ)
         VOR1=VPW(I,J,K,5)
         VOR2=VPW(I,J,K,6)
         VOR3=VPW(I,J,K,7)

         TT_TMP=0.
         IF (IVEL_CYLIND .EQ. 1) THEN
         TT_TMP=ATAN(YY/XX)
         IF (XX .LE. 0. .AND. YY .GE. 0.) THEN !2ND QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .LE. 0. .AND. YY .LE. 0.) THEN !3RD QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .GE. 0. .AND. YY .LE. 0.) THEN !4TH QUADRANT
           TT_TMP=TT_TMP+2.*PI
         ENDIF

         UR=UU*COS(TT_TMP)+VV*SIN(TT_TMP)
         UT=-UU*SIN(TT_TMP)+VV*COS(TT_TMP)

         UU=UR
         VV=UT
         ENDIF   

         IF (IVOR_CYLIND .EQ. 1) THEN
          IF (TT_TMP .EQ. 0.) THEN
         TT_TMP=ATAN(YY/XX)
         IF (XX .LE. 0. .AND. YY .GE. 0.) THEN !2ND QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .LE. 0. .AND. YY .LE. 0.) THEN !3RD QUADRANT
           TT_TMP=TT_TMP+PI
         ELSE IF (XX .GE. 0. .AND. YY .LE. 0.) THEN !4TH QUADRANT
           TT_TMP=TT_TMP+2.*PI
         ENDIF
         ENDIF

         WR=VOR1*COS(TT_TMP)+VOR2*SIN(TT_TMP)
         WT=-VOR1*SIN(TT_TMP)+VOR2*COS(TT_TMP)

         VOR1=WR
         VOR2=WT
         ENDIF         

         AL2=VLAMBDA2(I,J,K)
        
        IF (IPS .EQ. 0) THEN
         WRITE(42,44) XX,YY,ZZ,UU,VV,WW,PP,VOLF,BODY,VOR1,VOR2,VOR3,AL2
        ELSE IF (IPS .EQ. 1) THEN
         TEMP=T(I,J,K)
         WRITE(42,44) XX,YY,ZZ,UU,VV,WW,PP,VOLF,BODY,VOR1,VOR2,VOR3,AL2
     &  ,TEMP
        ENDIF

       ENDDO
       ENDDO
       ENDDO
        CLOSE(42)
 44     FORMAT(3F10.5,11E12.4)

       ELSE IF ( LVSGRID .EQ. 2 ) THEN
        OPEN(42,FILE=outfile)
        WRITE(42,*) 'VARIABLES="X","Y","Z","ALPHI_MIN"'
        WRITE(42,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
        DO K=1,N3FM
        DO J=1,N2FM
        DO I=1,N1FM
        	XX=XPF(I)
        	YY=YPF(J)
        	ZZ=ZPF(K)

       	 ALPHI_MIN=1000000.
       	 ALPHI_TMP=0.
         DO LLVS=1,NLVS
       	  IF (ALPHI_EXT(I,J,K,LLVS) .NE. 0.) THEN
       	   ALPHI_MIN=AMIN1(ABS(ALPHI_EXT(I,J,K,LLVS)),ALPHI_MIN)
            IF (ALPHI_MIN .EQ. ABS(ALPHI_EXT(I,J,K,LLVS))) 
     &                              ALPHI_TMP=ALPHI_EXT(I,J,K,LLVS)
          ENDIF
         ENDDO
          WRITE(42,45) XX,YY,ZZ,ALPHI_TMP
       ENDDO
       ENDDO
       ENDDO
        CLOSE(42)
 45     FORMAT(3F10.5,1000E12.4)

       ELSE IF (LVSGRID .EQ. 3) THEN

       	WRITE(*,*) 'LVSGRID = 3'
       CALL VORNLAMBDA2

!        IF (DENR .EQ. 1) THEN
!         IPSI=0	
!        ELSE
!         IPSI=1
!        ENDIF

       OPEN(42,FILE=outfile)
       WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","W","VOLF"
     & ,"WX","WY","WZ","lamda2","WX*DW/DX","WY*DW/DY","WZ*DW/DZ",
     & "dwdx","dwdy","VOR_DF1","VOR_DF2","VOR_DF3"'

        WRITE(42,*) 'ZONE I=',N1MS,',J=',N2MS,',K=',N3MS,',F=POINT'
        DO K=KBIG,KEND,KSKIP
        DO J=JBIG,JEND,JSKIP
        DO I=IBIG,IEND,ISKIP
         XX=XP(I)
         YY=YP(J)
         ZZ=ZP(K)
         UU=VPW(I,J,K,1)
         VV=VPW(I,J,K,2)
         WW=VPW(I,J,K,3)
         PP=VPW(I,J,K,4)
         VOLF=PSI_CN(I,J,K)
         BODY=FUNCBODY(XX,YY,ZZ)
         VOR1=VPW(I,J,K,5)
         VOR2=VPW(I,J,K,6)
         VOR3=VPW(I,J,K,7)
         AL2=VLAMBDA2(I,J,K)

       IF (IVEL_CYLIND .EQ. 0 .AND. IVOR_CYLIND .EQ. 0) THEN

         !HELICITY
!         HELICITY=(UU*VOR1+VV*VOR2+WW*VOR3)
!     &/(SQRT(UU**2+VV**2+WW**2)*SQRT(VOR1**2+VOR2**2+VOR3**2))

         !VORTEX_STRECHING TERM
       IF (FUNCBODY(XX,YY,ZZ) .GT. 0. .AND.
     &       FUNCBODY(XP(IPV(I)),YY,ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0. !NOT NEAR THE INTERFACE
     & .AND. FUNCBODY(XP(IMV(I)),YY,ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0.
     & .AND. FUNCBODY(XX,YP(JPV(J)),ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0.
     & .AND. FUNCBODY(XX,YP(JMV(J)),ZZ)*FUNCBODY(XX,YY,ZZ) .GT. 0.) THEN
     	
!!       IF (PSI_CN(I,J,K) .EQ. IPSI) THEN
!       	 VOR_CV1=-UU*(0.5*(VPW(IPV(I),J,K,7)+VPW(I,J,K,7))
!     &  	          -0.5*(VPW(I,J,K,7)+VPW(IMV(I),J,K,7)))*SSDX(I)
!       	 VOR_CV2=-VV*(0.5*(VPW(I,JPV(J),K,7)+VPW(I,J,K,7))
!     &  	          -0.5*(VPW(I,J,K,7)+VPW(I,JMV(J),K,7)))*SSDY(J)
!       	 VOR_CV3=-WW*(0.5*(VPW(I,J,KPV(K),7)+VPW(I,J,K,7))
!     &  	          -0.5*(VPW(I,J,K,7)+VPW(I,J,KMV(K),7)))*SSDZ(K)

         VOR_ST1=VOR1*(VPW(IPV(I),J,K,3)-VPW(IMV(I),J,K,3))
     & /(VDX(I)+VDX(IPV(I)))
         VOR_ST2=VOR2*(VPW(I,JPV(J),K,3)-VPW(I,JMV(J),K,3))
     & /(VDY(J)+VDY(JPV(J)))
         VOR_ST3=VOR3*(W(I,J,KPV(K))-W(I,J,K))*SSDZ(K)
         
       IF (PSI_CN(I,J,K) .EQ. 1) THEN

        VOR_DF1=VISR/DENR*((VPW(IPV(I),J,K,5)-VPW(I,J,K,5))*VVDX(IPV(I))
     &       	 -(VPW(I,J,K,5)-VPW(IMV(I),J,K,5))*VVDX(I))*SSDX(I)
        VOR_DF2=VISR/DENR*((VPW(I,JPV(J),K,6)-VPW(I,J,K,6))*VVDY(JPV(J))
     &      	 -(VPW(I,J,K,6)-VPW(I,JMV(J),K,6))*VVDY(J))*SSDY(J)
        VOR_DF3=VISR/DENR*((VPW(I,J,KPV(K),7)-VPW(I,J,K,7))*VVDZ(KPV(K))
     &       	 -(VPW(I,J,K,7)-VPW(I,J,KMV(K),7))*VVDZ(K))*SSDZ(K)
       ELSE
       	VOR_DF1=0.
       	VOR_DF2=0.
       	VOR_DF3=0.
       ENDIF
         
         VOR_TILTING=VOR_ST1+VOR_ST2
       IF (ABS(VOR_TILTING) .LT. 1.E-3 .AND. ABS(VOR_ST3) .LT.1.E-3)THEN
       	ST_OVER_TIL=0.
       ELSE
        ST_OVER_TIL=VOR_ST3/VOR_TILTING
       ENDIF
!         VOR_TILTING=(VPW(IPV(I),J,K,3)-VPW(IMV(I),J,K,3))
!     &/(VDX(I)+VDX(IPV(I)))!dwdx
!         ST_OVER_TIL=(VPW(I,JPV(J),K,3)-VPW(I,JMV(J),K,3))
!     &/(VDY(J)+VDY(JPV(J)))!dwdy

       ELSE
        VOR_CV1=0.
        VOR_CV2=0.
        VOR_CV3=0.
        VOR_ST1=0.
        VOR_ST2=0.
        VOR_ST3=0.
         	VOR_DF1=0.
         	VOR_DF2=0.
         	VOR_DF3=0.
!       ENDIF !IF (PSI_CN(I,J,K) .EQ. IPSI) THEN
       ENDIF ! IF (FUNCBODY(XX,YY,ZZ) .GT. 0. .AND.
       ELSE
        VOR_CV1=0.
        VOR_CV2=0.
        VOR_CV3=0.
        VOR_ST1=0.
        VOR_ST2=0.
        VOR_ST3=0.
       ENDIF

         WRITE(42,46) XX,YY,ZZ,WW,VOLF,VOR1,VOR2,VOR3,AL2
     &  ,VOR_ST1,VOR_ST2,VOR_ST3,VOR_TILTING,ST_OVER_TIL
     &  ,VOR_DF1,VOR_DF2,VOR_DF3

       ENDDO
       ENDDO
       ENDDO
        CLOSE(42)
 46     FORMAT(3F10.5,14E12.4)
 
       ELSE IF (LVSGRID .EQ. 4) THEN !FILTERING
       	
       	ALLOCATE(CONV(0:M1,0:M2,0:M3,6),CONV_FILTER(0:M1,0:M2,0:M3,6))
        ALLOCATE(ST(0:M1,0:M2,0:M3,9),ST_FILTER(0:M1,0:M2,0:M3,9))
        ALLOCATE(DIF(0:M1,0:M2,0:M3,6),DIF_FILTER(0:M1,0:M2,0:M3,6))
        ALLOCATE(U_FILTER(0:M1,0:M2,0:M3),V_FILTER(0:M1,0:M2,0:M3))
        ALLOCATE(W_FILTER(0:M1,0:M2,0:M3),PSI_FILTER(0:M1,0:M2,0:M3))
        ALLOCATE(TAU2(0:M1,0:M2,0:M3,6),TAU3(0:M1,0:M2,0:M3,6))

       	CALL FILTERING(VPW(0,0,0,1),U_FILTER)
       	CALL FILTERING(VPW(0,0,0,2),V_FILTER)
       	CALL FILTERING(VPW(0,0,0,3),W_FILTER)
       	CALL FILTERING(PSI_CN,PSI_FILTER)
       	
      !convection term
       	DO K=1,N3M
       	DO J=1,N2M
       	DO I=1,N1M
       	DEN=DENM+DEN_DIFF*PSI_FILTER(I,J,K)
       	CONV(I,J,K,1)=DEN*VPW(I,J,K,1)*VPW(I,J,K,1)
       	CONV(I,J,K,2)=DEN*VPW(I,J,K,1)*VPW(I,J,K,2)
       	CONV(I,J,K,3)=DEN*VPW(I,J,K,1)*VPW(I,J,K,3)
       	CONV(I,J,K,4)=DEN*VPW(I,J,K,2)*VPW(I,J,K,2)
       	CONV(I,J,K,5)=DEN*VPW(I,J,K,2)*VPW(I,J,K,3)
       	CONV(I,J,K,6)=DEN*VPW(I,J,K,3)*VPW(I,J,K,3)
        ENDDO
        ENDDO
        ENDDO
       
       	CALL FILTERING(CONV(0,0,0,1),CONV_FILTER(0,0,0,1))
       	CALL FILTERING(CONV(0,0,0,2),CONV_FILTER(0,0,0,2))
       	CALL FILTERING(CONV(0,0,0,3),CONV_FILTER(0,0,0,3))
       	CALL FILTERING(CONV(0,0,0,4),CONV_FILTER(0,0,0,4))
       	CALL FILTERING(CONV(0,0,0,5),CONV_FILTER(0,0,0,5))
       	CALL FILTERING(CONV(0,0,0,6),CONV_FILTER(0,0,0,6))

       	DO K=1,N3M
       	DO J=1,N2M
       	DO I=1,N1M
       	DEN=DENM+DEN_DIFF*PSI_FILTER(I,J,K)
        TAU2(I,J,K,1)=CONV_FILTER(I,J,K,1)
     &                          -DEN*U_FILTER(I,J,K)*U_FILTER(I,J,K)
        TAU2(I,J,K,2)=CONV_FILTER(I,J,K,2)
     &                          -DEN*U_FILTER(I,J,K)*V_FILTER(I,J,K)
        TAU2(I,J,K,3)=CONV_FILTER(I,J,K,3)
     &                          -DEN*U_FILTER(I,J,K)*W_FILTER(I,J,K)
        TAU2(I,J,K,4)=CONV_FILTER(I,J,K,4)
     &                          -DEN*V_FILTER(I,J,K)*V_FILTER(I,J,K)
        TAU2(I,J,K,5)=CONV_FILTER(I,J,K,5)
     &                          -DEN*V_FILTER(I,J,K)*W_FILTER(I,J,K)
        TAU2(I,J,K,6)=CONV_FILTER(I,J,K,6)
     &                          -DEN*W_FILTER(I,J,K)*W_FILTER(I,J,K)

        ENDDO
        ENDDO
        ENDDO
      !convection term
      
      !diffusion term
       	DO K=1,N3M
       	DO J=1,N2M
       	DO I=1,N1M
       	VIS=VISM+VIS_DIFF*PSI_FILTER(I,J,K)
       		
         ST(I,J,K,1)=(U(IPV(I),J,K)-U(I,J,K))*SSDX(I)        !DIF11
         ST(I,J,K,2)=0.5*((U(IPV(I),JPV(J),K)+U(I,JPV(J),K)) !DIF12
     &                    -(U(IPV(I),JMV(J),K)+U(I,JMV(J),K)))*SSDY(J)
         ST(I,J,K,3)=0.5*((U(IPV(I),J,KPV(K))+U(I,J,KPV(K))) !DIF13
     &                    -(U(IPV(I),J,KMV(K))+U(I,J,KMV(K))))*SSDZ(K)
         ST(I,J,K,4)=0.5*((V(IPV(I),JPV(J),K)+V(IPV(I),J,K)) !DIF21
     &                    -(V(IMV(I),JPV(J),K)+V(IMV(I),J,K)))*SSDX(I)
         ST(I,J,K,5)=(V(I,JPV(J),K)-V(I,J,K))*SSDY(J)        !DIF22
         ST(I,J,K,6)=0.5*((V(I,JPV(J),KPV(K))+V(I,J,KPV(K))) !DIF23
     &                    -(V(I,JPV(J),KMV(K))+V(I,J,KMV(K))))*SSDZ(K)
         ST(I,J,K,7)=0.5*((W(IPV(I),J,KPV(K))+W(IPV(I),J,K)) !DIF31
     &                    -(W(IMV(I),J,KPV(K))+W(IMV(I),J,K)))*SSDX(I)
         ST(I,J,K,8)=0.5*((W(I,JPV(J),KPV(K))+W(I,JPV(J),K)) !DIF32
     &                    -(W(I,JMV(J),KPV(K))+W(I,JMV(J),K)))*SSDY(J)
         ST(I,J,K,9)=(W(I,J,KPV(K))-W(I,J,K))*SSDZ(K)        !DIF33
         
         DIF(I,J,K,1)=VIS*(ST(I,J,K,1)+ST(I,J,K,1))
         DIF(I,J,K,2)=VIS*(ST(I,J,K,2)+ST(I,J,K,4))
         DIF(I,J,K,3)=VIS*(ST(I,J,K,3)+ST(I,J,K,7))
         DIF(I,J,K,4)=VIS*(ST(I,J,K,5)+ST(I,J,K,5))
         DIF(I,J,K,5)=VIS*(ST(I,J,K,6)+ST(I,J,K,8))
         DIF(I,J,K,6)=VIS*(ST(I,J,K,9)+ST(I,J,K,9))
        ENDDO
        ENDDO
        ENDDO
        
       	CALL FILTERING(ST(0,0,0,1),ST_FILTER(0,0,0,1))
       	CALL FILTERING(ST(0,0,0,2),ST_FILTER(0,0,0,2))
       	CALL FILTERING(ST(0,0,0,3),ST_FILTER(0,0,0,3))
       	CALL FILTERING(ST(0,0,0,4),ST_FILTER(0,0,0,4))
       	CALL FILTERING(ST(0,0,0,5),ST_FILTER(0,0,0,5))
       	CALL FILTERING(ST(0,0,0,6),ST_FILTER(0,0,0,6))
       	CALL FILTERING(ST(0,0,0,7),ST_FILTER(0,0,0,7))
       	CALL FILTERING(ST(0,0,0,8),ST_FILTER(0,0,0,8))
       	CALL FILTERING(ST(0,0,0,9),ST_FILTER(0,0,0,9))
       	
       	CALL FILTERING(DIF(0,0,0,1),DIF_FILTER(0,0,0,1))
       	CALL FILTERING(DIF(0,0,0,2),DIF_FILTER(0,0,0,2))
       	CALL FILTERING(DIF(0,0,0,3),DIF_FILTER(0,0,0,3))
       	CALL FILTERING(DIF(0,0,0,4),DIF_FILTER(0,0,0,4))
       	CALL FILTERING(DIF(0,0,0,5),DIF_FILTER(0,0,0,5))
       	CALL FILTERING(DIF(0,0,0,6),DIF_FILTER(0,0,0,6))
       	
       	DO K=1,N3M
       	DO J=1,N2M
       	DO I=1,N1M
       	VIS=VISM+VIS_DIFF*PSI_FILTER(I,J,K)

       	TAU3(I,J,K,1)=DIF_FILTER(I,J,K,1)
     &-VIS*(ST_FILTER(I,J,K,1)+ST_FILTER(I,J,K,1))
       	TAU3(I,J,K,2)=DIF_FILTER(I,J,K,2)
     &-VIS*(ST_FILTER(I,J,K,2)+ST_FILTER(I,J,K,4))
       	TAU3(I,J,K,3)=DIF_FILTER(I,J,K,3)
     &-VIS*(ST_FILTER(I,J,K,3)+ST_FILTER(I,J,K,7))
       	TAU3(I,J,K,4)=DIF_FILTER(I,J,K,4)
     &-VIS*(ST_FILTER(I,J,K,5)+ST_FILTER(I,J,K,5))
       	TAU3(I,J,K,5)=DIF_FILTER(I,J,K,5)
     &-VIS*(ST_FILTER(I,J,K,6)+ST_FILTER(I,J,K,8))
       	TAU3(I,J,K,6)=DIF_FILTER(I,J,K,6)
     &-VIS*(ST_FILTER(I,J,K,9)+ST_FILTER(I,J,K,9))
       	
        ENDDO
        ENDDO
        ENDDO
      !diffusion term

       OPEN(42,FILE=outfile)
       WRITE(42,'(A300)') 'VARIABLES="X","Y","Z","psi_cn","body"
     & ,"tau2","tau3"'

        WRITE(42,*) 'ZONE I=',N1M/2,',J=',N2M/2,',K=',N3M/2,',F=POINT'
        
        TAU2_AVG=0.
        TAU3_AVG=0.
        AA=0.

       	DO K=2,N3M,2
       	DO J=2,N2M,2
       	DO I=2,N1M,2
       	 XX=XP(I)
         YY=YP(J)
         ZZ=ZP(K)

         VOLF=PSI_FILTER(I,J,K)
         BODY=FUNCBODY(XX,YY,ZZ)

        TAU2_SUM=
     &  -(0.5*(TAU2(IPV(IPV(I)),J,K,3)-TAU2(IMV(IMV(I)),J,K,3))
     &     /(VDX(IPV(I))+VDX(I))
     &  +0.5*(TAU2(I,JPV(JPV(J)),K,3)-TAU2(I,JMV(JMV(J)),K,3))
     &     /(VDY(JPV(J))+VDY(J))
     &  +0.5*(TAU2(I,J,KPV(KPV(K)),3)-TAU2(I,J,KMV(KMV(K)),3))
     &     /(VDZ(KPV(K))+VDZ(K)))
        TAU3_SUM=
     &  +(0.5*(TAU3(IPV(IPV(I)),J,K,3)-TAU3(IMV(IMV(I)),J,K,3))
     &     /(VDX(IPV(I))+VDX(I))
     &  +0.5*(TAU3(I,JPV(JPV(J)),K,3)-TAU3(I,JMV(JMV(J)),K,3))
     &     /(VDY(JPV(J))+VDY(J))
     &  +0.5*(TAU3(I,J,KPV(KPV(K)),3)-TAU3(I,J,KMV(KMV(K)),3))
     &     /(VDZ(KPV(K))+VDZ(K)))
     
       IF (ABS(TAU3_SUM) .GE. 0.1) THEN
       	DVOL=SDX(I)*SDY(J)*SDZ(K)
        TAU2_AVG=TAU2_AVG+TAU2_SUM**2*DVOL
        TAU3_AVG=TAU3_AVG+TAU3_SUM**2*DVOL
        AA=AA+DVOL
       ENDIF

         WRITE(42,47) XX,YY,ZZ,VOLF,BODY,TAU2_SUM,TAU3_SUM
       	ENDDO
        ENDDO
        ENDDO
 47     FORMAT(3F10.5,6E12.4)
 
        WRITE(*,*) 'L2(TAU2)=',sqrt(TAU2_AVG/AA)
        WRITE(*,*) 'L2(TAU3)=',sqrt(TAU3_AVG/AA)
 
       	DEALLOCATE(CONV,CONV_FILTER,DIF,DIF_FILTER,ST,ST_FILTER)
       	DEALLOCATE(U_FILTER,V_FILTER,W_FILTER,PSI_FILTER)
        DEALLOCATE(TAU2,TAU3)
        
       ENDIF

        WRITE(*,*) '3-D DATE IS MADE!!'

        RETURN
        END
       
C===================================================================
      SUBROUTINE FILTERING(A0,TMP)
C===================================================================
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      
      REAL CFX1(M1M,-1:1),CFY1(M2M,-1:1),CFZ1(M3M,-1:1)
      REAL A0(0:M1,0:M2,0:M3),TMP(0:M1,0:M2,0:M3)

!     Simpson's rule COEFFICIENT
      DO I=1,N1M
      H1=VDX(IMV(I)) !FOR Z-CV
      H2=VDX(I)
      CFX1(I,-1)=(2.*H1-H2)/6./H1
      CFX1(I, 0)=((H1+H2)**2)/6./H1/H2
      CFX1(I, 1)=(2.*H2-H1)/6./H2
      ENDDO

      DO J=1,N2M
      H1=VDY(JMV(J)) !FOR Z-CV
      H2=VDY(J)
      CFY1(J,-1)=(2.*H1-H2)/6./H1
      CFY1(J, 0)=((H1+H2)**2)/6./H1/H2
      CFY1(J, 1)=(2.*H2-H1)/6./H2
      ENDDO

      DO K=1,N3M
      H1=SDZ(KMV(K)) !FOR Z-CV
      H2=SDZ(K)
      CFZ1(K,-1)=(2.*H1-H2)/6./H1
      CFZ1(K, 0)=((H1+H2)**2)/6./H1/H2
      CFZ1(K, 1)=(2.*H2-H1)/6./H2
      ENDDO
      
      !FILTERING
      !X-DIR
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      TMP(I,J,K)=CFX1(I,-1)*A0(IMV(I),J,K)
     &          +CFX1(I, 0)*A0(I,J,K)
     &          +CFX1(I, 1)*A0(IPV(I),J,K)
      ENDDO
      ENDDO
      ENDDO

      IF (IPX .NE. 1) THEN
      DO K=1,N3M
      DO J=1,N2M
      TMP(1,J,K)=A0(1,J,K)
      TMP(N1M,J,K)=A0(N1M,J,K)
      ENDDO
      ENDDO
      ENDIF
      
      !Y-DIR
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      A0(I,J,K)=CFZ1(K  ,-1)*TMP(I,J,KMV(K))
     &         +CFZ1(K  , 0)*TMP(I,J,K  )
     &         +CFZ1(K  , 1)*TMP(I,J,KPV(K))
      ENDDO
      ENDDO
      ENDDO
      
      IF (IPY .NE. 1) THEN
      DO J=1,N2M
      DO I=1,N1M
      A0(I,J,1  )=TMP(I,J,1  )
      A0(I,J,N3M)=TMP(I,J,N3M)
      ENDDO
      ENDDO
      ENDIF

      !Z-DIR
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
      TMP(I,J,K)=CFX1(J,-1)*A0(I,JMV(J),K)
     &          +CFX1(J, 0)*A0(I,J,K)
     &          +CFX1(J, 1)*A0(I,JPV(J),K)
      ENDDO
      ENDDO
      ENDDO

      IF (IPY .NE. 1) THEN
      DO K=1,N3M
      DO I=1,N1M
      TMP(I,1,K)=A0(J,1,K)
      TMP(I,N2M,K)=A0(I,N2M,K)
      ENDDO
      ENDDO
      ENDIF

        RETURN
        END
        
C=================================================================
C  SECTION FOR VORTEX IDENTIFICATION OF JEONG & HUSSAIN
C=================================================================
      SUBROUTINE VORNLAMBDA2
      USE FLOW_VAR
      
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      COMMON/VELOCT/U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),
     &              W(0:M1,0:M2,0:M3)
      COMMON/EDGEUVW/UEDGE(M3M),VEDGE(M3M),WEDGE(M3M)
      COMMON/CENTER/AAA(0:M1,0:M2,0:M3),UC(0:M1,0:M2,0:M3),
     &              VC(0:M1,0:M2,0:M3),WC(0:M1,0:M2,0:M3),
     &              PC(0:M1,0:M2,0:M3),VORX(0:M1,0:M2,0:M3),
     &              VORY(0:M1,0:M2,0:M3),VORZ(0:M1,0:M2,0:M3),
     &              VORMAG(0:M1,0:M2,0:M3)
      COMMON/LAMDA2F/VLAMBDA2(0:M1,0:M2,0:M3)

      REAL SS2(3,3),OM2(3,3),EE(3),ZZ(3,3)
      REAL SS(3,3),OM(3,3),TMAT(3,3,M1M,M2M),DD(3),EVEC(3,3)

C  GET VELOCITY AT CELL CENTER => BOUNDARY TREATMENT
C  ---------------------------
      IF (IPZ .NE. 1) THEN
       DO I=1,N1M   ! K - direction
       DO J=1,N2M
        JP=JPV(J)
        UC(I,J,N3)=0.5*(U(I,J,N3)+U(IPV(I),J,N3))
        VC(I,J,N3)=0.5*(V(I,J,N3)+V(I,JP,N3))
        WC(I,J,N3)=W(I,J,N3)

        UC(I,J,0)=0.5*(U(I,J,0)+U(IPV(I),J,0))
        VC(I,J,0)=0.5*(V(I,J,0)+V(I,JP,0))
        WC(I,J,0)=W(I,J,1)

       END DO
       END DO
      ENDIF

      IF (IPY .NE. 1) THEN
       DO K=1,N3M   ! J - direction
        KP=KPV(K)
       ! UPPER & LOWER BOUNDARY
       DO I=1,N1M
        UC(I,N2,K)=0.5*(U(I,N2,K)+U(IPV(I),N2,K))
        VC(I,N2,K)=V(I,N2,K)
        WC(I,N2,K)=0.5*(W(I,N2,K)+W(I,N2,KP))

        JS=0
        UC(I,0,K)=0.5*(U(I,0,K)+U(IPV(I),0,K))
        VC(I,0,K)=V(I,1,K)
        WC(I,0,K)=0.5*(W(I,0,K)+W(I,0,KP))
       END DO
       END DO
      ENDIF

      IF (IPX .NE. 1) THEN
       DO K=1,N3M   ! I - direction
        KP=KPV(K)
       DO J=1,N2M
        JP=JPV(J)
        UC(N1,J,K)=U(N1,J,K)
        VC(N1,J,K)=0.5*(V(N1,J,K)+V(N1,JP,K))
        WC(N1,J,K)=0.5*(W(N1,J,K)+W(N1,J,KP))

        UC(0,J,K)=U(1,J,K)
        VC(0,J,K)=0.5*(V(0,J,K)+V(0,JP,K))
        WC(0,J,K)=0.5*(W(0,J,K)+W(0,J,KP))
       END DO
       END DO
      ENDIF

C  TO COMPUTE VORTICITY & LAMBDA2
      DO 30 K=1,N3M
       KP1=KPV(K)
       KM1=KMV(K)

      DO J=1,N2M
       JP1=JPV(J)
       JM1=JMV(J)
c       IG=IGEOM(J)
      DO I=1,N1M
       IP1=IPV(I)
       IM1=IMV(I)
c       JG=JGEOM(I)

       SH11=(U(IP1,J,K)-U(I,J,K))*SSDX(I)
       SH22=(V(I,JP1,K)-V(I,J,K))*SSDY(J)
       SH33=(W(I,J,KP1)-W(I,J,K))*SSDZ(K)

       UN=0.5*(UC(I,J,K)*SDY(JP1)+UC(I,JP1,K)*SDY(J))*VVDY(JP1)
       US=0.5*(UC(I,JM1,K)*SDY(J)+UC(I,J,K)*SDY(JM1))*VVDY(J)
       SH12=(UN-US)*SSDY(J)

       UT=0.5*(UC(I,J,K)*SDZ(KP1)+UC(I,J,KP1)*SDZ(K))*VVDZ(KP1)
       UB=0.5*(UC(I,J,KM1)*SDZ(K)+UC(I,J,K)*SDZ(KM1))*VVDZ(K)
       SH13=(UT-UB)*SSDZ(K)

       VE=0.5*(VC(I,J,K)*SDX(IP1)+VC(IP1,J,K)*SDX(I))*VVDX(IP1)
       VW=0.5*(VC(IM1,J,K)*SDX(I)+VC(I,J,K)*SDX(IM1))*VVDX(I)
       SH21=(VE-VW)*SSDX(I)

       VT=0.5*(VC(I,J,K)*SDZ(KP1)+VC(I,J,KP1)*SDZ(K))*VVDZ(KP1)
       VB=0.5*(VC(I,J,KM1)*SDZ(K)+VC(I,J,K)*SDZ(KM1))*VVDZ(K)
       SH23=(VT-VB)*SSDZ(K)

       WE=0.5*(WC(I,J,K)*SDX(IM1)+WC(IP1,J,K)*SDX(I))*VVDX(IP1)
       WW=0.5*(WC(IM1,J,K)*SDX(I)+WC(I,J,K)*SDX(IM1))*VVDX(I)
       SH31=(WE-WW)*SSDX(I)

       WN=0.5*(WC(I,J,K)*SDY(JP1)+WC(I,JP1,K)*SDY(J))*VVDY(JP1)
       WS=0.5*(WC(I,JM1,K)*SDY(J)+WC(I,J,K)*SDY(JM1))*VVDY(J)
       SH32=(WN-WS)*SSDY(J)

       VORX(I,J,K)=SH32-SH23
       VORY(I,J,K)=SH13-SH31
       VORZ(I,J,K)=SH21-SH12
       VORMAG(I,J,K)=SQRT(VORX(I,J,K)**2+VORY(I,J,K)**2+
     &                    VORZ(I,J,K)**2)

       SS(1,1)=SH11
       SS(1,2)=0.5*(SH12+SH21)
       SS(1,3)=0.5*(SH13+SH31)
       SS(2,1)=0.5*(SH21+SH12)
       SS(2,2)=SH22
       SS(2,3)=0.5*(SH23+SH32)
       SS(3,1)=0.5*(SH31+SH13)
       SS(3,2)=0.5*(SH32+SH23)
       SS(3,3)=SH33
       OM(1,1)=0.
       OM(1,2)=0.5*(SH12-SH21)
       OM(1,3)=0.5*(SH13-SH31)
       OM(2,1)=0.5*(SH21-SH12)
       OM(2,2)=0.
       OM(2,3)=0.5*(SH23-SH32)
       OM(3,1)=0.5*(SH31-SH13)
       OM(3,2)=0.5*(SH32-SH23)
       OM(3,3)=0.

       TMAT(1,1,I,J)=SS(1,1)*SS(1,1)+SS(1,2)*SS(2,1)+SS(1,3)*SS(3,1)
     &              +OM(1,1)*OM(1,1)+OM(1,2)*OM(2,1)+OM(1,3)*OM(3,1)
       TMAT(1,2,I,J)=SS(1,1)*SS(1,2)+SS(1,2)*SS(2,2)+SS(1,3)*SS(3,2)
     &              +OM(1,1)*OM(1,2)+OM(1,2)*OM(2,2)+OM(1,3)*OM(3,2)
       TMAT(1,3,I,J)=SS(1,1)*SS(1,3)+SS(1,2)*SS(2,3)+SS(1,3)*SS(3,3)
     &              +OM(1,1)*OM(1,3)+OM(1,2)*OM(2,3)+OM(1,3)*OM(3,3)
       TMAT(2,2,I,J)=SS(2,1)*SS(1,2)+SS(2,2)*SS(2,2)+SS(2,3)*SS(3,2)
     &              +OM(2,1)*OM(1,2)+OM(2,2)*OM(2,2)+OM(2,3)*OM(3,2)
       TMAT(2,3,I,J)=SS(2,1)*SS(1,3)+SS(2,2)*SS(2,3)+SS(2,3)*SS(3,3)
     &              +OM(2,1)*OM(1,3)+OM(2,2)*OM(2,3)+OM(2,3)*OM(3,3)
       TMAT(3,3,I,J)=SS(3,1)*SS(1,3)+SS(3,2)*SS(2,3)+SS(3,3)*SS(3,3)
     &              +OM(3,1)*OM(1,3)+OM(3,2)*OM(2,3)+OM(3,3)*OM(3,3)
       TMAT(2,1,I,J)=TMAT(1,2,I,J)
       TMAT(3,1,I,J)=TMAT(1,3,I,J)
       TMAT(3,2,I,J)=TMAT(2,3,I,J)
      END DO
      END DO

      DO J=1,N2M
      DO I=1,N1M
       CALL JACOBI(TMAT(1,1,I,J),3,3,DD,EVEC,ITER)
       CALL SORTER(3,DD)
       VLAMBDA2(I,J,K)=0.5*( DD(2)-ABS(DD(2)) )
      END DO
      END DO

   30 CONTINUE

      RETURN
      END

C--------------------------------------------------------------------
      SUBROUTINE SORTER(N,RA)

      REAL RA(3),RRA
      INTEGER I,IR,J,L

      IF(N.LE.1) GOTO 30
      L=N/2+1
      IR=N
   10 CONTINUE
      IF(L .GT. 1) THEN
            L=L-1
            RRA=RA(L)
         ELSE
            RRA=RA(IR)
            RA(IR)=RA(1)
            IR=IR-1
            IF(IR.EQ.1)THEN
                 RA(1)=RRA
                 RETURN
            END IF
         END IF
      I=L
      J=L+L
   20 IF(J.LE.IR) THEN
         IF(J.LT.IR) THEN
              IF(RA(J).GT.RA(J+1))J=J+1
         END IF
         IF(RRA.GT.RA(J)) THEN
              RA(I)=RA(J)
              I=J
              J=J+J
         ELSE
              J=IR+1
         END IF
      GOTO 20
      END IF
      RA(I)=RRA

      GOTO 10
   30 RETURN
      END

C--------------------------------------------------------------
C  ROUTINE TO COMPUTE EIGENVALUES OF MATRIX
C    - ROUTINE FROM NUMERICAL RECIPE IN FORTRAN
C--------------------------------------------------------------
      SUBROUTINE jacobi(a,n,np,d,v,nrot)
      INTEGER n,np,nrot,NMAX
      REAL a(np,np),d(np),v(np,np)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j
      REAL c,g,h,s,sm,t,tau,theta,tresh,b(NMAX),z(NMAX)
      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.
11      continue
        v(ip,ip)=1.
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))

14        continue
15      continue
        if(sm.eq.0.)return
        if(i.lt.4)then
          tresh=0.2*sm/n**2
        else
          tresh=0.
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+
     *g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h

              else
                theta=0.5*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.)t=-t
              endif
              c=1./sqrt(1+t**2)
              s=t*c
              tau=s/(1.+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)

                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)

                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      pause 'too many iterations in jacobi'
      return
      END

C*******************************************************************
      SUBROUTINE CAL_VOL(W,PSI_CN)
C*******************************************************************
      USE FLOW_VAR
      
      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR

      REAL W(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)
      
        OPEN(201,FILE='0POST.DAT',POSITION='APPEND')
!        WRITE(201,*) 'VARIABLES="TIME","VOL_M","VOL_P"'

      QM1=0.
      QM2=0.
      AVOL1=0.
      AVOL2=0.    
      DO K=1,N3M
      DO J=1,N2M
      DO I=1,N1M
        IF (FUNCBODY(XP(I),YP(J),Z(K)) .GE. 0.) THEN
          WC=0.5*(W(I,J,K)+W(I,J,KPV(K)))
          DVOL=SDX(I)*SDY(J)*VDZ(K)

!          IF (PSI_CN(I,J,K) .EQ. 0.) THEN
!          QM1=QM1+WC*DVOL
!          AVOL1=AVOL1+DVOL
!          ELSE
!          QM2=QM2+WC*DVOL
!          AVOL2=AVOL2+DVOL
!          ENDIF

          QM1=QM1+WC*(1.-PSI_CN(I,J,K))*DVOL
          QM2=QM2+WC*PSI_CN(I,J,K)*DVOL
          AVOL1=AVOL1+(1.-PSI_CN(I,J,K))*DVOL
          AVOL2=AVOL2+PSI_CN(I,J,K)*DVOL
        ENDIF
      ENDDO
      ENDDO
      ENDDO

        Q1=QM1/ZL
        Q2=QM2/ZL
        Q_TOT=Q1+Q2

        IF (AVOL1 .NE. 0.) THEN
         BULK1=QM1/AVOL1
        ELSE
         BULK1=0.
        ENDIF
        
        IF (AVOL2 .NE. 0.) THEN
         BULK2=QM2/AVOL2
        ELSE
         BULK2=0.
        ENDIF

C=====SAVE
        WRITE(201,200) TIME,AVOL1,AVOL2,Q1,Q2,Q_TOT,BULK1,BULK2
        WRITE(*,202) TIME,AVOL1,AVOL2
 200    FORMAT(F12.8,7ES15.6)
 202    FORMAT('TIME=',F12.8,' VOL1=',ES14.6,' VOL2=',ES14.6)

        CLOSE(201)

        RETURN
        END