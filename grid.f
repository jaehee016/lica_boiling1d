C*****************************************************************
C       Grid Generation Program (gridall.f)
C
C       Uniform or Non-uniform Using a HYPERBOLIC TANGENT Function
C
C       cf90 -o grid grid3d.f
C*****************************************************************
C       N1 & N2 & N3 : MUST BE AN ODD NUMBER

        PARAMETER(MGRID=1000000)
        PARAMETER(MSECT=100)
        COMMON /COORD/ XYZ(MGRID,6)
        DIMENSION IDIR(MSECT),XYZI(MSECT),XYZF(MSECT),NGI(MSECT)
     &           ,NGF(MSECT),IOPT(MSECT),FAC1(MSECT),FAC2(MSECT)
        CHARACTER*4 DUMMY
        REAL SIZE(MGRID,6)

C       IDIR() : GRID GENERATION DIRECTION
C       XYZI(): INITIAL COORDINATE
C       XYZF(): FINAL COORDINATE
C       NGI(): INITIAL GRID INDEX
C       NGF(): FINAL GRID INDEX
C       IOPT(): 1(UNIFORM), 2(GEOMETRICAL PROGRESSION)
C               3(HYPERBOLIC TANGENT)
C       FAC1(): ER or GAM
C       FAC2(): YESEXP or XSI

C------ read input file
        OPEN(11,FILE='grid.in')
        READ(11,*) DUMMY
        READ(11,*) ISPARE1,ISPARE2,ISPARE3,ISPARE4
        READ(11,*) DUMMY
        READ(11,*) DUMMY
        READ(11,*) NSECT1,N1,N2,N3,XL,YL,ZL
        READ(11,*) DUMMY
        DO L=1,NSECT1
        READ(11,*) IDIR(L),XYZI(L),XYZF(L),NGI(L),NGF(L),IOPT(L)
     &            ,FAC1(L),FAC2(L)
        ENDDO
        READ(11,*) DUMMY
        READ(11,*) DUMMY
        READ(11,*) NSECT2,N1G,N2G,N3G,XLG,YLG,ZLG
        READ(11,*) DUMMY
        DO L=NSECT1+1,NSECT1+NSECT2
        READ(11,*) IDIR(L),XYZI(L),XYZF(L),NGI(L),NGF(L),IOPT(L)
     &            ,FAC1(L),FAC2(L)
        ENDDO
            
!            WRITE(*,*)
!            RADIUS_R=XL
!            WRITE(*,139) XL/2
!            WRITE(*,140) YL/2
! 139   FORMAT('######RADIUS_R ALONG X IS ',F12.5)
! 140   FORMAT('######RADIUS_R ALONG Y IS ',F12.5)


        NSECT=NSECT1+NSECT2



       NGI_NEW1=0
       NGI_NEW2=0
       NGI_NEW3=0
       NGI_NEW4=0
       NGF_NEW1=0
       NGF_NEW2=0
       NGF_NEW3=0
       NGF_NEW4=0

C------ grid generation
        DO 300 L=1,NSECT
          !FOR IBM PIPE_B
        IF ( IDIR(L) .EQ. 1 .AND. (ISPARE1 .NE. 0) ) THEN
          IF (NGI(L) .EQ. 1) THEN
            NGI(L)=NGI(L)+ISPARE1
            NGI_NEW1=NGI(L)
          ENDIF
          IF (NGF(L) .EQ. N1) THEN
             NGF(L)=NGF(L)-ISPARE1
             NGF_NEW1=NGF(L)
          ENDIF
        ENDIF
        IF ( IDIR(L) .EQ. 2 .AND. (ISPARE2 .NE. 0) ) THEN
          IF (NGI(L) .EQ. 1) THEN
            NGI(L)=NGI(L)+ISPARE2
            NGI_NEW2=NGI(L)
          ENDIF
          IF (NGF(L) .EQ. N2) THEN
            NGF(L)=NGF(L)-ISPARE2
            NGF_NEW2=NGF(L)
          ENDIF
        ENDIF
        IF ( IDIR(L) .EQ. 4 .AND. (ISPARE3 .NE. 0) ) THEN
          IF (NGI(L) .EQ. 1) THEN
            NGI(L)=NGI(L)+ISPARE3
            NGI_NEW3=NGI(L)
          ENDIF
          IF (NGF(L) .EQ. N1G) THEN
             NGF(L)=NGF(L)-ISPARE3
             NGF_NEW3=NGF(L)
          ENDIF
        ENDIF
        IF ( IDIR(L) .EQ. 5 .AND. (ISPARE4 .NE. 0) ) THEN
          IF (NGI(L) .EQ. 1) THEN
            NGI(L)=NGI(L)+ISPARE4
            NGI_NEW4=NGI(L)
          ENDIF
          IF (NGF(L) .EQ. N2G) THEN
            NGF(L)=NGF(L)-ISPARE4
            NGF_NEW4=NGF(L)
          ENDIF
        ENDIF
          !FOR IBM PIPE_E

        IF (IOPT(L).EQ.1) THEN
           CALL UNIFORM(IDIR(L),XYZI(L),XYZF(L),NGI(L),NGF(L))
        ELSE IF (IOPT(L).EQ.2) THEN
           CALL GEOPRO(IDIR(L),XYZI(L),XYZF(L),NGI(L),NGF(L)
     &                ,FAC1(L),FAC2(L))
        ELSE
           CALL HYPTAN(IDIR(L),XYZI(L),XYZF(L),NGI(L),NGF(L)
     &                ,FAC1(L),FAC2(L))
        ENDIF

          !FOR IBM PIPE_B
        IF ( IDIR(L) .EQ. 1 .AND. ISPARE1 .NE. 0 ) THEN
         IF (NGI(L) .EQ. NGI_NEW1) THEN
           ID=IDIR(L)
           NGI(L)=NGI(L)-ISPARE1
           DO IS=ISPARE1,1,-1
          XYZ(NGI(L)+IS-1,ID)=XYZ(NGI(L)+IS,ID)
     &                  -(XYZ(NGI(L)+IS+1,ID)-XYZ(NGI(L)+IS,ID))**2
     &                        /(XYZ(NGI(L)+IS+2,ID)-XYZ(NGI(L)+IS+1,ID))
           ENDDO
         ENDIF
        IF (NGF(L) .EQ. NGF_NEW1) THEN
          ID=IDIR(L)
          NGF(L)=NGF(L)+ISPARE1
          DO IS=ISPARE1,1,-1
        XYZ(NGF(L)-IS+1,ID)=XYZ(NGF(L)-IS,ID)
     &            +(XYZ(NGF(L)-IS,ID)-XYZ(NGF(L)-IS-1,ID))**2
     &                      /(XYZ(NGF(L)-IS-1,ID)-XYZ(NGF(L)-IS-2,ID))
          ENDDO
         ENDIF
        
        ENDIF
        
        
       IF ( IDIR(L) .EQ. 2 .AND. ISPARE2 .NE. 0 ) THEN
         IF (NGI(L) .EQ. NGI_NEW2) THEN
          ID=IDIR(L)
          NGI(L)=NGI(L)-ISPARE2
          DO IS=ISPARE2,1,-1
        XYZ(NGI(L)+IS-1,ID)=XYZ(NGI(L)+IS,ID)
     &                  -(XYZ(NGI(L)+IS+1,ID)-XYZ(NGI(L)+IS,ID))**2
     &                        /(XYZ(NGI(L)+IS+2,ID)-XYZ(NGI(L)+IS+1,ID))
          ENDDO
        ENDIF
         IF (NGF(L) .EQ. NGF_NEW2) THEN
          ID=IDIR(L)
          NGF(L)=NGF(L)+ISPARE2
          DO IS=ISPARE2,1,-1
        XYZ(NGF(L)-IS+1,ID)=XYZ(NGF(L)-IS,ID)
     &            +(XYZ(NGF(L)-IS,ID)-XYZ(NGF(L)-IS-1,ID))**2
     &                       /(XYZ(NGF(L)-IS-1,ID)-XYZ(NGF(L)-IS-2,ID))
          ENDDO
        ENDIF
       ENDIF
        
        
       IF ( IDIR(L) .EQ. 4 .AND. ISPARE3 .NE. 0 ) THEN
        IF (NGI(L) .EQ. NGI_NEW3) THEN
          ID=IDIR(L)
          NGI(L)=NGI(L)-ISPARE3
          DO IS=ISPARE3,1,-1
        XYZ(NGI(L)+IS-1,ID)=XYZ(NGI(L)+IS,ID)
     &                  -(XYZ(NGI(L)+IS+1,ID)-XYZ(NGI(L)+IS,ID))**2
     &                        /(XYZ(NGI(L)+IS+2,ID)-XYZ(NGI(L)+IS+1,ID))
          ENDDO
         ENDIF
        IF (NGF(L) .EQ. NGF_NEW3) THEN
          ID=IDIR(L)
          NGF(L)=NGF(L)+ISPARE3
          DO IS=ISPARE3,1,-1
        XYZ(NGF(L)-IS+1,ID)=XYZ(NGF(L)-IS,ID)
     &            +(XYZ(NGF(L)-IS,ID)-XYZ(NGF(L)-IS-1,ID))**2
     &                       /(XYZ(NGF(L)-IS-1,ID)-XYZ(NGF(L)-IS-2,ID))
          ENDDO
        ENDIF
       ENDIF
        
        
        IF ( IDIR(L) .EQ. 5 .AND. ISPARE4 .NE. 0 ) THEN
        IF (NGI(L) .EQ. NGI_NEW4) THEN
          ID=IDIR(L)
          NGI(L)=NGI(L)-ISPARE4
          DO IS=ISPARE4,1,-1
        XYZ(NGI(L)+IS-1,ID)=XYZ(NGI(L)+IS,ID)
     &                  -(XYZ(NGI(L)+IS+1,ID)-XYZ(NGI(L)+IS,ID))**2
     &                        /(XYZ(NGI(L)+IS+2,ID)-XYZ(NGI(L)+IS+1,ID))
          ENDDO
        ENDIF
       IF (NGF(L) .EQ. NGF_NEW4) THEN
          ID=IDIR(L)
          NGF(L)=NGF(L)+ISPARE4
          DO IS=ISPARE4,1,-1
        XYZ(NGF(L)-IS+1,ID)=XYZ(NGF(L)-IS,ID)
     &            +(XYZ(NGF(L)-IS,ID)-XYZ(NGF(L)-IS-1,ID))**2
     &                       /(XYZ(NGF(L)-IS-1,ID)-XYZ(NGF(L)-IS-2,ID))
          ENDDO
        ENDIF
       ENDIF
          !FOR IBM PIPE_E
 300    CONTINUE

!        DO 410 I=1,N1-1
! 410       SIZE(I,1)=XYZ(I+1,1)-XYZ(I,1)
!        DO 420 I=1,N2-1
! 420       SIZE(I,2)=XYZ(I+1,2)-XYZ(I,2)
!        DO 430 I=1,N3-1
! 430       SIZE(I,3)=XYZ(I+1,3)-XYZ(I,3)
!        SIZE(N1,1)=SIZE(N1-1,1)
!        SIZE(N2,2)=SIZE(N2-1,2)
!        SIZE(N3,3)=SIZE(N3-1,3)

C------ print out
C       3-D
        OPEN(11,FILE='grid.dat')
        WRITE(11,111) N1,N2,N3
        WRITE(11,112) XL,YL,ZL
        WRITE(11,113) (XYZ(I,1),I=1,N1)
        WRITE(11,113) (XYZ(J,2),J=1,N2)
        WRITE(11,113) (XYZ(K,3),K=1,N3)
        WRITE(11,111) N1G,N2G,N3G
        WRITE(11,112) XLG,YLG,ZLG
        WRITE(11,113) (XYZ(I,4),I=1,N1G)
        WRITE(11,113) (XYZ(J,5),J=1,N2G)
        WRITE(11,113) (XYZ(K,6),K=1,N3G)
        CLOSE(11)
 111    FORMAT(3I23)
 112    FORMAT(3F23.15)
 113    FORMAT(5ES23.15)

!        OPEN(31,FILE='grid_xyzmp.dat')
!        WRITE(31,*)'XMP'
!        DO I=1,N1-1
!        WRITE(31,'(I5,F13.5)')I,0.5*(XYZ(I,1)+XYZ(I+1,1))
!        ENDDO
!        WRITE(31,*)'YMP'
!        DO J=1,N2-1
!        WRITE(31,'(I5,F13.5)')J,0.5*(XYZ(J,2)+XYZ(J+1,2))
!        ENDDO
!        WRITE(31,*)'ZMP'
!        DO K=1,N3-1
!        WRITE(31,'(I5,F13.5)')K,0.5*(XYZ(K,3)+XYZ(K+1,3))
!        ENDDO
!        CLOSE(31)
!
!        OPEN(31,FILE='grid_del.dat')
!        WRITE(31,*)'XMP'
!        DO I=1,N1-1
!        WRITE(31,'(I5,2F15.7)')I,XYZ(I,1),XYZ(I+1,1)-XYZ(I,1)
!        ENDDO
!        WRITE(31,*)'YMP'
!        DO J=1,N2-1
!        WRITE(31,'(I5,2F15.7)')J,XYZ(J,2),XYZ(J+1,2)-XYZ(J,2)
!        ENDDO
!        WRITE(31,*)'ZMP'
!        DO K=1,N3-1
!        WRITE(31,'(I5,2F15.7)')K,XYZ(K,3),XYZ(K+1,3)-XYZ(K,3)
!        ENDDO
!        CLOSE(31)

!        OPEN(12,FILE='grid.tec')
!C        WRITE(12,*) 'VARIABLES="X","Y"'
!C        WRITE(12,*) 'ZONE I=',N1,',J=',N2,',F=POINT'
!c        WRITE(12,121) ((XYZ(I,1),XYZ(J,2),I=1,N1),J=1,N2)
!c        WRITE(12,121) ((XYZ(I,1),XYZ(K,3),I=1,N1),K=1,N3)
!
!        WRITE(12,*) 'VARIABLES="X","Y","Z"'
!        WRITE(12,*) 'ZONE I=',N1,',J=',N2,',K=',N3,',F=POINT'
!        WRITE(12,121) (((XYZ(I,1),XYZ(J,2),XYZ(K,3)
!     &               ,I=1,N1),J=1,N2),K=1,N3)
! 121    FORMAT(3E13.5)
!        CLOSE(12)

        N1M=N1-1
        N2M=N2-1
        N3M=N3-1
!        ISKIP=1!4
!        OPEN(21,FILE='xy.tec')
!        WRITE(21,*) 'VARIABLES="x","y"'
!        WRITE(21,*) 'ZONE I=',N1M/ISKIP,',J=',N2M/ISKIP,',F=POINT'
!        WRITE(21,121) ((XYZ(I,1),XYZ(J,2),I=ISKIP,N1M,ISKIP)
!     &                                   ,J=ISKIP,N2M,ISKIP)
!
!        OPEN(22,FILE='xz.tec')
!        WRITE(22,*) 'VARIABLES="x","z"'
!        WRITE(22,*) 'ZONE I=',N1M/ISKIP,',J=',N3M/ISKIP,',F=POINT'
!        WRITE(22,121) ((XYZ(I,1),XYZ(K,3),I=ISKIP,N1M,ISKIP)
!     &                                   ,K=ISKIP,N3M,ISKIP)
!
! 121    FORMAT(3E13.5)

!        OPEN(13,FILE='grid.x')
!        WRITE(13,*) '   i   x    dx'
!        WRITE(13,131) (I,XYZ(I,1),SIZE(I,1),I=1,N1)
!        CLOSE(13)
! 131    FORMAT(I4,2F15.8)
!
!       OPEN(14,FILE='grid.y')
!       WRITE(14,*) '   j   y    dy'
!        WRITE(14,131) (J,XYZ(J,2),SIZE(J,2),J=1,N2)
!        CLOSE(14)
!
!        OPEN(15,FILE='grid.z')
!        WRITE(15,*) '   k   z    dz'
!        WRITE(15,131) (K,XYZ(K,3),SIZE(K,3),K=1,N3)
!        CLOSE(15)

        STOP
        END


C*****************************************************************
        SUBROUTINE UNIFORM(IDIR,XYZI,XYZF,NGI,NGF)
C*****************************************************************
C       Uniform Distribution
        PARAMETER(MGRID=1000000)
        COMMON /COORD/ XYZ(MGRID,6)

        DC=(XYZF-XYZI)/(NGF-NGI)
        DO 100 N=NGI,NGF
           XYZ(N,IDIR)=XYZI+(N-NGI)*DC
 100    CONTINUE
        
        WRITE(*,11) IDIR,DC
 11     FORMAT(I2,' Uniform ',F8.5)
 
c        OPEN(16,FILE='gridresults.dat',ACCESS='APPEND')
c        WRITE(16,12)IDIR,DC
c 12     FORMAT(I2,' Uniform',F8.5)
c        CLOSE(16)

        RETURN
        END

C*****************************************************************
        SUBROUTINE GEOPRO(IDIR,XYZI,XYZF,NGI,NGF,ER,YESEXP)
C*****************************************************************
C       Non-uniform Distribution Using a Geomerical Progression
C       Assume : given expansion (compression) ratio and grid number
C       ER: EXPANSION RATIO OR INVERSE OF COMPREESION RATIO
C       YESEXP: 1.(EXPANSION), 0.(COMPRESSION)

        PARAMETER(MGRID=1000000)
        COMMON /COORD/ XYZ(MGRID,6)

        IF (ER.EQ.1.) THEN
           PRINT*, 'Expansion Ratio is 1 !!!'
           STOP
        ENDIF
        IF (YESEXP.NE.1.) ER=1./ER

        DSI=(XYZF-XYZI)*(1.-ER)/(1.-ER**(NGF-NGI))
        DSF=DSI*ER**(NGF-NGI-1)

        DO 100 N=NGI+1,NGF-1
           XYZ(N,IDIR)=XYZI+DSI*(1.-ER**(N-NGI))/(1.-ER)
 100    CONTINUE
        XYZ(NGI,IDIR)=XYZI
        XYZ(NGF,IDIR)=XYZF

        WRITE(*,11) IDIR,DSI,DSF
 11     FORMAT(I2,' Geom. Progress.',F10.7,' (beg)',F10.7,' (end)')

c        OPEN(16,FILE='gridresults.dat',ACCESS='APPEND')
c        WRITE(16,13)IDIR,DSI,DSF
c 13     FORMAT(I2,' Geom. Progress.',F8.5,' (beg)',F8.5,' (end)')
c        CLOSE(16)


        RETURN
        END

C*****************************************************************
        SUBROUTINE HYPTAN(IDIR,XYZI,XYZF,NGI,NGF,GAM,XSI)
C*****************************************************************
C       Non-uniform Distribution Using a Hyperbolic Tangent Function
C       GAM: SLOPE (expansion ratio increase with GAM)
C           (expansion ratio gets larger with GAM)
C       XSI: INFLECTION POINT (grid clustering position)
C            1.0(flat plate), 0.5 (channel)
C            0.0(grid compression, opposite to the case of 1.0)

        PARAMETER(MGRID=1000000)
        COMMON /COORD/ XYZ(MGRID,6)

C------ when XSI=0.5 or 1.0
        IF (XSI .NE. 0.) THEN
        DO 100 N=NGI,NGF
           SEG=FLOAT(N-NGI)/FLOAT(NGF-NGI)
           XYZ(N,IDIR)=XYZI+(XYZF-XYZI)*XSI*(1.-TANH(GAM*(XSI-SEG))
     &                      /TANH(GAM*XSI))
 100    CONTINUE
        GO TO 333
        ENDIF

C------ when XSI=0.0, using a different function
        XSI=1.
        DO 200 N=NGI,NGF
           SEG=FLOAT(N-NGI)/FLOAT(NGF-NGI)
           XYZ(N,IDIR)=XYZI+(XYZF-XYZI)*XSI*TANH(GAM*SEG)
     &                      /TANH(GAM*XSI)
 200    CONTINUE

C------ print out
 333    DS0=XYZ(NGI+1,IDIR)-XYZ(NGI,IDIR)
        DS1=XYZ(NGI+(NGF-NGI)/2,IDIR)-XYZ(NGI+(NGF-NGI)/2-1,IDIR)
        DS2=XYZ(NGF,IDIR)-XYZ(NGF-1,IDIR)
        WRITE(*,11) IDIR,DS0,DS1,DS2
 11     FORMAT(I2,' Hypertan',F10.7,' (beg)',F10.7,' (mid)',
     &            F10.7,' (end)')

c        OPEN(16,FILE='gridresults.dat',ACCESS='APPEND')
c        Write(16,14)IDIR,DS0,DS1,DS2
c 14     FORMAT(I2,' Hypertan',F8.5,' (beg)',F8.5,' (mid)',F8.5,' (end)')
c        CLOSE(16)

        RETURN
        END


!C*******************************************************************
!        SUBROUTINE ADAPTIVE_GRID
!C*******************************************************************
!        INCLUDE 'param.h'
!        INCLUDE 'geom.h'
!       REAL DXCEN(100),DYCEN(100),DZCEN(100)
!       
!       CRI=0.5
!
!      DO 1000 ITER_GRID=1,1000
!
!        LX=0
!        LX=0
!        LZ=0
!       DO K=1,N3
!       DO J=1,N2
!       DO I=1,N1
!
!         !X-DIR
!        DX=ABS( FUNCBODY(XMP(I),YMP(J),ZMP(K),IS)
!     &                    -FUNCBODY(XMP(I-1),YMP(J),ZMP(K),IS) )*SSDX(I)
!       IF ( DX .GE. CRI ) THEN
!         LX=LX+1
!         DXCEN(LX)=0.5*(XMP(I)+XMP(I-1))
!       ENDIF
!
!        !Y-DIR
!        DY=ABS( FUNCBODY(XMP(I),YMP(J),ZMP(K),IS)
!     &                          -FUNCBODY(XMP(I),YMP(J-1),ZMP(K),IS) ) 
!       IF ( DY .GE. CRI ) THEN
!         LY=LY+1
!         DYCEN(LY)=0.5*(YMP(J)+YMP(J-1))
!       ENDIF
!
!        !Z-DIR
!        DZ=ABS( FUNCBODY(XMP(I),YMP(J),ZMP(K),IS)
!     &                          -FUNCBODY(XMP(I),YMP(J),ZMP(K-1),IS) ) 
!       IF ( DZ .GE. CRI ) THEN
!         LZ=LZ+1
!         DZCEN(LZ)=0.5*(ZMP(K)+ZMP(K-1))
!       ENDIF
!
!       ENDDO
!       ENDDO
!       ENDDO
!       
!        CALL MAKEGEOM(DXCEN,DYCEN,DZCEN,LX,LY,LZ)
!       
! 1000  CONTINUE
!
!       RETURN
!       END
!
!C*******************************************************************
!        SUBROUTINE MAKEGEOM(DXCEN,DYCEN,DZCEN,LX,LY,LZ)
!C*******************************************************************
!       REAL DXCEN(100),DYCEN(100),DZCEN(100)
!
!      OPEN(55,FILE='geom.h')
!      WRITE(55,*)'        COMMON/DIM/N1,N1M,N2,N2M,N3,N3M'
!      WRITE(55,*)'        COMMON/INDEXPM/IPV(0:M1),JPV(M2M),KPV(0:M3)
!     &                ,IMV(0:M1),JMV(M2M),KMV(0:M3)'
!      WRITE(55,*)'        COMMON/FIX/FIXJL(M2M),FIXJU(M2M)
!     &               ,FIXIL(M1M),FIXIU(M1M) ,FIXKL(M3M),FIXKU(M3M)'
!      WRITE(55,*)'        COMMON/SCALES/XL,YL,ZL'
!      WRITE(55,*)'        COMMON/GEOMINPUT/ND,NXF,NXB,NY'
!      WRITE(55,*)'        COMMON/INPUT/XFG,XBG,XFS,XBS,YG,YS'
!      WRITE(55,*)'        COMMON/COORD/X(0:M1),Y(0:M2),Z(0:M3)'
!      WRITE(55,*)'        COMMON/POSITION/XMP(0:M1),YMP(0:M2),ZMP(0:M3)'
!      WRITE(55,*)'        COMMON/VAR/SDX(0:M1),SDY(0:M2),SDZ(0:M3)
!     &            ,VDX(M1),VDY(M2),VDZ(M3)'
!      WRITE(55,*)'        COMMON/VVAR/SSDX(0:M1),SSDY(0:M2),SSDZ(0:M3)
!     &             ,VVDX(M1),VVDY(M2),VVDZ(M3)'
!      WRITE(55,*)'        COMMON/ROTATE/OMEG'
!      WRITE(55,*)'        COMMON/FINT/FINTX1(M1M),FINTX2(M1M)
!     &              ,FINTY1(M2M),FINTY2(M2M),FINTZ1(M3M),FINTZ2(M3M)'
!      CLOSE(55)
!
!       RETURN
!       END
!       
!C*******************************************************************
!        FUNCTION FUNCBODY(X,Y,Z,IS)
!C*******************************************************************
!C
!C       Function for an immersed body
!C       Required condition to use secant method :
!C          1. FUNCBODY(X,Y,Z)<0  inside the body
!C             FUNCBODY(X,Y,Z)=0  at the body
!C             FUNCBODY(X,Y,Z)>0  outside the body
!C          2. FUNCBODY Must be continuous function
!C          *It is sufficient to satisfy only 1. to finding forcing
!C           point but 2. is critical to use secant method
!C       Ex.
!C          FUNCBODY=X**2+Y**2+Z**2-0.5**2       ! for sphere
!C          FUNCBODY=X**2+Y**2-0.5**2            ! for infinite cylinder
!
!        FUNCBODY=-(Y**2.+Z**2.-0.5**2.)    ! for PIPE
!!        FUNCBODY=X**2.+Y**2.-0.5**2.      ! for CYLINDER
!
!        RETURN
!        END