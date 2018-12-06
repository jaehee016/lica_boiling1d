!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
      SUBROUTINE Lagprt(vol2,clvs,xx,yy,zz)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!      
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR
      
      USE LVS_VAR
      USE LVS_GEOM_VAR
      !use lsg_var !  하나 만드는게 관리가 쉬울 듯, 지금은 lvs_var 에 기생
      
      implicit none
      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)       
      ! REAL VOL2
      REAL XX,YY,zz,VOLF,VOL2,VOLI
      REAL X_CEN,Y_CEN,Z_CEN
      INTEGER CLVS
      integer i,j,k
 

       NPRT=NPRT+1
       WRITE(*,*) 'NUMBER OF LAGPRT=',NPRT
       

       X_CEN=XX
       Y_CEN=YY
       Z_CEN=ZZ 
      
      xlagprt(Nprt,1)=x_cen
      xlagprt(Nprt,2)=y_cen
      xlagprt(Nprt,3)=z_cen
      xlagprt(Nprt,4)=vol2
      !XLAGPRT(N,IDIR) IDIR=1,2,3 : POSITION, IDIR=4: VOLUME 
      !VLAGPRT(N,IDIR) IDIR=1,2,3 : DIRECTION 

       WRITE(*,*)'CHECK_VAL1=',xlagprt(Nprt,1)
       WRITE(*,*)'CHECK_VAL2=',xlagprt(Nprt,2)
       WRITE(*,*)'CHECK_VAL3=',xlagprt(Nprt,3)
       WRITE(*,*)'CHECK_VAL4=',xlagprt(Nprt,4)

      RETURN
      END
      
      
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!      
      SUBROUTINE LAGTRANS(u,v,w,VOL_TOT_ORI)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!            
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR
      
      USE LVS_VAR
      USE LVS_GEOM_VAR
      use LVS_COUPLING
      !use lsg_var !  하나 만드는게 관리가 쉬울 듯, 지금은 lvs_var 에 기생
      
      USE TWO_PHASE_PROPERTY
      
      implicit none
      REAL*8 PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)       
      REAL*8 VOL_TOT_ORI(MLVS)
      REAL*8 XX,YY,ZZ,VOLF,VOL2,VOLI
      REAL*8 X_CEN,Y_CEN,Z_CEN
      REAL*8 U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      ! REAL UF(M1F,M2F,M3F),VF(M1F,M2F,M3F),WF(M1F,M2F,M3F) 
      REAL*8 XPRT,YPRT,ZPRT,UPRT,VPRT,WPRT
     &             ,UPPRT,VPPRT,WPPRT
            !INTEGER CLVS
      INTEGER*8 I,J,K,LPRT,II,JJ,KK,IGEO,JGEO,KGEO,LL,IZ,IPRT,ipprt
      INTEGER*8 Iii,jjj,kkk,ib,jb,kb,lvs1,nn1
      INTEGER*8 ICC,IP1,IP2,
     &          JCC,JP1,JP2,
     &          KCC,KP1,KP2
      REAL*8  UX1,UY1,UZ1
     &     ,UX2,UY2,UZ2
     &     ,UX3,UY3,UZ3
     &     ,UX4,UY4,UZ4
     &     ,UX5,UY5,UZ5
     &     ,UX6,UY6,UZ6
     &     ,UX7,UY7,UZ7
     &     ,UX8,UY8,UZ8
      REAL*8 FXGEO,FYGEO,FZGEO   
      ! REAL*8 DUDX(3,10000,3)
      REAL*8 RER(3),CD(3),DRAG(3),DUDT1(3),DUDT2(3)
      REAL*8 buoyancy,frho,pi,rad, funcbody
      INTEGER*8 IDG1,IDG2,IDG3,IDG4,IDG5,IDG6
      CHARACTER*12 tname
      real*8 x1,x2,y1,y2,z1,z2,x0,y0,z0
      
      nn1=n_max-3
      
      if (nprt .eq. 0) goto 9090
      PI=acos(-1.)
      
      idg1=(laghist+1)/100000
      idg2=(laghist+1-idg1*100000)/10000
      idg3=(laghist+1-idg1*100000-idg2*10000)/1000
      idg4=(laghist+1-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(laghist+1-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6= laghist+1-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
      tname='LAGPRT'//char(idg1+48)//char(idg2+48)//
     &      char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
     
      
      ! initialize
      xprt=0.
      yprt=0.
      zprt=0.

  947 continue
      do lprt=1,nprt
      ! if the particle has gone into the lvs, xlagprt(lprt,4)=0
      if (xlagprt(lprt,4) .eq. 0.) then 
      ! write(*,*)lprt
      do iprt=lprt,nprt-1
      ipprt=iprt+1
      xlagprt(iprt,1)=xlagprt(ipprt,1)
      xlagprt(iprt,2)=xlagprt(ipprt,2)
      xlagprt(iprt,3)=xlagprt(ipprt,3)
      xlagprt(iprt,4)=xlagprt(ipprt,4)
      ulagprt(iprt,1)=ulagprt(ipprt,1)
      ulagprt(iprt,2)=ulagprt(ipprt,2)
      ulagprt(iprt,3)=ulagprt(ipprt,3)
      ulagprt(iprt,4)=ulagprt(ipprt,4)
      ulagprt(iprt,5)=ulagprt(ipprt,5)
      ulagprt(iprt,6)=ulagprt(ipprt,6)
      vlagprt(iprt,1)=vlagprt(ipprt,1)
      vlagprt(iprt,2)=vlagprt(ipprt,2)
      vlagprt(iprt,3)=vlagprt(ipprt,3)
      dvdt(iprt,1)=dvdt(ipprt,1)
      dvdt(iprt,2)=dvdt(ipprt,2)
      dvdt(iprt,3)=dvdt(ipprt,3)
      dvdt(iprt,4)=dvdt(ipprt,4)
      dvdt(iprt,5)=dvdt(ipprt,5)
      dvdt(iprt,6)=dvdt(ipprt,6)
      enddo
      nprt=nprt-1
      goto 947
      endif
            ! write(*,*)lprt, xlagprt(lprt,4)
      enddo

      ! before start, all the particle indices are pulled up.      
      
      do lprt=1,nprt
      ! write(*,*) 'particle number=',lprt

      
      xprt=xlagprt(lprt,1) ! position of k step
      yprt=xlagprt(lprt,2) ! position of k step
      zprt=xlagprt(lprt,3) ! position of k step
       ! write(*,*)'                                  '
       ! write(*,*)'check',xprt,yprt,zprt
  ! 616 format(3e13.5)
      
       ! WRITE(*,*)'START, NPRT=',LPRT
       ! WRITE(*,*)'CHECK_VAL1=',xlagprt(Lprt,1)
       ! WRITE(*,*)'CHECK_VAL2=',xlagprt(Lprt,2)
       ! WRITE(*,*)'CHECK_VAL3=',xlagprt(Lprt,3)
       ! WRITE(*,*)'CHECK_VAL4=',xlagprt(Lprt,4)
       ! WRITE(*,*)'GO ON!---------------------'      
      
      DO 100 II=1,N1M-1
      IF (XP(II+1) .GT. XPRT) THEN
      IGEO=II
      FXGEO=(XPRT-XP(II))/(XP(II+1)-XP(II))
      GOTO 110
      ENDIF
  100 CONTINUE 
  110 CONTINUE  
      
      DO 101 JJ=1,N2M-1
      IF (YP(JJ+1) .GT. YPRT) THEN 
      JGEO=JJ
      FYGEO=(YPRT-YP(JJ))/(YP(JJ+1)-YP(JJ))
      GOTO 111
      ENDIF
  101 CONTINUE 
  111 CONTINUE     
      
      DO 102 KK=1,N3M-1
      IF (ZP(KK+1) .GT. ZPRT) THEN
      KGEO=KK
      FZGEO=(ZPRT-ZP(KK))/(ZP(KK+1)-ZP(KK))
      GOTO 112
      ENDIF
      IF ((ZPRT .GT. ZP(N3M)) .and. (ZPRT .lt. ZL) )THEN
      KGEO=N3M
      FZGEO=(ZPRT-ZP(N3M))/(ZP(N3)-ZP(KK))
      endif
      
  102 CONTINUE 
  112 CONTINUE    
      
       ICC=IGEO
       IP1=1+ICC
       IP2=1+IP1
       x0=x(icc)
       x1=x(ip1)
       x2=x(ip2)
       JCC=JGEO
       JP1=1+JCC
       JP2=1+JP1
       y0=y(jcc)
       y1=y(jp1)
       y2=y(jp2)
       KCC=KGEO
       KP1=1+KCC
       KP2=1+KP1
       z0=z(kcc)
       z1=z(kp1)
       z2=z(kp2)
      IF ((ZPRT .GT. ZP(N3M)) .and. (ZPRT .lt. ZL) )THEN
      KCC=KGEO
      KP1=KPV(N3M)
      KP2=KP1+1
       z0=z(kcc)
       z1=z(kp1)
       z2=z(kp2)
      endif 
      ! write(*,*)'index',icc,jcc,kcc
!------ UX
! UX1~UX8: cell-center streamwise velocity       
      UX1=(U(ICC,JCC,KCC)+U(IP1,JCC,KCC))*0.5
      UX2=(U(IP1,JCC,KCC)+U(IP2,JCC,KCC))*0.5      
      UX3=(U(ICC,JP1,KCC)+U(IP1,JP1,KCC))*0.5
      UX4=(U(IP1,JP1,KCC)+U(IP2,JP1,KCC))*0.5      
      UX5=(U(ICC,JCC,KP1)+U(IP1,JCC,KP1))*0.5
      UX6=(U(IP1,JCC,KP1)+U(IP2,JCC,KP1))*0.5      
      UX7=(U(ICC,JP1,KP1)+U(IP1,JP1,KP1))*0.5
      UX8=(U(IP1,JP1,KP1)+U(IP2,JP1,KP1))*0.5          
      
      !ux1,ux2
      if     ( funcbody(x1,y0,z0) .le. 0.) then
            UX1=(U(ICC,JCC,KCC)+0.)*0.5     
           if (funcbody(x2,y0,z0) .le. 0.) then
            ux2=0.
           else
            ux2=(0. + u(ip2,jcc,kcc))*0.5
           endif
      else
           if (funcbody(x2,y0,z0) .le. 0.) then
            ux2=u(ip1,jcc,kcc)*0.5
           endif     
      endif
      !ux3,ux4
      if ( funcbody(x1,y1,z0) .le. 0.) then    
           if (funcbody(x0,y1,z0) .le. 0.) then 
            UX3=0.!(U(ICC,JP1,KCC)+U(IP1,JP1,KCC))*0.5
           else 
            ux3=u(icc,jp1,kcc)*0.5
           endif
            if (funcbody(x2,y1,z0) .le. 0.) then
            UX4=0.
           else 
            ux4=u(ip2,jp1,kcc)*0.5
           endif
      else
           if (funcbody(x0,y1,z0) .le. 0.) then
            UX3=U(Ip1,JP1,KCC)*0.5 !+U(IP1,JP1,KCC))*0.5
           endif
           if (funcbody(x2,y1,z0) .le. 0.) then 
            UX4=U(IP1,JP1,KCC)*0.5
           endif           
      endif
      !UX5,UX6
      if ( funcbody(x1,y0,z1) .le. 0.) then    
           if (funcbody(x0,y0,z1) .le. 0.) then 
            UX5=0.!               UX5=(U(ICC,JCC,KP1)+U(IP1,JCC,KP1))*0.5
           else 
            ux5=U(ICC,JCC,KP1)*0.5
           endif
            if (funcbody(x2,y0,z1) .le. 0.) then
            UX6=0.!      UX6=(U(IP1,JCC,KP1)+U(IP2,JCC,KP1))*0.5
           else 
            ux6=U(IP2,JCC,KP1)*0.5!      UX6=(U(IP1,JCC,KP1)+U(IP2,JCC,KP1))*0.5
           endif
      else
           if (funcbody(x0,y0,z1) .le. 0.) then
            UX5= U(IP1,JCC,KP1)*0.5
           endif
           if (funcbody(x2,y0,z1) .le. 0.) then 
            UX6= U(IP1,JCC,KP1)*0.5
           endif           
      endif      
      !     UX7=(U(ICC,JP1,KP1)+U(IP1,JP1,KP1))*0.5
      !     UX8=(U(IP1,JP1,KP1)+U(IP2,JP1,KP1))*0.5  
      if ( funcbody(x1,y1,z1) .le. 0.) then    
           if (funcbody(x0,y1,z1) .le. 0.) then
            UX7=0.!               UX5=(U(ICC,JCC,KP1)+U(IP1,JCC,KP1))*0.5
           else 
            ux7=U(ICC,Jp1,KP1)*0.5
           endif
            if (funcbody(x2,y1,z1) .le. 0.) then 
            UX8=0.!      UX6=(U(IP1,JCC,KP1)+U(IP2,JCC,KP1))*0.5
           else 
            ux8=U(IP2,Jp1,KP1)*0.5!      UX6=(U(IP1,JCC,KP1)+U(IP2,JCC,KP1))*0.5
           endif
      else
           if (funcbody(x0,y1,z1) .le. 0.) then 
            UX7= U(IP1,Jp1,KP1)*0.5
           endif
           if (funcbody(x2,y1,z1) .le. 0.) then 
            UX8= U(IP1,Jp1,KP1)*0.5
           endif           
      endif
      
      
      
       UPRT=(1.-FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UX1
     &     +(   FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UX2
     &     +(1.-FXGEO)*(   FYGEO)*(1.-FZGEO)*UX3
     &     +(   FXGEO)*(   FYGEO)*(1.-FZGEO)*UX4      
     &     +(1.-FXGEO)*(1.-FYGEO)*(   FZGEO)*UX5
     &     +(   FXGEO)*(1.-FYGEO)*(   FZGEO)*UX6
     &     +(1.-FXGEO)*(   FYGEO)*(   FZGEO)*UX7
     &     +(   FXGEO)*(   FYGEO)*(   FZGEO)*UX8     
     
      dudx(1,lprt,1)=(1.-FYGEO)*(1.-FZGEO)*(UX2-UX1)/(XP(II+1)-XP(II))
     &              +(   FYGEO)*(1.-FZGEO)*(UX4-UX3)/(XP(II+1)-XP(II))
     &              +(1.-FYGEO)*(   FZGEO)*(UX6-UX5)/(XP(II+1)-XP(II))
     &              +(   FYGEO)*(   FZGEO)*(UX8-UX7)/(XP(II+1)-XP(II))
        
      dudx(2,lprt,1)=(1.-FXGEO)*(1.-FZGEO)*(UX3-UX1)/(YP(JJ+1)-YP(JJ)) 
     &              +(   FXGEO)*(1.-FZGEO)*(UX4-UX2)/(YP(JJ+1)-YP(JJ)) 
     &              +(1.-FXGEO)*(   FZGEO)*(UX7-UX5)/(YP(JJ+1)-YP(JJ)) 
     &              +(   FXGEO)*(   FZGEO)*(UX8-UX6)/(YP(JJ+1)-YP(JJ)) 
     
      dudx(3,lprt,1)=(1.-FXGEO)*(1.-FYGEO)*(UX5-UX1)/(ZP(KK+1)-ZP(KK)) 
     &              +(   FXGEO)*(1.-FYGEO)*(UX6-UX2)/(ZP(KK+1)-ZP(KK)) 
     &              +(1.-FXGEO)*(   FYGEO)*(UX7-UX3)/(ZP(KK+1)-ZP(KK)) 
     &              +(   FXGEO)*(   FYGEO)*(UX8-UX4)/(ZP(KK+1)-ZP(KK))

      
      
       IF (JCC.NE.0) THEN
       UY1=(V(ICC,JCC,KCC)+V(ICC,JP1,KCC))*0.5
       UY2=(V(IP1,JCC,KCC)+V(IP1,JP1,KCC))*0.5
       UY3=(V(ICC,JP1,KCC)+V(ICC,JP2,KCC))*0.5
       UY4=(V(IP1,JP1,KCC)+V(IP1,JP2,KCC))*0.5
       UY5=(V(ICC,JCC,KP1)+V(ICC,JP1,KP1))*0.5
       UY6=(V(IP1,JCC,KP1)+V(IP1,JP1,KP1))*0.5
       UY7=(V(ICC,JP1,KP1)+V(ICC,JP2,KP1))*0.5
       UY8=(V(IP1,JP1,KP1)+V(IP1,JP2,KP1))*0.5
       
       ! UY1=(V(ICC,JCC,KCC)+V(ICC,JP1,KCC))*0.5
! UY3=(V(ICC,JP1,KCC)+V(ICC,JP2,KCC))*0.5
      if     ( funcbody(x0,y1,z0) .le. 0.) then
            uy1=(v(ICC,JCC,KCC)+0.)*0.5     
           if (funcbody(x0,y2,z0) .le. 0.) then
            uy3=0.
           else
            uy3=(0. + v(icc,jp2,kcc))*0.5
           endif
      else
           if (funcbody(x2,y0,z0) .le. 0.) then
            uy3=v(icc,jp1,kcc)*0.5
           endif     
      endif
! UY2=(V(IP1,JCC,KCC)+V(IP1,JP1,KCC))*0.5
! UY4=(V(IP1,JP1,KCC)+V(IP1,JP2,KCC))*0.5
      if ( funcbody(x1,y1,z0) .le. 0.) then    
           if (funcbody(x1,y0,z0) .le. 0.) then 
            uy2=0.!(v(ICC,JP1,KCC)+v(IP1,JP1,KCC))*0.5
           else 
            uy2=v(ip1,jcc,kcc)*0.5
           endif
            if (funcbody(x1,y2,z0) .le. 0.) then
            uy4=0.
           else 
            uy4=v(ip1,jp2,kcc)*0.5
           endif
      else
           if (funcbody(x1,y0,z0) .le. 0.) then
            uy2=v(Ip1,JP1,KCC)*0.5 !+v(IP1,JP1,KCC))*0.5
           endif
           if (funcbody(x1,y2,z0) .le. 0.) then 
            uy4=v(IP1,JP1,KCC)*0.5
           endif           
      endif
! UY5=(V(ICC,JCC,KP1)+V(ICC,JP1,KP1))*0.5
! UY7=(V(ICC,JP1,KP1)+V(ICC,JP2,KP1))*0.5
      if ( funcbody(x0,y1,z1) .le. 0.) then    
           if (funcbody(x0,y0,z1) .le. 0.) then 
            uy5=0.!               uy5=(v(ICC,JCC,KP1)+v(IP1,JCC,KP1))*0.5
           else 
            uy5=v(ICC,JCC,KP1)*0.5
           endif
            if (funcbody(x0,y2,z1) .le. 0.) then
            uy7=0.!      uy6=(v(IP1,JCC,KP1)+v(IP2,JCC,KP1))*0.5
           else 
            uy7=v(Icc,Jp2,KP1)*0.5!      uy6=(v(IP1,JCC,KP1)+v(IP2,JCC,KP1))*0.5
           endif
      else
           if (funcbody(x0,y0,z1) .le. 0.) then
            uy5= v(Icc,Jp1,KP1)*0.5
           endif
           if (funcbody(x0,y2,z1) .le. 0.) then 
            uy7= v(Icc,Jp1,KP1)*0.5
           endif           
      endif      
! UY6=(V(IP1,JCC,KP1)+V(IP1,JP1,KP1))*0.5
! UY8=(V(IP1,JP1,KP1)+V(IP1,JP2,KP1))*0.5
      if ( funcbody(x1,y1,z1) .le. 0.) then    
           if (funcbody(x1,y0,z1) .le. 0.) then
            uy6=0.!               uy5=(v(ICC,JCC,KP1)+v(IP1,JCC,KP1))*0.5
           else 
            uy6=v(Ip1,Jcc,KP1)*0.5
           endif
            if (funcbody(x2,y1,z1) .le. 0.) then 
            uy8=0.!      uy6=(v(IP1,JCC,KP1)+v(IP2,JCC,KP1))*0.5
           else 
            uy8=v(Ip1,Jp2,KP1)*0.5!      uy6=(v(IP1,JCC,KP1)+v(IP2,JCC,KP1))*0.5
           endif
      else
           if (funcbody(x1,y0,z1) .le. 0.) then 
            uy6= v(IP1,Jp1,KP1)*0.5
           endif
           if (funcbody(x2,y1,z1) .le. 0.) then 
            uy8= v(IP1,Jp1,KP1)*0.5
           endif           
      endif
       
       
       ELSE
!C##############
       UY1= V(ICC,1,KCC)
       UY2= V(IP1,1,KCC)
       UY3=(V(ICC,1,KCC)+V(ICC,2,KCC))*0.5
       UY4=(V(IP1,1,KCC)+V(IP1,2,KCC))*0.5
       UY5= V(ICC,1,KP1) 
       UY6= V(IP1,1,KP1) 
       UY7=(V(ICC,1,KP1)+V(ICC,2,KP1))*0.5
       UY8=(V(IP1,1,KP1)+V(IP1,2,KP1))*0.5
       ENDIF    
       

       
       VPRT=(1.-FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UY1
     &     +(   FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UY2
     &     +(1.-FXGEO)*(   FYGEO)*(1.-FZGEO)*UY3
     &     +(   FXGEO)*(   FYGEO)*(1.-FZGEO)*UY4
     &     +(1.-FXGEO)*(1.-FYGEO)*(   FZGEO)*UY5
     &     +(   FXGEO)*(1.-FYGEO)*(   FZGEO)*UY6
     &     +(1.-FXGEO)*(   FYGEO)*(   FZGEO)*UY7
     &     +(   FXGEO)*(   FYGEO)*(   FZGEO)*UY8

      dudx(1,lprt,2)=(1.-FYGEO)*(1.-FZGEO)*(UY2-UY1)/(XP(II+1)-XP(II))
     &              +(   FYGEO)*(1.-FZGEO)*(UY4-UY3)/(XP(II+1)-XP(II))
     &              +(1.-FYGEO)*(   FZGEO)*(UY6-UY5)/(XP(II+1)-XP(II))
     &              +(   FYGEO)*(   FZGEO)*(UY8-UY7)/(XP(II+1)-XP(II))
                                             
      dudx(2,lprt,2)=(1.-FXGEO)*(1.-FZGEO)*(UY3-UY1)/(YP(JJ+1)-YP(JJ)) 
     &              +(   FXGEO)*(1.-FZGEO)*(UY4-UY2)/(YP(JJ+1)-YP(JJ)) 
     &              +(1.-FXGEO)*(   FZGEO)*(UY7-UY5)/(YP(JJ+1)-YP(JJ)) 
     &              +(   FXGEO)*(   FZGEO)*(UY8-UY6)/(YP(JJ+1)-YP(JJ)) 
                                           
      dudx(3,lprt,2)=(1.-FXGEO)*(1.-FYGEO)*(UY5-UY1)/(ZP(KK+1)-ZP(KK)) 
     &              +(   FXGEO)*(1.-FYGEO)*(UY6-UY2)/(ZP(KK+1)-ZP(KK)) 
     &              +(1.-FXGEO)*(   FYGEO)*(UY7-UY3)/(ZP(KK+1)-ZP(KK)) 
     &              +(   FXGEO)*(   FYGEO)*(UY8-UY4)/(ZP(KK+1)-ZP(KK))
     
     

     
!C------ UT
!C UT1~UT8: cell-center azimuthal velocity
       UZ1=(W(ICC,JCC,KCC)+W(ICC,JCC,KP1))*0.5
       UZ2=(W(IP1,JCC,KCC)+W(IP1,JCC,KP1))*0.5
       UZ3=(W(ICC,JP1,KCC)+W(ICC,JP1,KP1))*0.5
       UZ4=(W(IP1,JP1,KCC)+W(IP1,JP1,KP1))*0.5
       UZ5=(W(ICC,JCC,KP1)+W(ICC,JCC,KP2))*0.5
       UZ6=(W(IP1,JCC,KP1)+W(IP1,JCC,KP2))*0.5
       UZ7=(W(ICC,JP1,KP1)+W(ICC,JP1,KP2))*0.5
       UZ8=(W(IP1,JP1,KP1)+W(IP1,JP1,KP2))*0.5
       WPRT=(1.-FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UZ1
     &     +(   FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UZ2
     &     +(1.-FXGEO)*(   FYGEO)*(1.-FZGEO)*UZ3
     &     +(   FXGEO)*(   FYGEO)*(1.-FZGEO)*UZ4
     &     +(1.-FXGEO)*(1.-FYGEO)*(   FZGEO)*UZ5
     &     +(   FXGEO)*(1.-FYGEO)*(   FZGEO)*UZ6
     &     +(1.-FXGEO)*(   FYGEO)*(   FZGEO)*UZ7
     &     +(   FXGEO)*(   FYGEO)*(   FZGEO)*UZ8
     
      dudx(1,lprt,3)=(1.-FYGEO)*(1.-FZGEO)*(UZ2-UZ1)/(XP(II+1)-XP(II))
     &              +(   FYGEO)*(1.-FZGEO)*(UZ4-UZ3)/(XP(II+1)-XP(II))
     &              +(1.-FYGEO)*(   FZGEO)*(UZ6-UZ5)/(XP(II+1)-XP(II))
     &              +(   FYGEO)*(   FZGEO)*(UZ8-UZ7)/(XP(II+1)-XP(II))
                                            
      dudx(2,lprt,3)=(1.-FXGEO)*(1.-FZGEO)*(UZ3-UZ1)/(YP(JJ+1)-YP(JJ)) 
     &              +(   FXGEO)*(1.-FZGEO)*(UZ4-UZ2)/(YP(JJ+1)-YP(JJ)) 
     &              +(1.-FXGEO)*(   FZGEO)*(UZ7-UZ5)/(YP(JJ+1)-YP(JJ)) 
     &              +(   FXGEO)*(   FZGEO)*(UZ8-UZ6)/(YP(JJ+1)-YP(JJ)) 
                                        
      dudx(3,lprt,3)=(1.-FXGEO)*(1.-FYGEO)*(UZ5-UZ1)/(ZP(KK+1)-ZP(KK)) 
     &              +(   FXGEO)*(1.-FYGEO)*(UZ6-UZ2)/(ZP(KK+1)-ZP(KK)) 
     &              +(1.-FXGEO)*(   FYGEO)*(UZ7-UZ3)/(ZP(KK+1)-ZP(KK)) 
     &              +(   FXGEO)*(   FYGEO)*(UZ8-UZ4)/(ZP(KK+1)-ZP(KK))

                                                                  !지금은..나중엔 (k+1)
      ulagprt(lprt,1)=uPRt  ! k위치에서 k+1의 속도를 interpolate, 4은 끝에 k위치에서 k의 속도 interpolate
      ulagprt(lprt,2)=vPRt  ! k위치에서 k+1의 속도를 interpolate, 5은 끝에 k위치에서 k의 속도 interpolate
      ulagprt(lprt,3)=wPRt  ! k위치에서 k+1의 속도를 interpolate, 6은 끝에 k위치에서 k의 속도 interpolate
      !________________________________________________________!
      
      if ((ulagprt(lprt,4) .eq. 0.) .and.
     &    (ulagprt(lprt,5) .eq. 0.) .and.
     &    (ulagprt(lprt,6) .eq. 0.))then
      
      VLAGPRT(LPRT,1)=UPRT+1.D-6    ! first step after becoming particle
      VLAGPRT(LPRT,2)=VPRT+1.D-6    ! 완전히 같게 하면 뒤에서 탈이 나서...floating divided by zero
      VLAGPRT(LPRT,3)=WPRT+1.D-6       
      ULAGPRT(LPRT,4)=UPRT    ! first step after becoming particle
      ULAGPRT(LPRT,5)=VPRT
      ULAGPRT(LPRT,6)=WPRT  
      endif
      
      RAD=(XLAGPRT(LPRT,4)*0.75/PI)**(1./3.)
      iii=icou2(icc)
      jjj=jcou2(jcc)
      kkk=kcou2(kcc)
      ! write(*,*)'another index', iii,jjj,kkk
      ! write(*,*)'?????why??????'
      ! write(*,*)uPRt,ulagprt(lprt,4)
      ! write(*,*)vPRt,ulagprt(lprt,5)
      ! write(*,*)wPRt,ulagprt(lprt,6)
      ! write(*,*)'??????????????'
      !in case the bubble goes into the lvs but not in the activated region                                                
      if (psif(iii,jjj,kkk) .eq. 1.) then 
      ! write(*,*)iii,jjj,kkk,'surrounding is water'
      FRHO=DENR/(1.+0.5*DENR)
      elseif (psif(iii,jjj,kkk) .eq. 0.) then 
      write(*,*)iii,jjj,kkk,'surrounding is air'
      LVS1=MASK_BUB(icc,jcc,kcc)
      VOL_TOT_ORI(LVS1)=VOL_TOT_ORI(LVS1)+XLAGPRT(LPRT,4)
      XLAGPRT(LPRT,4)=0.
      CALL GLOBAL_MASS_CORRECTION(NN1,LVS1,VOL_TOT_ORI,dtconst)
      GOTO 444
      ! frho=2d0/3d0
      endif
      ! BUOYANCY
      BUOYANCY=(1.-1.5*FRHO)*GRAVITY

      RER=0.
      CD=0.
      DRAG=0.
      
      ! ! ! ! STOKES DRAG  
      ! ! ! ur=dsqrt(ulagprt(lprt,4)**2+ulagprt(lprt,5)**2+ulagprt(lprt,6)**2)
      ! ! ! vr=dsqrt(vlagprt(lprt,1)**2+vlagprt(lprt,2)**2+vlagprt(lprt,3)**2)
      ! ! ! if (psif(iii,jjj,kkk) .eq. 1.) then 
      ! ! ! RER=2.*RAD*rep*dABS(ur-vr)
      ! ! ! else
      ! ! ! write(*,*)iii,jjj,kkk,'surrounding is air-2'
      ! ! ! rer=2.*RAD*dABS(ur-vr)
      ! ! ! endif   
     
      ! ! ! CD=1+0.15*RER**(0.687)+0.0175/(1.+42500.*RER**(-1.16))  
      
      DO LL=1,3
      ! STOKES DRAG  
      ! Kim (2009) for stability
      cd(ll)=24d0/rep*(1d0+0.15d0*rer(ll)**0.687d0)
      drag(ll)=frho*0.375d0*cd(ll)/rad
     & *abs(ULAGPRT(LPRT,LL+3)-VLAGPRT(LPRT,LL))
     & *(ULAGPRT(LPRT,LL+3)-VLAGPRT(LPRT,LL))
       ! RER(LL)=(2./rep)*RAD*ABS(ULAGPRT(LPRT,LL+3)-VLAGPRT(LPRT,LL))
       ! if (rer(ll) .lt. 1d3) then 
       ! fff=1d0*((rer(ll)**(2d0/3d0))/6d0)
       ! elseif (rer(ll) .lt. 3d5) then
       ! fff=0.0183d0*rer(ll)
       ! elseif (rer(ll) .gt. 3d5) then
       ! fff=0.0183d0*300000d0
       ! endif
       ! drag(ll)=(1d0/rep)*frho*4.5d0*((1d0/rad)**2d0)
      ! &      *(ULAGPRT(LPRT,LL+3)-VLAGPRT(LPRT,LL))
      ! &      *fff
      
      ! INERTIA
      DUDT1(LL)=FRHO*((ULAGPRT(LPRT,LL)-ULAGPRT(LPRT,LL+3))/DT
     &               + ULAGPRT(LPRT,4)*DUDX(4,LPRT,LL)
     &               + ULAGPRT(LPRT,5)*DUDX(5,LPRT,LL)
     &               + ULAGPRT(LPRT,6)*DUDX(6,LPRT,LL))     
      !ADDED MASS: ROTATIONAL EFFECT NEGLECTED FOR THE MOMENT
      DUDT2(LL)=0.5*FRHO*((ULAGPRT(LPRT,LL)-ULAGPRT(LPRT,LL+3))/DT
     &                   + VLAGPRT(LPRT,1)*DUDX(4,LPRT,LL)
     &                   + VLAGPRT(LPRT,2)*DUDX(5,LPRT,LL)
     &                   + VLAGPRT(LPRT,3)*DUDX(6,LPRT,LL))
      
      DUDX(4,LPRT,LL)=DUDX(1,LPRT,LL)
      DUDX(5,LPRT,LL)=DUDX(2,LPRT,LL)
      DUDX(6,LPRT,LL)=DUDX(3,LPRT,LL)
      ! if (lprt .eq. 15) then
!       write(*,*) '                        '!lprt,ulagprt(lprt,ll+3),RER(LL),vlagprt(lprt,ll)
!       write(*,*) ulagprt(lprt,ll+3),vlagprt(lprt,ll)
!       write(*,*) '-----', rer(ll)!,cd(ll)  
!       write(*,*) 'radius',rad
!       write(*,*) 'drag',ll,drag(ll)
!       write(*,*) 'uper(-,456),vlag',ulagprt(lprt,ll+3),vlagprt(lprt,ll)
!       write(*,*) 'uper(-,123),vlag',ulagprt(lprt,ll  ),vlagprt(lprt,ll)
!       write(*,*) 'what is wrong',DUDX(1,LPRT,LL)
!       write(*,*) 'what is wrong',DUDX(2,LPRT,LL)
!       write(*,*) 'what is wrong',DUDX(3,LPRT,LL)       
!       write(*,*) 'dudt1,2',      DUDT1(ll),DUDT2(ll)
!       write(*,*) 'dvdt',         dvdt(lprt,ll)
       ! endif      
      
      ENDDO
      
      !DVDT=ABDPRT
      DVDT(LPRT,1)=(          DRAG(1)+DUDT1(1)+DUDT2(1))
      DVDT(LPRT,2)=(          DRAG(2)+DUDT1(2)+DUDT2(2))
      DVDT(LPRT,3)=(-BUOYANCY+DRAG(3)+DUDT1(3)+DUDT2(3))

      IF ((DVDT(LPRT,4) .EQ. 0.) .AND.
     &    (DVDT(LPRT,5) .EQ. 0.) .AND.
     &    (DVDT(LPRT,6) .EQ. 0.))THEN
      ! FIRST STEP AFTER BECOMING PARTICLE : explicit euler
      UPPRT=VLAGPRT(LPRT,1)+DT*DVDT(LPRT,1)
      VPPRT=VLAGPRT(LPRT,2)+DT*DVDT(LPRT,2)
      WPPRT=VLAGPRT(LPRT,3)+DT*DVDT(LPRT,3) 
      else
      ! if not, RK3
      UPPRT=VLAGPRT(LPRT,1)
     &      +DT*(GAMMA(MSUB)*DVDT(LPRT,1)+RO(MSUB)*DVDT(LPRT,4))
      VPPRT=VLAGPRT(LPRT,2)
     &      +DT*(GAMMA(MSUB)*DVDT(LPRT,2)+RO(MSUB)*DVDT(LPRT,5))
      WPPRT=VLAGPRT(LPRT,3)
     &      +DT*(GAMMA(MSUB)*DVDT(LPRT,3)+RO(MSUB)*DVDT(LPRT,6))      
      endif
      
         ! write(*,*)'--------------------------------------------'
         ! write(*,*)DVDT(LPRT,1)
         ! write(*,*)DVDT(LPRT,2)
         ! write(*,*)DVDT(LPRT,3)
         ! write(*,*)upprt,uprt
         ! write(*,*)vpprt,vprt
         ! write(*,*)wpprt,wprt     
         ! write(*,*)'--------------------------------------------'
      
      !TO USE AT NEXT STEP
      DVDT(LPRT,4)=DVDT(LPRT,1)
      DVDT(LPRT,5)=DVDT(LPRT,2)
      DVDT(LPRT,6)=DVDT(LPRT,3)

      !________________________________________________________!
      !   below part doesn't need to be changed. dx_p/dt=V_p; RK3 step
      !  Vp(n)
      !  Vp(*)=VlagPRT(lprt,1)
      !  Vp(**)=uprt,vprt,wprt;
      !__Vp(n)=GAMMA(MSUB)*DT*uPRT+RO(MSUB)*DT*VlagPRT(lprt,1)_! (particle 위치정보가 k+1 이 아니라 k 임..
      
      XPRT=XPRT+dt*(GAMMA(MSUB)*UPPRT+RO(MSUB)*VLAGPRT(LPRT,1))
      YPRT=YPRT+dt*(GAMMA(MSUB)*VPPRT+RO(MSUB)*VLAGPRT(LPRT,2))
      ZPRT=ZPRT+dt*(GAMMA(MSUB)*WPPRT+RO(MSUB)*VLAGPRT(LPRT,3))

      IF ((ZPRT .GT. ZL) .AND. (IPZ .EQ. 1)) then 
      Iz= zprt/ZL
      write(*,*) 'z periodic',iz,zprt,zl
      ZPRT=ZPRT-iz*ZL
      write(*,*) 'z periodic_zprt=',zprt
      elseif ((zprt .lt. 0.) .and. (ipz .eq. 1)) then 
      zprt=zl+zprt
      endif 
      
      ! write(*,*)lprt,xprt,yprt,zprt,funcbody(xprt,yprt,zprt)      
      if (FUNCBODY(XPRT,YPRT,ZPRT) .LT. 0.) THEN   ! 벽 밖으로 나가면 그만큼 튕겨 안쪽으로 넣기
!      write(*,*)'out!!',lprt,xprt,yprt,zprt,funcbody(xprt,yprt,zprt)
      XPRT=XPRT-2d0*(dSQRT(XPRT**2d0+YPRT**2d0)-0.5d0*YL)
     &        *(XPRT/dSQRT(XPRT**2d0+YPRT**2d0))
      YPRT=YPRT-2d0*(dSQRT(XPRT**2d0+YPRT**2d0)-0.5d0*YL)
     &        *(YPRT/dSQRT(XPRT**2d0+YPRT**2d0))
!        write(*,*)'out,then in',
!     &             lprt,xprt,yprt,zprt,funcbody(xprt,yprt,zprt)
      ENDIF

      !-----------------------------------------------------------------
      !-----------------------------------------------------------------
      !-----------여기에 -particle이 lvs과 겹치지 않는지 확인하는-------
      !-----------서브루틴을 불러온다. 만약 겹치면 그녀석은 제명...-----
      !-----------그러면 순서가 또 개판이 되는데...---------------------
      ! call LAGLVS(XPRT,YPRT,ZPRT,LPRT) MOVED TO TWO_PHASE.F DUE TO VOL_TOT_ORI
      !-----------------------------------------------------------------
      
      XLAGPRT(LPRT,1)=XPRT
      XLAGPRT(LPRT,2)=YPRT
      XLAGPRT(LPRT,3)=ZPRT
      
      VLAGPRT(LPRT,1)=UPPRT
      VLAGPRT(LPRT,2)=VPPRT
      VLAGPRT(LPRT,3)=WPPRT      
      

  444 CONTINUE
       enddo
!C------ print out (Tecplot format)
!C       all particles

      IF ((MOD(NTIME,NPRINT).EQ.0).and. (msub .eq. 3) ) then !.and. (msub .eq. 3) )THEN
        OPEN(51,FILE=tname//'.DAT')
      WRITE(51,510) laghist+1, NPRT
      DO LPRT=1,NPRT
      WRITE(51,511) xlagprt(Lprt,1),xlagprt(Lprt,2),
     &              xlagprt(Lprt,3),xlagprt(Lprt,4),
     &              vlagprt(Lprt,1),vlagprt(Lprt,2),
     &              vlagprt(Lprt,3)
      endDO
      ENDIF
 510   FORMAT(2I8)     
 511   FORMAT(7E13.5)   

 
      IF ((MOD(NTIME,NPRINT).EQ.0) .and. (msub .eq. 3) ) THEN
        OPEN(52,FILE=tname//'tec.DAT')
      WRITE(52,*) 'variables="X","Y","Z","vol"'
      DO LPRT=1,NPRT
      WRITE(52,511) xlagprt(Lprt,1),xlagprt(Lprt,2),
     &              xlagprt(Lprt,3),xlagprt(Lprt,4)
      endDO
      ENDIF    
      
      ! DO LPRT=1,NPRT 
      ! WRITE(*,*) xlagprt(Lprt,1),xlagprt(Lprt,2),
      ! &              xlagprt(Lprt,3),xlagprt(Lprt,4),
      ! &              vlagprt(Lprt,1),vlagprt(Lprt,2),
      ! &              vlagprt(Lprt,3)
      ! endDO     
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
! for the next step.

      do lprt=1,nprt
      xprt=xlagprt(lprt,1) ! position of k+1 step
      yprt=xlagprt(lprt,2) ! position of k+1 step
      zprt=xlagprt(lprt,3) ! position of k+1 step
  
      DO 200 II=1,N1M-1
      IF (XP(II+1) .GT. XPRT) THEN
      IGEO=II
      FXGEO=(XPRT-XP(II))/(XP(II+1)-XP(II))
      GOTO 210
      ENDIF
  200 CONTINUE 
  210 CONTINUE  
      
      DO 201 JJ=1,N2M-1
      IF (YP(JJ+1) .GT. YPRT) THEN 
      JGEO=JJ
      FYGEO=(YPRT-YP(JJ))/(YP(JJ+1)-YP(JJ))
      GOTO 211
      ENDIF
  201 CONTINUE 
  211 CONTINUE     
      
      DO 202 KK=1,N3m-1
      IF (ZP(KK+1) .GT. ZPRT) THEN
      KGEO=KK
      FZGEO=(ZPRT-ZP(KK))/(ZP(KK+1)-ZP(KK))
      GOTO 212
      ENDIF
      IF ((ZPRT .GT. ZP(N3M)) .and. (ZPRT .lt. ZL) )THEN
      KGEO=N3M
      FZGEO=(ZPRT-ZP(N3M))/(ZP(N3)-ZP(KK))
      endif
  202 CONTINUE 
  212 CONTINUE    
      
       ICC=IGEO
       IP1=1+ICC
       IP2=1+IP1
       JCC=JGEO
       JP1=1+JCC
       JP2=1+JP1
       KCC=KGEO
       KP1=1+KCC
       KP2=1+KP1
      
!------ UX
! UX1~UX8: cell-center streamwise velocity
      UX1=(U(ICC,JCC,KCC)+U(IP1,JCC,KCC))*0.5
      UX2=(U(IP1,JCC,KCC)+U(IP2,JCC,KCC))*0.5
      UX3=(U(ICC,JP1,KCC)+U(IP1,JP1,KCC))*0.5
      UX4=(U(IP1,JP1,KCC)+U(IP2,JP1,KCC))*0.5
      UX5=(U(ICC,JCC,KP1)+U(IP1,JCC,KP1))*0.5
      UX6=(U(IP1,JCC,KP1)+U(IP2,JCC,KP1))*0.5
      UX7=(U(ICC,JP1,KP1)+U(IP1,JP1,KP1))*0.5
      UX8=(U(IP1,JP1,KP1)+U(IP2,JP1,KP1))*0.5
       UPRT=(1.-FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UX1
     &     +(   FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UX2
     &     +(1.-FXGEO)*(   FYGEO)*(1.-FZGEO)*UX3
     &     +(   FXGEO)*(   FYGEO)*(1.-FZGEO)*UX4      
     &     +(1.-FXGEO)*(1.-FYGEO)*(   FZGEO)*UX5
     &     +(   FXGEO)*(1.-FYGEO)*(   FZGEO)*UX6
     &     +(1.-FXGEO)*(   FYGEO)*(   FZGEO)*UX7
     &     +(   FXGEO)*(   FYGEO)*(   FZGEO)*UX8
      
!C------ UR
       IF (JCC.NE.0) THEN
       UY1=(V(ICC,JCC,KCC)+V(ICC,JP1,KCC))*0.5
       UY2=(V(IP1,JCC,KCC)+V(IP1,JP1,KCC))*0.5
       UY3=(V(ICC,JP1,KCC)+V(ICC,JP2,KCC))*0.5
       UY4=(V(IP1,JP1,KCC)+V(IP1,JP2,KCC))*0.5
       UY5=(V(ICC,JCC,KP1)+V(ICC,JP1,KP1))*0.5
       UY6=(V(IP1,JCC,KP1)+V(IP1,JP1,KP1))*0.5
       UY7=(V(ICC,JP1,KP1)+V(ICC,JP2,KP1))*0.5
       UY8=(V(IP1,JP1,KP1)+V(IP1,JP2,KP1))*0.5
       ELSE
!C##############
!C CAUTION: UR is not defined when J=0
       UY1= V(ICC,1,KCC)
       UY2= V(IP1,1,KCC)
       UY3=(V(ICC,1,KCC)+V(ICC,2,KCC))*0.5
       UY4=(V(IP1,1,KCC)+V(IP1,2,KCC))*0.5
       UY5= V(ICC,1,KP1) 
       UY6= V(IP1,1,KP1) 
       UY7=(V(ICC,1,KP1)+V(ICC,2,KP1))*0.5
       UY8=(V(IP1,1,KP1)+V(IP1,2,KP1))*0.5
       ENDIF
       
       VPRT=(1.-FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UY1
     &     +(   FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UY2
     &     +(1.-FXGEO)*(   FYGEO)*(1.-FZGEO)*UY3
     &     +(   FXGEO)*(   FYGEO)*(1.-FZGEO)*UY4
     &     +(1.-FXGEO)*(1.-FYGEO)*(   FZGEO)*UY5
     &     +(   FXGEO)*(1.-FYGEO)*(   FZGEO)*UY6
     &     +(1.-FXGEO)*(   FYGEO)*(   FZGEO)*UY7
     &     +(   FXGEO)*(   FYGEO)*(   FZGEO)*UY8      
!C------ UT
!C UT1~UT8: cell-center azimuthal velocity
       UZ1=(W(ICC,JCC,KCC)+W(ICC,JCC,KP1))*0.5
       UZ2=(W(IP1,JCC,KCC)+W(IP1,JCC,KP1))*0.5
       UZ3=(W(ICC,JP1,KCC)+W(ICC,JP1,KP1))*0.5
       UZ4=(W(IP1,JP1,KCC)+W(IP1,JP1,KP1))*0.5
       UZ5=(W(ICC,JCC,KP1)+W(ICC,JCC,KP2))*0.5
       UZ6=(W(IP1,JCC,KP1)+W(IP1,JCC,KP2))*0.5
       UZ7=(W(ICC,JP1,KP1)+W(ICC,JP1,KP2))*0.5
       UZ8=(W(IP1,JP1,KP1)+W(IP1,JP1,KP2))*0.5
       WPRT=(1.-FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UZ1
     &     +(   FXGEO)*(1.-FYGEO)*(1.-FZGEO)*UZ2
     &     +(1.-FXGEO)*(   FYGEO)*(1.-FZGEO)*UZ3
     &     +(   FXGEO)*(   FYGEO)*(1.-FZGEO)*UZ4
     &     +(1.-FXGEO)*(1.-FYGEO)*(   FZGEO)*UZ5
     &     +(   FXGEO)*(1.-FYGEO)*(   FZGEO)*UZ6
     &     +(1.-FXGEO)*(   FYGEO)*(   FZGEO)*UZ7
     &     +(   FXGEO)*(   FYGEO)*(   FZGEO)*UZ8 
!C----save interpolated fluid velocity at the new particle position     
      ulagprt(lprt,4)=uPRt  !k+1위치에서 k+1 속도를 interpolate
      ulagprt(lprt,5)=vPRt  !k+1위치에서 k+1 속도를 interpolate
      ulagprt(lprt,6)=wPRt  !k+1위치에서 k+1 속도를 interpolate      
       enddo
      
 9090 continue
      
      RETURN
      END

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!      
      SUBROUTINE READLAG(ihist)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
      USE LVS_VAR
      USE LVS_GEOM_VAR
      !use lsg_var !  하나 만드는게 관리가 쉬울 듯, 지금은 lvs_var 에 기생
      INTEGER IDG1,IDG2,IDG3,IDG4,IDG5,IDG6,ihist
      CHARACTER*12 tname

      idg1=(ihist)/100000
      idg2=(ihist-idg1*100000)/10000
      idg3=(ihist-idg1*100000-idg2*10000)/1000
      idg4=(ihist-idg1*100000-idg2*10000-idg3*1000)/100
      idg5=(ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100)/10
      idg6= ihist-idg1*100000-idg2*10000-idg3*1000-idg4*100-idg5*10
      tname='LAGPRT'//char(idg1+48)//char(idg2+48)//
     &      char(idg3+48)//char(idg4+48)//char(idg5+48)//char(idg6+48)
     
      write(*,*) tname
      OPEN(51,FILE=tname//'.DAT')
      READ(51,510) laghist, NPRT
      DO LPRT=1,NPRT
      READ(51,511) xlagprt(Lprt,1),xlagprt(Lprt,2),
     &             xlagprt(Lprt,3),xlagprt(Lprt,4),
     &             vlagprt(Lprt,1),vlagprt(Lprt,2),
     &             vlagprt(Lprt,3)
      endDO
      write(*,*)nprt,'particles read'
 510   FORMAT(2I8)  
 511    FORMAT(7E13.5)      
      
      return
      end
    
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
      SUBROUTINE BREAKUP(LLVS,VOL_TOT_ORI)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX!
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR
      
      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE
      
      REAL D_CRI,TEMP                                         
      REAL DX,DY,DZ,DVOL,DIVOL
      REAL VOL_TEMP
      INTEGER MASK_TAG(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      INTEGER      TAG(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      INTEGER LEVEL,LLEVEL,LVS1,LLVS,L
      INTEGER I,J,K,II,JJ,KK,ISS,JSS,KSS,I2,J2,K2
      INTEGER Itmp,jtmp,ktmp,II2,JJ2,KK2
      ! INTEGER IQ,JQ,KQ,IQQ,JQQ,KQQ

      INTEGER MASK_LVS(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      INTEGER   MASK_S(-2:M1F+3,-2:M2F+3,-2:M3F+3)   
      INTEGER   ICHECK(-2:M1F+3,-2:M2F+3,-2:M3F+3)   
      INTEGER N,NN,NN1,NN2,NM,NC1
      INTEGER NUMS,NUMS_INITIAL,NUMS_CRI      
      INTEGER NUMS_CHECK    

      INTEGER NUMA_TMP(N_MAX)
      INTEGER I_TMP(MF_BAND,N_MAX)
      INTEGER J_TMP(MF_BAND,N_MAX)
      INTEGER K_TMP(MF_BAND,N_MAX)
      
      REAL XX,YY,ZZ
      REAL ALPHI_TMP(-2:M1F+3,-2:M2F+3,-2:M3F+3)      
      ! REAL ALPHI_2XT(-2:M1F+3,-2:M2F+3,-2:M3F+3)     
      REAL ALPHI_3XT(-2:M1F+3,-2:M2F+3,-2:M3F+3)     
      REAL VOL_TOT_ORI(MLVS)
      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)   
      REAL PSIF2(M1L:M1U,M2L:M2U,M3L:M3U)   
      ! REAL ALP1(-2:M1F+3,-2:M2F+3,-2:M3F+3)      
      ! REAL ALP2(-2:M1F+3,-2:M2F+3,-2:M3F+3)      
      ! REAL ALPT(-2:M1F+3,-2:M2F+3,-2:M3F+3) 
      INTEGER max_Tag,max_maskT
      INTEGER NS,NC,IS,IC
      INTEGER I_S(MF_BAND),I_C(MF_BAND),
     &        J_S(MF_BAND),J_C(MF_BAND),
     &        K_S(MF_BAND),K_C(MF_BAND) 
      real XXX(1000),YYY(1000),ZZZ(1000)
      
      INTEGER nlvs1
      REAL VOL_TMPS(0:1000), VOL_TOT_TMP,VOLF,VOLI
      
      nn1=n_max-3
      nn2=n_max-1 

!XXXXXXXXXXXXXXXXXXXXXXXXXXX dVOLUME & dLENGTH XXXXXXXXXXXXXXXXXXXXXXXXX
        IF (N3FM .EQ. 1) THEN
        DVOL=SDXF*SDYF
        ELSE
        DVOL=SDXF*SDYF*SDZF
        ENDIF
        
        DX=10.
        DY=10.
        DZ=10.
        DO K=1,N3M
        DO J=1,N2M
        DO I=1,N1M
         DZ=MIN(DZ,SDZ(k))
         DX=MIN(DX,SDX(I))
         DY=MIN(DY,SDY(j))    
        ENDDO
        ENDDO
        ENDDO 
        DIVOL =27.*DX*DY*DZ !8.*(DX*DY*DZ)=64*sdxf*sdyf*sdzf
               
        D_CRI=1.2*MAX(SDXF,SDYF,SDZF)  
!XXXXXXXXXXXXXXXXXXXXXXXXXX dVOLUME & dLENGTH XXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
!XXXXXXXXXXXXX INITIALIZE VARIABLES & SAVE ORIGINAL ALPHI_EXT XXXXXXXXXX
!$OMP PARALLEL DO      
       DO K=M3L,M3U 
       DO J=M2L,M2U        
       DO I=M1L,M1U 
       PSIF_MTP(I,J,K)=0.
       PSIF2(I,J,K)=0.
       ENDDO
       ENDDO
       ENDDO           

      CALL GLOBAL_MASS_CORRECTION(NN1,LLVS,VOL_TOT_ORI,dtconst)   ! 이거 안하면 손실 된 채로 시작한다. 
      
      ! CALCULATE PSIF_MTP FOR THE CURRENT LLVS: 
      ! ALPHI로 작업할 수 없는 이유는, 모든 영역에 정의되지 않음 (BUBBLE 완전 안쪽이라던지)
      ! ALPHI_EXT 값을 사용해서 구하기 때문에 비우기 전에 구해놔야 함
      ! 그래서 IOUTER=1 : 모든 영역에 대해서 구해야하니까...      
      CALL CAL_PSIF_MTP(1,NN2,LLVS,PSIF_MTP)        
!$OMP PARALLEL DO      
       DO K=M3L,M3U 
       DO J=M2L,M2U        
       DO I=M1L,M1U 
       PSIF2(I,J,K)=PSIF_MTP(I,J,K)
       ENDDO
       ENDDO
       ENDDO  
 
 
!$OMP PARALLEL DO            
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       MASK_TAG(I,J,K)=0       !
            TAG(I,J,K)=0       ! 구분짓는 작업 중 지나갔는지 확인하는 TAG
       MASK_LVS(I,J,K)=0       ! 작업하는 ORIGINAL LVS 의 ACTIVATED GRID
         MASK_S(I,J,K)=0       !
         
       ALPHI_TMP(I,J,K)=ALPHI_EXT(I,J,K,1)
       ALPHI_3XT(I,J,K)=0.
       ALPHI_EXT(i,j,k,1)=0.
      ENDDO
      ENDDO
      ENDDO
         
      
      ! ALPHI_EXT=ALPHI(LLVS), TMP VALUES ARE MADE  
!$OMP PARALLEL DO private(N,NN,I,J,K)      
      DO N=1,N_MAX
      NUMA_TMP(N)=NUMA(N,LLVS)  
      DO NN=1,NUMA(N,LLVS)     
      I_TMP(NN,N)=I_B(NN,N,LLVS)
      J_TMP(NN,N)=J_B(NN,N,LLVS)
      K_TMP(NN,N)=K_B(NN,N,LLVS) 
      I=I_B(NN,N,LLVS)
      J=J_B(NN,N,LLVS)
      K=K_B(NN,N,LLVS)     
      MASK_LVS(I,J,K)=1 ! MASK_LVS COVERS THE ORIGINAL LVS FIELD
      ENDDO
      ENDDO 

!XXXXXXXXXXXXXX INITIALIZE VARIABLES & SAVE ORIGINAL ALPHI_EXT XXXXXXXXXX      
!XXXXXXXXXXXXXXXXXXXXXXXXX TAGGING DIFFERENT LEVEL XXXXXXXXXXXXXXXXXXXXXX 
      LEVEL=1
      NS=0                     
      NC=0  
       
      XX=0.
      YY=0.
      ZZ=0.
      VOLF=0.
      vol_tot_tmp=0.     

       
      DO 503 N=1,N_MAX         
      DO 503 NN=1,NUMA_TMP(N) ! TEMP. VAL. because the llvs index changes with lvs1
      I=I_TMP(NN,N)           ! 일단 모든 ACTIVATED 된 영역에 대해 한다.
      J=J_TMP(NN,N)           ! neighbor를 tag 하기 때문에 activated 만 해도 됨
      K=K_TMP(NN,N)            
   
       
       IF ((PSIF_MTP(I,J,K) .LT. 0.999) .AND. (TAG(I,J,K) .EQ. 0)) THEN !만약 INTERFACE (+내부)인데 TAG=0
       TAG(I,J,K)=LEVEL  ! 그 INTERFACE CELL은 LEVEL 로 TAG
       NS=NS+1           ! 그 INTERFACE CELL COUNT
       I_S(NS)=I         ! NS 번째 I_S의 값으로 LIBRARY화
       J_S(NS)=J
       K_S(NS)=K
       
          VOLF=(1.-PSIF_MTP(I,J,K))*SDXF*SDYF*SDZF
          XX=XX+XPF(I)*VOLF
          YY=YY+YPF(J)*VOLF
          ZZ=ZZ+ZPF(K)*VOLF
          vol_tot_tmp=vol_tot_tmp+VOLF 
        
  501  CONTINUE  
        ! 정의된 NS CELL들에 대해서 +1,-1 씩 주변을 확인 
        DO IS=1,NS
        I2=I_S(IS)
        J2=J_S(IS)
        K2=K_S(IS)

!----------------------------------------------------------------------NS cell 에 대해 NC를 정의         
        DO 502 KK2=K2-1,K2+1  
        KK=KK2
          IF (IPZ .EQ. 1) THEN
          IF (KK .EQ. 0) KK=N3FM
          IF (KK .EQ. N3F) KK=1  
          ENDIF          
        DO 502 JJ2=J2-1,J2+1
        JJ=JJ2   
          IF (IPY .EQ. 1) THEN
          IF (JJ .EQ. 0) JJ=N2FM
          IF (JJ .EQ. N2F) JJ=1  
          ENDIF      
        DO 502 II2=I2-1,I2+1
        II=II2
          IF (IPX .EQ. 1) THEN
          IF (II .EQ. 0) II=N1FM
          IF (II .EQ. N1F) II=1  
          ENDIF 
      
      ! TEST 중인 NS 그리드의 주변으로 INTERFACE인데 TAG=0이면 NC 로 카운트하고 TAG.
      IF ((PSIF_MTP(II,JJ,KK) .LT. 0.999).AND.(TAG(II,JJ,KK).EQ.0)) THEN
      TAG(II,JJ,KK)=LEVEL
      NC=NC+1
      I_C(NC)=II
      J_C(NC)=JJ
      K_C(NC)=KK
          
          VOLF=(1.-PSIF_MTP(II,JJ,KK))*SDXF*SDYF*SDZF
          XX=XX+XPF(II)*VOLF
          YY=YY+YPF(JJ)*VOLF
          ZZ=ZZ+ZPF(KK)*VOLF
          vol_tot_tmp=vol_tot_tmp+VOLF
 
      ENDIF
      
  502   CONTINUE
        ENDDO! IS=1,NS 이걸 다 돌아야 I_S=I_C 를 정의하지..
      ! NS에 대해 모두 수행하면, 이제 ...
!----------------------------------------------------------------------      
      !NC 를 NS 에 넣어준다.           
      DO IC=1,NC
      I_S(IC)=I_C(IC)
      J_S(IC)=J_C(IC)
      K_S(IC)=K_C(IC)
      ENDDO
      NS=NC
      NC=0     
      
      IF (NS .eq. 0) THEN
        vol_tmps(level)= vol_tot_tmp  
        IF (VOL_TOT_TMP .NE. 0.) THEN
        VOLI=1./VOL_TOT_TMP
        XXX(level)=XX*VOLI
        YYY(level)=YY*VOLI
        ZZZ(level)=ZZ*VOLI
        
        XX=0.
        YY=0.
        ZZ=0.
        VOLF=0.
        vol_tot_tmp=0. ! 여기는 하나의 level 판별이 끝나야 다시 올 수 있음      
        ENDIF
        LEVEL=LEVEL+1
         GOTO 503
      ELSE ! (NS .ne. 0) 
         GOTO 501
      ENDIF         
      ENDIF ! FIRST TAG RELATED IF
  503 CONTINUE ! the big tagging loop
                   
      LEVEL=LEVEL-1 ! level+1 at the last step, even if there is nothing to tag with (level+1)       
      
!XXXXXXXXXXXXXXXXXXXXXXXXXXX IF NO LEVEL-SET BREAKUP XXXXXXXXXXXXXXXXXXXX 
      IF (LEVEL .LE. 1) then
      !write(*,*)'----------no llvs breakup------------'
!$OMP PARALLEL DO private(I,J)             
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
      ALPHI_EXT(i,j,k,1)=ALPHI_TMP(I,J,K) ! if there is no breakup, restore the alphi_ext
      ENDDO
      ENDDO
      ENDDO  
      vol_tmps(level)=0.
      vol_tot_tmp=0.      
      GOTO 777 
      ENDIF            
!XXXXXXXXXXXXXXXXXX IF THERE IS LEVEL-SET BREAKUP XXXXXXXXXXXXXXXXXXXXXXX      
      write(*,*) '----------LVS BREAKUP---------------- '      
      NLVS1=NLVS+(level-1)
      ! NLVS=NLVS+(LEVEL-1)
      write(*,*) 'BREAKUP OF LVS=',LLVS       
      write(*,*) 'NLVS          =',NLVS     
      do llevel=1,level
      write(*,*) 'level and volume',llevel,vol_tmps(llevel)
      enddo
!XXXXXXXXXXXXXXXXXXXXXXXXX TAGGING DIFFERENT LEVEL XXXXXXXXXXXXXXXXXXXXXX       
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX    

        DO 888 LLEVEL=1,LEVEL !LEVEL,1,-1 !
        !INITIALISE FOR LLEVEL>1
!$OMP PARALLEL DO private(I,J)              
        DO K=-2,N3F+3
        DO J=-2,N2F+3
        DO I=-2,N1F+3
         ALPHI_EXT(I,J,K,1)=0.
         MASK_TAG(I,J,K)=0
        ENDDO
        ENDDO
        ENDDO

       WRITE(*,*) LLEVEL,':THE LEVEL WORKING ON OUT OF', LEVEL 
       
       IF (LLEVEL .EQ. 1) LVS1=LLVS
       IF (LLEVEL .GE. 2) then
         DO L=1,NLVS
         IF (LVSON(L) .NE. 1) THEN 
         LVS1=L  !       LVS1=NLVS1+(LLEVEL-1)  
         LVSON(LVS1)=1
         write(*,*) lvs1,'was empty'
         goto 504
         ENDIF
         ENDDO
         NLVS=NLVS+1
         LVS1=NLVS
         LVSON(LVS1)=1
         write(*,*) 'none was empty, lvs1=',nlvs
  504  continue
       endif
      
       !여기서 미리 VOL_TOT_TMP로 판별하면 넘 작은 애들은 BAND_GENERATION안해도 됨
      IF (VOL_TMPS(LLEVEL) .EQ. 0.) THEN       
        WRITE(*,*)'ZERO VOLUME NOT COUNT'         
        if (lvs1.eq. nlvs)  NLVS=NLVS-1
        LVSON(LVS1)=0
        GOTO 888
      ELSE IF((VOL_TMPS(LLEVEL).NE.0.).AND.
     &                              (VOL_TMPS(LLEVEL).LT.DIVOL)) THEN 
       WRITE(*,*) 'TOTAL BUBBLE VOL = ',LVS1,VOL_TMPS(LLEVEL)   
       WRITE(*,*) 'CALL LAGPRT'
       WRITE(*,*)  XXX(LLEVEL),YYY(LLEVEL),ZZZ(LLEVEL)       
       CALL LAGPRT(VOL_TMPS(LLEVEL),LVS1,
     &                          XXX(LLEVEL),YYY(LLEVEL),ZZZ(LLEVEL))
       if (lvs1 .eq. nlvs) nlvs=nlvs-1
       LVSON(LVS1)=0
       GOTO 888
       ENDIF 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!       
       NUMS_INITIAL=0                              
!$OMP PARALLEL DO private(N,NN,I,J,K) 
!$OMP&reduction(+:NUMS_INITIAL)
      DO N=1,N_MAX       
      DO NN=1,NUMA_TMP(N)
      I=I_TMP(NN,N)
      J=J_TMP(NN,N)
      K=K_TMP(NN,N)   
       IF ((ABS(ALPHI_TMP(I,J,K)) .Le. D_CRI) .AND. ! 이게 엄청 가까운 상태에서 갈라진 곳을 찾으니까
     &   ((TAG(I,J,K)   .EQ. LLEVEL) .OR.             ! 저쪽의 꼬투리가 여기까지 오는 거 같..?
     &    (TAG(I-1,J,K) .EQ. LLEVEL) .OR.            ! IF INSIDE D_CRI AND TAG=LLEVEL AROUND IT
     &    (TAG(I+1,J,K) .EQ. LLEVEL) .OR.            ! 이거 주변의 27개로 하면 안됨
     &    (TAG(I,J-1,K) .EQ. LLEVEL) .OR.
     &    (TAG(I,J+1,K) .EQ. LLEVEL) .OR.
     &    (TAG(I,J,K-1) .EQ. LLEVEL) .OR.
     &    (TAG(I,J,K+1) .EQ. LLEVEL))  
     &                              )   THEN  !YOU HAVE TO MAKE ALPHI_2XT=0. EXCEPTION OUT OF THE BANDS--> nn,n,llvs mode로 계산
       NUMS_INITIAL=NUMS_INITIAL+1  
       MASK_TAG(I,J,K)=LLEVEL     ! THIS SIMPLIFIES THE ABOVE CRITERION
       !ALPHI_3XT(I,J,K)=ALPHI_2XT(I,J,K)
       ENDIF
      ENDDO
      ENDDO
       
      ! DISTRIBUTION IS ARBITRARY AT THE MOMENT 
      NUMS_CRI=NUMS_INITIAL/(N_MAX-1)
       write(*,*) 'numsinit, numscri',nums_initial, nums_cri
      ! initialize 
      do i=1,n_max
      NUMA(i,LVS1)=0
      enddo
      
      NUMS_INITIAL=0
      N=1     
      DO K=1,N3FM
      DO J=1,N2FM
      DO I=1,N1FM
       IF (MASK_TAG(I,J,K) .EQ. LLEVEL) THEN  
       NUMS_INITIAL=NUMS_INITIAL+1           
       I_B(NUMS_INITIAL,N,LVS1)=I
       J_B(NUMS_INITIAL,N,LVS1)=J
       K_B(NUMS_INITIAL,N,LVS1)=K
       IF (NUMS_INITIAL .EQ. NUMS_CRI) THEN
       if (n .eq. n_max) goto 771
       NUMA(N,LVS1)=NUMS_CRI  ! INITIALLY COUNTED ACTIVATED CELLS PER BAND
       NUMS_INITIAL=0
       N=N+1
       ENDIF
      ENDIF
  771 continue
      ENDDO
      ENDDO
      ENDDO
      
      NUMA(N,LVS1)=NUMS_INITIAL ! THE REMAININGS
      write(*,*)'final numsinit',NUMS_INITIAL ! THE REMAININGS
  
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
        DO NN=1,NUMA(N,lvs1)
         I=I_B(NN,N,LVS1)
         J=J_B(NN,N,LVS1)
         K=K_B(NN,N,LVS1)
       IF(    (ALPHI_TMP(I,J,K)*ALPHI_TMP(IPF(I),J,K) .LT. 0.)    
     &   .OR. (ALPHI_TMP(I,J,K)*ALPHI_TMP(IMF(I),J,K) .LT. 0.)
     &   .OR. (ALPHI_TMP(I,J,K)*ALPHI_TMP(I,JPF(J),K) .LT. 0.)
     &   .OR. (ALPHI_TMP(I,J,K)*ALPHI_TMP(I,JMF(J),K) .LT. 0.)     
     &   .OR. (ALPHI_TMP(I,J,K)*ALPHI_TMP(I,J,KPF(K)) .LT. 0.)     
     &   .OR. (ALPHI_TMP(I,J,K)*ALPHI_TMP(I,J,KMF(K)) .LT. 0.)
     &   .OR. (ABS(PSIF2(I,J,K)-0.5) .LT. 0.5)   ! 이거 PSIF2로 해야지..? 그래야 LEVEL 계속 돌면서도 유지가 됨..?지금까지는 ALPHI_TMP가 잘 잡아줬나봐? 이거 해도 BAND=0뜨면 안되니까 뒤에 처리 필요함 
     &                       ) THEN !!.OR. (ABS(PSIF_MTP(I,J,K)-0.5) .LT. 0.5)) 
           NUMS=NUMS+1
           I_B(NUMS,1,LVS1)=I
           J_B(NUMS,1,LVS1)=J
           K_B(NUMS,1,LVS1)=K
           MASK_S(I,J,K)=1
           ALPHI(NN,1,LVS1)  =ALPHI_TMP(I,J,K)
           ALPHI_EXT(I,J,K,1)=ALPHI_TMP(I,J,K)
           ALPHI_3XT(I,J,K)  =ALPHI_TMP(I,J,K)
       ELSE
           MASK_S(I,J,K)=2 ! else, not first band
       ENDIF
        ENDDO
       ENDDO      
       NUMA(1,LVS1)=NUMS  ! resetting the initial band cells  
       !NUMS_CHECK=NUMA(1,LVS1)
       write(*,*)'nums of first band', nums
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!C===============================MAKE S-TUBE============================C
      !2018-01-03_ IF NUMS=0; THERE IS SOMETHING THERE, BUT LVS CANNOT CATCH IT OUT
      !           THE VOLUME HAS TO BE DETECTABLE IN ORDER TO BE TRANSFERRED TO LAGRANGIAN PARTICLE
      !           THUS, IF NUMS=0, OTHER METHODS NEEDED....
      !2018-01-13_ MAJOR ERROR WITH THE VOLUME CONSERVATION WAS DUE TO THE ABSENCE OF 
      !          VOLUME CORRECTION. THE ORIGINAL PSIF_MTP SUM IS STORED, AND USED FOR THE LAGPRTS
      !          ONLY THE ONES THAT REQUIRES LVS-FUNCT ARE RE-GENERATED HERE
      !IF (NUMS_CHECK .NE. 0) THEN !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
          DO 1003 N=2,N_MAX
          NUMS=0
          NM=N-1
          DO 1004 NN=1,NUMA(NM,LVS1)  !AT S-TUBE
          I=I_B(NN,NM,LVS1)
          J=J_B(NN,NM,LVS1)
          K=K_B(NN,NM,LVS1)  

         DO 1005 K2=K-1,K+1
            KK=K2
            IF (IPZ .EQ. 1) THEN
             IF (KK .EQ. 0) KK=N3FM
             IF (KK .EQ. N3F) KK=1  
            ELSE
             IF ((K2 .EQ. 1) .OR. (K2 .EQ. N3F)) goto 1005
            ENDIF
         DO 1006 J2=J-1,J+1
            JJ=J2
            IF (IPY .EQ. 1) THEN
             IF (JJ .EQ. 0) JJ=N2FM
             IF (JJ .EQ. N2F) JJ=1  
            ELSE
             IF ((J2 .EQ. 1) .OR. (J2 .EQ. N2F)) goto 1006
            ENDIF
          DO 1007 I2=I-1,I+1
            II=I2
            IF (IPX .EQ. 1) THEN
             IF (II .EQ. 0) II=N1FM
             IF (II .EQ. N1F) II=1  
            ELSE
             IF ((I2 .EQ. 1) .OR. (I2 .EQ. N1FM))  goto 1007
            ENDIF            
           IF ( MASK_S(II,JJ,KK) .NE. 1 ) THEN ! MAKE NEW BANDS!!
          NUMS=NUMS+1
          I_B(NUMS,N,LVS1)=II
          J_B(NUMS,N,LVS1)=JJ
          K_B(NUMS,N,LVS1)=KK
          MASK_S(II,JJ,KK)=1       
       ENDIF
 1007   CONTINUE
 1006   CONTINUE
 1005   CONTINUE
 1004   CONTINUE
           NUMA(N,LVS1)=NUMS
 1003   CONTINUE
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ICHECK=0
      
!C===============================MAKE S-TUBE============================C
        DO 1103 N=1,N_MAX
         DO 1104 NN=1,NUMA(N,LVS1)  !AT S-TUBE
          I=I_B(NN,N,LVS1)
          J=J_B(NN,N,LVS1)
          K=K_B(NN,N,LVS1)  
          
         DO 1105 K2=K-1,K+1
            KK=K2
            IF (IPZ .EQ. 1) THEN
             IF (KK .EQ. 0) KK=N3FM
             IF (KK .EQ. N3F) KK=1  
            ELSE
             IF ((K2 .EQ. 1) .OR.(K2 .EQ. N3F)) goto 1105
            ENDIF
         DO 1106 J2=J-1,J+1
            JJ=J2
            IF (IPY .EQ. 1) THEN
             IF (JJ .EQ. 0) JJ=N2FM
             IF (JJ .EQ. N2F) JJ=1  
            ELSE
             IF ((J2 .EQ. 1) .OR. (J2 .EQ. N2F))  goto 1106
            ENDIF
          DO 1107 I2=I-1,I+1
            II=I2
            IF (IPX .EQ. 1) THEN
             IF (II .EQ. 0) II=N1FM
             IF (II .EQ. N1F) II=1 
            ELSE
             IF ((I2 .EQ. 1) .OR. (I2 .EQ. N1FM)) goto 1107
            ENDIF
       IF (ICHECK(II,JJ,KK) .EQ. 0) THEN
       IF ( ((TAG(II,JJ,KK) .NE. LLEVEL) .AND. (TAG(II,JJ,KK) .NE.0))
     & .or. (n.ge.10)) THEN ! IF IT IS ON OTHER LEVEL or AT THE LAST BAND..? 
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
        ENDIF!IF ( K2 .EQ. 0) THEN
      ENDIF
       
       IF (ALPHI_3XT(I,J,K) .GE. 0.)  THEN
       ALPHI_3XT(II,JJ,KK)=ALPHI_3XT(I,J,K)+SQRT(DX**2+DY**2+DZ**2) 
       ICHECK(II,JJ,KK)=1
       ELSE                               
       ALPHI_3XT(II,JJ,KK)=ALPHI_3XT(I,J,K)-SQRT(DX**2+DY**2+DZ**2)
       ICHECK(II,JJ,KK)=1       
       ENDIF
       ! WRITE(*,*)I,J,K  
       
       ELSE
       ALPHI_3XT(II,JJ,KK)=ALPHI_TMP(II,JJ,KK)   
       ICHECK(II,JJ,KK)=1       
       ENDIF !IF( (MASK_S(II,JJ,KK).EQ. 0) .OR. (N .GE. 11) ) THEN
       ENDIF 
 1107   CONTINUE
 1106   CONTINUE
 1105   CONTINUE
 1104  CONTINUE
!C==========================makes a I_B,J_B,K_B=========================C
 1103  CONTINUE

       NUMA_MAX=0
       DO N=1,N_MAX
         NUMA_MAX=MAX(NUMA_MAX,NUMA(N,LVS1))
       ENDDO
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
          ALPHI(NN,N,LVS1)=ALPHI_3XT(I,J,K)
        ENDDO
        ENDDO

!$OMP PARALLEL DO    
        DO K=1,N3F
        DO J=1,N2F
        DO I=1,N1F
        ALPHI_EXT(I,J,K,1)=0
        ENDDO
        ENDDO
        ENDDO
        
        DO N=1,N_MAX
!$OMP PARALLEL DO PRIVATE(I,J,K)
        DO NN=1,NUMA(N,LVS1)
         I=I_B(NN,N,LVS1)
         J=J_B(NN,N,LVS1)
         K=K_B(NN,N,LVS1)                                      
        ALPHI_EXT(I,J,K,1)=ALPHI(NN,N,LVS1)
        ! COLOR(I,J,K)=LVS1
        ENDDO
        ENDDO        
        
        VOL_TOT_ORI(LVS1)=VOL_TMPS(LLEVEL)
        CALL GLOBAL_MASS_CORRECTION(NN1,LVS1,VOL_TOT_ORI,dtconst)                  
        WRITE(*,*) 'BAND_gen-TOTAL BUBBLE VOL = ',LVS1,VOL_TOT_ORI(LVS1)             

      !PSIF DOES NOT NEED OT BE CHANGED HERE
      !COLLISION ASSESSED AFTER THE SUBROUTINE USING THE PSIF
      !PSIF IS SET THERE.                                                                  
       
 888  CONTINUE ! ENDDO ! this is end of llevel=1, level


!$OMP PARALLEL DO            
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       ALPHI_EXT(I,J,K,1)=0.
      ENDDO
      ENDDO
      ENDDO          
      ! WE DO NOT NEED TO CALCULATE THE ALPHI AGAIN            
        DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)        
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
          ALPHI_EXT(I,J,K,1)=ALPHI(NN,N,LLVS) 
        ENDDO
        ENDDo
        
        
        



  777 CONTINUE
 
      RETURN
      END
 
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE LAGLVS(U,V,W,VOL_TOT_ORI)
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR      
      USE LVS_VAR
      USE LVS_GEOM_VAR
      USE LVS_COUPLING
      
      IMPLICIT NONE
      
      INTEGER*8   I,J,K,II,JJ,KK,IG,JG,KG
      INTEGER*8   I2,J2,K2
      INTEGER*8   IP1,IM1,JP1,JM1,KP1,KM1
      INTEGER*8   IL,IH,JL,JH,KL,KH,iilvs
      INTEGER*8   LPRT,NUM_CHECK,NUMC,IPRT
      INTEGER*8   ILVS(30),iTAG(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      INTEGER*8   LVS1,LVS2,N,NN,NN1,NN2
      
      REAL*8      VOL_TOT_ORI(MLVS)
      REAL*8      U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL*8      XPRT,YPRT,ZPRT
      REAL*8      RAD,PI      
      REAL*8      GX,GY,GZ,AGI,DALPHI,GXD,DIS
      
      PI=ACOS(-1D0)   
      NN1=N_MAX-3
      NN2=N_MAX-1
      
      DO 999 LPRT=1,NPRT
      NUM_CHECK=0
      XPRT=XLAGPRT(LPRT,1)
      YPRT=XLAGPRT(LPRT,2)
      ZPRT=XLAGPRT(LPRT,3)
      RAD=0.5D0*(6D0*XLAGPRT(LPRT,4)/PI)**(1D0/3D0)

!$OMP PARALLEL DO private(I,J)
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       ALPHI_EXT(I,J,K,1)=0.
       iTAG(i,j,k)=0
      ENDDO
      ENDDO
      ENDDO
 


      ! write(*,*)'laglvs', xprt,yprt,zprt
 
      !FIND WHERE THE PARTICLE IS: TO CHECK MASK_BUB  
      ! AT THE REFINED GRID, IG, JG, KG : CORRESPONDING NNS GRID, I2, J2, K2
      DO 100 II=1,N1FM-1
      IF (XF(II+1) .GT. XPRT) THEN
      IG=II
      GOTO 110
      ENDIF
  100 CONTINUE 
  110 CONTINUE   
      DO 101 JJ=1,N2FM-1
      IF (YF(JJ+1) .GT. YPRT) THEN 
      JG=JJ
      GOTO 111
      ENDIF
  101 CONTINUE 
  111 CONTINUE   
      DO 102 KK=1,N3FM-1
      IF (ZF(KK+1) .GT. ZPRT) THEN
      KG=KK
      GOTO 112
      ENDIF
  102 CONTINUE 
  112 CONTINUE
          
      ! write(*,*)'index',ig,jg,kg
      
      I2=ICOUMP_VEL(IG)
      J2=JCOUMP_VEL(JG)
      K2=KCOUMP_VEL(KG)      
      !FIGURE OUT WHAT LVS IS CLOSE TO THE PARTICLE
      DO K=K2-1,K2+1
      kk=k
       if (ipz .eq. 1) then 
        if (kk .eq. 0 ) kk = n3m
        if (kk .eq. n3) kk = 1
       endif
      DO J=J2-1,J2+1 
      jj=j
        if ( jj .eq. 0) jj = 1
        if ( jj .eq. n2) jj=n2m
      DO I=I2-1,I2+1
      ii=i
        if (ii .eq. 0) ii=1
        if (ii .eq. n1) ii=n1m
       
       IF (MASK_BUB(iI,jJ,kK) .NE. 0) THEN ! EVERYTING STARTS IF THERE IS SOMETHING THERE.
       NUM_CHECK=NUM_CHECK+1      
       ILVS(NUM_CHECK)=MASK_BUB(iI,jJ,kK)    ! NUMCHECK[0,27]   
        ! write(*,*) mask_bub(i,j,k)
       ENDIF
      ENDDO
      ENDDO
      ENDDO
      
      
      ! IF THERE IS NO LVS CLOSE, NEXT PARTICLE
      IF (NUM_CHECK .EQ. 0) GOTO 999
      ! IF NOT, MAKE AN ALPHI FIELD WITH LVS TAGGED TO THE CELL.--------
      IF (NUM_CHECK .GT. 1) THEN     
      DO NUMC=1,NUM_CHECK
      iilvs=ilvs(num_check)
      if ( lvson(iilvs) .eq. 1) then 
      LVS2=ILVS(NUMC)
      else
      goto 999
      endif
      ! write(*,*) lvs2
!$OMP PARALLEL DO private(NN,I,J,K)       
       DO N=1,NN2
       DO NN=1,NUMA(N,LVS2)
       I=I_B(NN,N,LVS2)
       J=J_B(NN,N,LVS2)
       K=K_B(NN,N,LVS2)
        IF (iTAG(I,J,K) .EQ. 0) THEN 
        ALPHI_EXT(I,J,K,1)=ALPHI(NN,N,LVS2)
        iTAG(I,J,K)=LVS2
        ELSE                                                        
        ALPHI_EXT(I,J,K,1)=DMIN1(DABS(ALPHI_EXT(I,J,K,1)),            ! THERE CAN BE AN EXCEPTION
     &                           DABS(ALPHI(NN,N,LVS2)))             ! BUT IT IS A MATTER OF ~SDXF
        IF (ALPHI_EXT(I,J,K,1) .EQ. ALPHI(NN,N,LVS2)) iTAG(I,J,K)=LVS2! IT CAN BE ADDED TO EITHER. 
        ENDIF !TAG-RELATED
       ENDDO
       ENDDO
      ENDDO    
      ELSE
      iTAG(IG,JG,KG)=ILVS(1)      
      ENDIF
      
      ! THE LVS OF INTEREST IS THEREFORE...
      if (itag(ig,jg,kg) .eq. 0) then 
      do k=kg-3,kg+3
      do j=jg-3,jg+3
      do i=ig-3,ig+3
      LVS1=iTAG(I,J,K) ! THE LVS WITH THE MIN ALPHI AT THE PARTICLE CONTAINING CELL CENTER
      enddo
      enddo
      enddo
      else
      lvs1=itag(ig,jg,kg)
      endif
      
      if (lvs1 .eq. 0) goto 999
       ! write(*,*)'is lvs1 problem?',lvs1
!$OMP PARALLEL DO private(I,J)
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       ALPHI_EXT(I,J,K,1)=0.
      ENDDO
      ENDDO
      ENDDO
!$OMP PARALLEL DO private(NN,I,J,K)      
      DO N=1,NN2
      DO NN=1,NUMA(N,LVS1)
      I=I_B(NN,N,LVS1)
      J=J_B(NN,N,LVS1)
      K=K_B(NN,N,LVS1)
       ALPHI_EXT(I,J,K,1)=ALPHI(NN,N,LVS1)
      ENDDO
      ENDDO
      
      
      
      
      ! PARTICLE (X-X_CM).(N)=D(ALPHI_CM)
        IP1=IPF(IG)
        IM1=IMF(IG)
        JP1=JPF(JG)
        JM1=JMF(JG)
        KP1=KPF(KG)
        KM1=KMF(KG)   
      GX=0.5D0*( ALPHI_EXT(IP1,JG,KG,1)-ALPHI_EXT(IM1,JG,KG,1) )*SSDXF
      GY=0.5D0*( ALPHI_EXT(IG,JP1,KG,1)-ALPHI_EXT(IG,JM1,KG,1) )*SSDYF
      GZ=0.5D0*( ALPHI_EXT(IG,JG,KP1,1)-ALPHI_EXT(IG,JG,KM1,1) )*SSDZF     
      IF (GX .EQ. 0. .AND. GY .EQ. 0. .AND. GZ .EQ. 0. ) THEN 
      GOTO 999 ! NOT ANYWHERE NEAR ALPHI OAO
      ELSE
      AGI=1D0/DSQRT(GX**2D0 + GY**2D0 + GZ**2D0)     
      ENDIF
      DALPHI=AGI*((XPRT-XPF(IG))*GX+(YPRT-YPF(JG))*GY+(ZPRT-ZPF(KG))*GZ)
      GXD=ALPHI_EXT(IG,JG,KG,1)+DALPHI
      !ASSESS WHETER THE PARTICLE IS CLOSE ENOUGH, AND IF IT IS, ADD TO LVS1, ERASE FROM PRT.
      IF (GXD .LE. RAD) THEN 

      ! THIS IS FOR ADDING ALPHI_EXT OF THE PARTICLE
      IL=ICOU1(I2-1) ! LEVEL SET INDICES OF THE ASSESSED REGION
      IH=ICOU2(I2+1)      
      JL=JCOU1(J2-1)
      JH=JCOU2(J2+1)      
      KL=KCOU1(K2-1)
      KH=KCOU2(K2+1)
      
      DO K=KL,KH
      DO J=JL,JH
      DO I=IL,IH    
      DIS=DSQRT((XPRT-XPF(I))**2+(YPRT-YPF(J))**2+(ZPRT-ZPF(K))**2)
      ALPHI_EXT(I,J,K,1)= DMIN1(ALPHI_EXT(I,J,K,1),DIS-RAD) ! SHOULD BE CHANGED FOR DROPLET
      ENDDO
      ENDDO
      ENDDO
      
      VOL_TOT_ORI(LVS1)=VOL_TOT_ORI(LVS1)+XLAGPRT(LPRT,4)
      XLAGPRT(LPRT,4)=0.
      CALL GLOBAL_MASS_CORRECTION(NN1,LVS1,VOL_TOT_ORI,dtconst)
      WRITE(*,*) 'LAG TO LVS(N): VOL = ',LVS1,VOL_TOT_ORI(LVS1)

      
      DO K=K2-1,K2+1
      DO J=J2-1,J2+1
      DO I=I2-1,I2+1
      DIS=DSQRT((XPRT-X(I))**2+(YPRT-YP(J))**2+(ZPRT-ZP(K))**2)    
      IF (DIS .LE. RAD) U(I,J,K)=VLAGPRT(LPRT,1)
      DIS=DSQRT((XPRT-XP(I))**2+(YPRT-Y(J))**2+(ZPRT-ZP(K))**2)     
      IF (DIS .LE. RAD) V(I,J,K)=VLAGPRT(LPRT,2)
      DIS=DSQRT((XPRT-XP(I))**2+(YPRT-YP(J))**2+(ZPRT-Z(K))**2)      
      IF (DIS .LE. RAD) W(I,J,K)=VLAGPRT(LPRT,3)      
      ENDDO
      ENDDO
      ENDDO

      ! save to the lvs...
!$OMP PARALLEL DO private(NN,I,J,K)      
      DO N=1,NN2
      DO NN=1,NUMA(N,LVS1)
      I=I_B(NN,N,LVS1)
      J=J_B(NN,N,LVS1)
      K=K_B(NN,N,LVS1)
       ALPHI(NN,N,LVS1)=ALPHI_EXT(I,J,K,1)
      ENDDO
      ENDDO      
      ENDIF 
      ! IF NOT CLOSE ENOUGHT, JUST END. 
      
      
  999 CONTINUE    
      return 
      end
      

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE THIN_CRI(UF,VF,WF,LLVS)
!일단은 transport uf,vf,wf deallocate 전에 criterion 을 계산해야해서, 
! transport 마지막에 넣었는데, 이건 다른 lvs 사이의 resolution은 안됨
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      USE PARAM_VAR
      USE FLOW_GEOM_VAR
      USE FLOW_VAR      
      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE
!``````````````````````````````````````````````````````````````````````      
      REAL*8    UF(M1F,M2F,M3F),VF(M1F,M2F,M3F),WF(M1F,M2F,M3F)
      REAL*8    GX,GY,GZ,AGI,UC,VC,WC,DIS
      REAL*8    UNORM(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      REAL*8    NORM(-2:M1F+3,-2:M2F+3,-2:M3F+3,3)
      REAL*8    XTAG(-2:M1F+3,-2:M2F+3,-2:M3F+3,3)
      REAL*8    PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)      
      INTEGER*8 LLVS
      INTEGER*8 I,J,K
      INTEGER*8 IP1,JP1,KP1
      INTEGER*8 IM1,JM1,KM1
      INTEGER*8 N,NN,NN1,NN2
      INTEGER*8  TAG(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      INTEGER*8 TTAG(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      INTEGER*8 NCOUNT
!``````````````````````````````````````````````````````````````````````
      REAL*8    X_TEST,Y_TEST,Z_TEST,CR
      INTEGER*8 ITEST ,IGEO ,JGEO ,KGEO
      INTEGER*8 NPTEST,IGEO2,JGEO2,KGEO2
      INTEGER*8 II,JJ,KK,II2,JJ2,KK2,IT
      
      
      
      NN1=N_MAX-3
      NN2=N_MAX-1
      
!$OMP PARALLEL DO      
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
      TAG(I,J,K)=0
      TTAG(I,J,K)=0
      
      ENDDO
      ENDDO
      ENDDO
      
!$OMP PARALLEL DO private(I,J)
       DO K=1,N3FM
       DO J=1,N2FM
       DO I=1,N1FM
         PSIF_MTP(I,J,K)=0.        
       ENDDO
       ENDDO
       ENDDO
      ! 인간은 같은 실수를 반복한다 ㅜㅠ
      call cal_psif_mtp(1,NN1,LLVS,PSIF_MTP)


!$OMP PARALLEL DO private(N,NN,I,J,K,IP1,IM1,JP1,JM1,KP1,KM1)
!$OMP&private(GX,GY,GZ,AGI,UC,VC,WC)     
      DO N=1,NN1-1!3~4 까지 해야할 것 같은데 혹시 모르니까 전부다 해주자...
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
        ! FOR THE INTERFACIAL CELLS, DIRECTLY ADJACENT TO ALPHI=0, ONLY.
        ! ALPHI <0 ONLY: NORMAL VECTOR TOWARDS THE CELL CENTER IS INWARD TO THE INTERFACE
       IF (   (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(IP1,J,K,1) .LT. 0D0)    
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(IM1,J,K,1) .LT. 0D0)
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,JP1,K,1) .LT. 0D0)
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,JM1,K,1) .LT. 0D0)     
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,J,KP1,1) .LT. 0D0)     
     &   .OR. (ALPHI_EXT(I,J,K,1)*ALPHI_EXT(I,J,KM1,1) .LT. 0D0)
     &                                                         ) THEN  
      IF (ALPHI_EXT(I,J,K,1) .LT. 0D0) TAG(I,J,K)=1 !SHOULD CHANGE TO .GT. FOR DROPLETS!! 
      GX=0.5D0*( ALPHI_EXT(IP1,J,K,1)-ALPHI_EXT(IM1,J,K,1) )*SSDXF
      GY=0.5D0*( ALPHI_EXT(I,JP1,K,1)-ALPHI_EXT(I,JM1,K,1) )*SSDYF
      GZ=0.5D0*( ALPHI_EXT(I,J,KP1,1)-ALPHI_EXT(I,J,KM1,1) )*SSDZF     
      AGI=1D0/DSQRT(GX**2D0 + GY**2D0 + GZ**2D0)
      
      NORM(I,J,K,1)=GX*AGI
      NORM(I,J,K,2)=GY*AGI
      NORM(I,J,K,3)=GZ*AGI
      
      UC=0.5D0*( UF(I,J,K)+UF(IP1,J,K) )
      VC=0.5D0*( VF(I,J,K)+VF(I,JP1,K) )
      WC=0.5D0*( WF(I,J,K)+WF(I,J,KP1) )
      
      UNORM(I,J,K)=( UC*GX + VC*GY + WC*GZ )*AGI
      
      ENDIF
      ENDDO 
      ENDDO
      
      DIS=0.25D0*DMIN1(SDXF,SDYF,SDZF)
      
! !$OMP PARALLEL DO 
! !$OMP&private(X_TEST,Y_TEST,Z_TEST,ITEST,I,J,NPTEST)
! !$OMP&private(II,JJ,KK,II2,JJ2,KK2,IT)
! !$OMP&private(IGEO,JGEO,KGEO,IGEO2,JGEO2,KGEO2)
! !$OMP&private(CR)
      DO 99 K=1,N3FM
      DO 99 J=1,N2FM
      DO 99 I=1,N1FM
      NPTEST=0
      
      ! FOR AN INTERFACE CELL
      IF (TAG(I,J,K) .EQ. 1) THEN       
      DO ITEST=1,200 ! ? 
      
      X_TEST=XPF(I)-ITEST*DIS*NORM(I,J,K,1) ! MINUS, BECAUSE THE TAGGING SHOULD GO INWARD (NORM > 0) . + IF DROPLET
      Y_TEST=YPF(J)-ITEST*DIS*NORM(I,J,K,2)
      Z_TEST=ZPF(K)-ITEST*DIS*NORM(I,J,K,3)  
            IF ((IPZ .EQ. 1) .AND. (Z_TEST .LT. ZL  )) Z_TEST=Z_TEST-ZL
      DO 100 II=1,N1FM-1
      IF (XF(II+1) .GT. X_TEST) THEN
      IGEO=II
      GOTO 110
      ENDIF
  100 CONTINUE 
  110 CONTINUE        
      DO 101 JJ=1,N2FM-1
      IF (YF(JJ+1) .GT. Y_TEST) THEN 
      JGEO=JJ
      GOTO 111
      ENDIF
  101 CONTINUE 
  111 CONTINUE          
      DO 102 KK=1,N3FM-1
      IF (ZF(KK+1) .GT. Z_TEST) THEN
      KGEO=KK
      GOTO 112
      ENDIF    
  102 CONTINUE 
  112 CONTINUE 

      !AS SOON AS ALPHI>0; NPTEST=ITEST. IF ITEST> 200?********************  
      IF (psif_mtp(IGEO,JGEO,KGEO) .gt. 0.5)  THEN 
      NPTEST=ITEST
      ! CALCULATE THE CONTRACTION AT THE START AND FINISHING CELL  
      CR=UNORM(I,J,K)+UNORM(IGEO,JGEO,KGEO)
      IF ((NPTEST .LT. 30) .AND. (CR .LE. 0D0) )   THEN 
      write(*,600)i,j,k
      write(*,601)NORM(I,J,K,1),NORM(I,J,K,2),NORM(I,J,K,3)
  600 format(3I5) 
  601 format(3F15.8) 
            TTAG(I,J,K)=1 
      DO IT=1,NPTEST 
        X_TEST=XPF(I)-IT*DIS*NORM(I,J,K,1) 
        Y_TEST=YPF(J)-IT*DIS*NORM(I,J,K,2)
        Z_TEST=ZPF(K)-IT*DIS*NORM(I,J,K,3)
        IF ((IPZ .EQ. 1) .AND. (Z_TEST .LT. ZL  )) Z_TEST=Z_TEST-ZL
            DO 200 II2=1,N1FM-1
            IF (XF(II2+1) .GT. X_TEST) THEN
            IGEO2=II2          
            GOTO 210           
            ENDIF              
  200       CONTINUE           
  210       CONTINUE                                    
            DO 201 JJ2=1,N2FM-1 
            IF (YF(JJ2+1) .GT. Y_TEST) THEN
            JGEO2=JJ2          
            GOTO 211           
            ENDIF              
  201       CONTINUE           
  211       CONTINUE                                    
            DO 202 KK2=1,N3FM-1 
            IF (ZF(KK2+1) .GT. Z_TEST) THEN
            KGEO2=KK2
            GOTO 212
            ENDIF      
  202       CONTINUE 
  212       CONTINUE   
        TTAG(IGEO2,JGEO2,KGEO2)=1        
      ENDDO!IT=1,NPTEST  
      ! AFTER TAGGING, DON'T LINGER IN THE ITEST=1,200, AND GO OUT TO SEEK OTHER CELL
      GOTO 99
      ELSE     ! IF ( NPTEST > 25 )
      GOTO 99  ! NOTHING HAPPENS, LOOK FOR THE NEXT CELL
      ENDIF    ! IF (NPTEST .LT. 25) : 0.25*SDXF 기준으로 해서 NPTEST < 24 이면 6dXG 정도       
      ENDIF
      ENDDO! ITEST=1,200   
  
      ENDIF    ! IF (TAG(I,J,K) .EQ. 1) THEN 
  99  CONTINUE ! GET ANOTHER INTERFACIAL CELL
      
      OPEN(61,FILE='checkTAG.dat')
      WRITE(61,*) 'VARIABLES="X","Y","Z","psif","TTAG"'
      WRITE(61,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
      DO K=1,N3FM
      DO J=1,N2FM
      DO I=1,N1FM
        WRITE(61,149) XPF(I),YPF(J),ZPF(K),psif_mtp(i,j,k),TTAG(i,j,k)
      ENDDO
      ENDDO
      ENDDO
       CLOSE(61)
 149  FORMAT(4F15.8,1I3)   
      stop
         
      
      RETURN 
      END
      

!

!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      SUBROUTINE LVSALLTOONE
!1LVS 으로 하려는 시도 TEST 위해 여러개의 LVS 를 하나로 묶어본다.
!XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
      USE FLOW_VAR

      USE PARAM_VAR
      USE FLOW_GEOM_VAR

      USE LVS_VAR
      USE LVS_GEOM_VAR
      
      IMPLICIT NONE

      ! REAL U(0:M1,0:M2,0:M3),V(0:M1,0:M2,0:M3),W(0:M1,0:M2,0:M3)
      REAL PSI_CN(0:M1,0:M2,0:M3)

      REAL PSIF_MTP(M1L:M1U,M2L:M2U,M3L:M3U)
      REAL VOL_TOT_ORI(MLVS)

      INTEGER BUBCOL(MLVS)
      INTEGER BUBCOL2(MLVS)

       INTEGER MASK_S(-2:M1F+3,-2:M2F+3,-2:M3F+3)     
      INTEGER I,J,K,N,NN,LLVS,NN1,NN2
      INTEGER II,JJ,KK,I2,J2,K2,NM
      INTEGER LVS_START,LVS_END,NUMS,NUMS_INITIAL,NUMS_CRI
      REAL XX,YY,ZZ,ALPHI_TMP,DVOL,CRI
      INTEGER ISKIP,JSKIP,KSKIP
      REAL PI
      REAL ASSIGN_LVS
      REAL ALPHI_EXT2(-2:M1F+3,-2:M2F+3,-2:M3F+3)
      REAL ALPHI_EXT3(-2:M1F+3,-2:M2F+3,-2:M3F+3)
       REAL DX,DY,DZ     
      !ALLOCATE( ALPHI_EXT(-2:M1F+3,-2:M2F+3,-2:M3F+3,1))

      

       NN1=N_MAX-3
       NN2=N_MAX-1
       
       CRI=1.2*MAX(SDXF,SDYF,SDZF) 
!$OMP PARALLEL DO private(I,J,K)
      DO K=-2,N3F+3
      DO J=-2,N2F+3
      DO I=-2,N1F+3
       ALPHI_EXT(I,J,K,1)=0D0
       ALPHI_EXT2(I,J,K)= 0D0
       ALPHI_EXT3(I,J,K)= 0D0
      ENDDO
      ENDDO
      ENDDO
      
      
      ! save all the level-sets into one lvs
      DO 776 LLVS=1,NLVS !MULTIPLE_LEVEL_SET_APPROACH

       DO N=1,N_max
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
         I=I_B(NN,N,LLVS)
         J=J_B(NN,N,LLVS)
         K=K_B(NN,N,LLVS)
         ALPHI_EXT(I,J,K,1)=ALPHI(NN,N,LLVS)
        ENDDO
       ENDDO
       !CALL ALPHI_BC(LLVS)
       
  776 CONTINUE

         
  
      !THIS IS FOR EQUALLY DISTRIBUTE THE INITIAL LEVEL-SET TO EACH NUMBER OF BAND.
      NUMS_INITIAL=0
!$OMP PARALLEL DO private(I,J,K,ALPHI_TMP)
!$OMP&reduction(+:NUMS_INITIAL)
      DO 1000 K=1,N3FM
      DO 1000 J=1,N2FM
      DO 1000 I=1,N1FM 
        ALPHI_TMP=ALPHI_EXT(I,J,K,1)
        IF( ABS(ALPHI_TMP) .LE. CRI ) NUMS_INITIAL=NUMS_INITIAL+1
 1000  CONTINUE
          NUMS_CRI=NUMS_INITIAL/(N_MAX-1)  !-1 IS FOR SAFTY.
      write(*,*) 'nums_cri=', nums_cri
      NUMS_INITIAL=0
      N=1
C!$OMP PARALLEL DO private(I,J,K,XX,YY,ZZ,ALPHI_TMP)
C!$OMP&reduction(+:NUMS_INITIAL)
      DO 1001 K=1,N3FM
      DO 1001 J=1,N2FM
      DO 1001 I=1,N1FM 
        XX=XPF(I)
        YY=YPF(J)
        ZZ=ZPF(K)
        ALPHI_TMP=ALPHI_EXT(I,J,K,1)
        IF( ABS(ALPHI_TMP) .LE. CRI ) THEN
          NUMS_INITIAL=NUMS_INITIAL+1
          I_B(NUMS_INITIAL,N,1)=I
          J_B(NUMS_INITIAL,N,1)=J
          K_B(NUMS_INITIAL,N,1)=K
          ALPHI_EXT(I,J,K,1)=ALPHI_TMP
          IF (NUMS_INITIAL .EQ. NUMS_CRI) THEN
            NUMA(N,1)=NUMS_CRI
            NUMS_INITIAL=0
            N=N+1
          ENDIF
        ENDIF
        alphi_ext2(i,j,k)=alphi_ext(i,j,k,1)
 1001  CONTINUE
          NUMA(N,1)=NUMS_INITIAL !REMAIN
!CCCCCC------------------------FIRST-STEP LVSINIT---------------------CCC

!CCCCCC----------------------------BAND_INIT--------------------------CCC
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
        DO NN=1,NUMA(N,1)
            I=I_B(NN,N,1)
            J=J_B(NN,N,1)
            K=K_B(NN,N,1)
       IF(    (ALPHI_EXT2(I,J,K)*ALPHI_EXT2(IPF(I),J,K) .LT. 0.)    
     &   .OR. (ALPHI_EXT2(I,J,K)*ALPHI_EXT2(IMF(I),J,K) .LT. 0.)
     &   .OR. (ALPHI_EXT2(I,J,K)*ALPHI_EXT2(I,JPF(J),K) .LT. 0.)
     &   .OR. (ALPHI_EXT2(I,J,K)*ALPHI_EXT2(I,JMF(J),K) .LT. 0.)     
     &   .OR. (ALPHI_EXT2(I,J,K)*ALPHI_EXT2(I,J,KPF(K)) .LT. 0.)     
     &   .OR. (ALPHI_EXT2(I,J,K)*ALPHI_EXT2(I,J,KMF(K)) .LT. 0.)
     &                                                            ) THEN
           NUMS=NUMS+1
           I_B(NUMS,1,1)=I
           J_B(NUMS,1,1)=J
           K_B(NUMS,1,1)=K
           MASK_S(I,J,K)=1

       ELSE
           MASK_S(I,J,K)=2
       ENDIF

       ENDDO
      ENDDO
           NUMA(1,1)=NUMS
!C===============================MAKE S-TUBE============================C

        DO 1003 N=2,N_MAX
          NUMS=0
!C==========================makes a I_B,J_B,K_B=========================C
!c!$OMP PARALLEL DO private(I,J,K)  !MASK_S IS SHOULD BE SYNCRONIZED EVERY STEP.. HIGH COST. ORDER IS IMPORTANT..
          NM=N-1
         DO 1004 NN=1,NUMA(NM,1)  !AT S-TUBE
                  I=I_B(NN,NM,1)
                  J=J_B(NN,NM,1)
                  K=K_B(NN,NM,1)  

         DO 1005 K2=K-1,K+1
            KK=K2
            IF (IPZ .EQ. 1) THEN
             IF (KK .EQ. 0) KK=N3FM
             IF (KK .EQ. N3F) KK=1  
            ENDIF
         DO 1006 J2=J-1,J+1
            JJ=J2
            IF (IPY .EQ. 1) THEN
             IF (JJ .EQ. 0) JJ=N2FM
             IF (JJ .EQ. N2F) JJ=1  
            ENDIF
          DO 1007 I2=I-1,I+1
            II=I2
            IF (IPX .EQ. 1) THEN
             IF (II .EQ. 0) II=N1FM
             IF (II .EQ. N1F) II=1  
            ENDIF

!       IF ( MASK_S(II,JJ,KK) .EQ. 1 ) ->JUST GO THROUGH
       IF ( MASK_S(II,JJ,KK) .NE. 1 ) THEN
          NUMS=NUMS+1
          I_B(NUMS,N,1)=II
          J_B(NUMS,N,1)=JJ
          K_B(NUMS,N,1)=KK
          MASK_S(II,JJ,KK)=1  
      ENDIF !IF ( MASK_S(II,JJ,KK) .NE. 1 )

 1007   CONTINUE
 1006   CONTINUE
 1005   CONTINUE

 1004  CONTINUE
           NUMA(N,1)=NUMS
!C==========================makes a I_B,J_B,K_B=========================C
 1003  CONTINUE


      LVSON=0     
      LLVS=1
      LVSON(1)=1
      nlvs=1     
      
      DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,1)
            I=I_B(NN,N,1)
            J=J_B(NN,N,1)
            K=K_B(NN,N,1)
         ALPHI_EXT(I,J,K,1)=alphi_ext2(i,j,k)
         alphi_ext3(i,j,k)=alphi_ext(i,j,k,1)
        ENDDO
       ENDDO
!---------------------------------------------------------------------------      

      
        CALL CAL_PSIF_MTP(1,N_MAX,1,PSIF_MTP)
        
      DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,1)
         I=I_B(NN,N,1)
         J=J_B(NN,N,1)
         K=K_B(NN,N,1)
         PSIF(I,J,K)=PSIF_MTP(I,J,K)
        ENDDO
       ENDDO         

       OPEN(146,FILE='N#alphi_contour_AFTER_BANDG.DAT')
       WRITE(146,*) 'VARIABLES="X","Y","Z","alphi","PSIF"'
      WRITE(146,*) 'ZONE I=',N1FM,',J=',N2FM,',K=',N3FM,',F=POINT'
      DO K=1,N3FM
      DO J=1,N2FM
      DO I=1,N1FM
      WRITE(146,147) XPF(I),YPF(J),ZPF(K),
     &               alphi_ext3(i,j,k),PSIF_MTP(I,J,K)
      ENDDO
      ENDDO
      ENDDO          
       CLOSE(146)
 147  FORMAT(5F15.8)
!---------------------------------------------------------------------------

       NUMA_MAX=0
       DO N=1,N_MAX
         NUMA_MAX=max(NUMA_MAX,NUMA(N,1))
       ENDDO
        WRITE(*,*) 'NUMA_MAX/MF_BAND=',FLOAT(NUMA_MAX)/FLOAT(MF_BAND) ! MF_BAND 는 총 그리드 수의 어떤 PORTION이 되도록 설정
       IF (NUMA_MAX .GE. MF_BAND) THEN                                ! 활성화 된 CELL 의 개수가 위의 PORTION보다 크면 문제. 
        WRITE(*,*) 'MF_BAND IS NOT ENOUGH!!'
        WRITE(*,*) 'SIMULATION IS STOPPED!!'
        STOP
       ENDIF
!CCCCCC----------------------------BAND_INIT--------------------------CCC

      
       DO N=1,N_MAX
!$OMP PARALLEL DO private(I,J,K)
        DO NN=1,NUMA(N,LLVS)
            I=I_B(NN,N,LLVS)
            J=J_B(NN,N,LLVS)
            K=K_B(NN,N,LLVS)
            ALPHI(NN,N,LLVS)=ALPHI_EXT(I,J,K,1)
        ENDDO
       ENDDO
        stop
      RETURN
      END