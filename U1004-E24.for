      SUBROUTINE UEXTERNALDB(LOP,LRESTART,TIME,DTIME,KSTEP,KINC)
C     line 6 13~17 107 188~198 302 303 need to be modified for different models.
      INCLUDE 'ABA_PARAM.INC'
C
      INTEGER :: IM1,IM2,I1,I2
      PARAMETER (IM1 = 720,IM2 = 803)! IM1(IM2) is the number of the elements(nodes) in the SBPML domain.
      DIMENSION :: AInf_NOC_Na(IM1,2),AInf_NOC_Nb(IM1,2),
     1 Aele_first(IM1,2),AXofN_pm(IM2,2),Amaterials(IM1,1)
      common /matrix/ AInf_NOC_Na,AInf_NOC_Nb,
     1 Aele_first,AXofN_pm,Amaterials
      
      IF(lop.eq.0) THEN!The path of the SBPML elements' information created by the programe 'Generate matching layer parameters with one click'.
          OPEN(101,file='C:\ABAQUS\demo1\read\Inf_NOC_Na.dat')
          OPEN(102,file='C:\ABAQUS\demo1\read\Inf_NOC_Nb.dat')
          OPEN(103,file='C:\ABAQUS\demo1\read\ele_first.dat')
          OPEN(104,file='C:\ABAQUS\demo1\read\XofN_pm.dat')
          OPEN(105,file='C:\ABAQUS\demo1\read\materials.dat')
          DO I1=1,IM1
              read(101,*),AInf_NOC_Na(I1,:)
              read(102,*),AInf_NOC_Nb(I1,:)
              read(103,*),Aele_first(I1,:)
              read(105,*),Amaterials(I1,:)
          ENDDO
          DO I2 = 1,IM2
              read(104,*),AXofN_pm(I2,:)
          ENDDO
          Close(101)
          Close(102)
          Close(103)
          Close(104)
          Close(105)
      ENDIF
      
      
      RETURN
      END
C
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,
     3     NPREDF,LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,
     4     PERIOD)
C     CPE4与PML混合
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION RHS(NDOFEL,1),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(NPROPS),COORDS(MCRD,NNODE),
     2     U(NDOFEL),DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),
     4     DDLMAG(MDLOAD,*),PREDEF(2,NPREDF,NNODE),LFLAGS(*),
     5     JPROPS(*),SRESID(NDOFEL),AK(NDOFEL,NDOFEL),AM(NDOFEL,NDOFEL),
     6     CC(NDOFEL,NDOFEL)
      IF (NSVARS.eq.1240.D0) THEN
          !此单元为PML
          IF(KSTEP.EQ.1.AND.KINC.EQ.1) THEN
              call MutexInit( 1 )      ! initialize Mutex #1
              call MutexLock( 1 )      ! lock Mutex #1
              CALL FORMMATRX(AK,AM,CC,NDOFEL,NNODE,NSVARS,MCRD,
     1    MLVARX,MDLOAD,NPREDF,JELEM,COORDS,PROPS)
              call MutexUnlock( 1 )   ! unlock Mutex #1
              CALL STOREMATRICES(SVARS,NSVARS,AK,AM,CC,NDOFEL)
          ELSE
              CALL READMATRICES(SVARS,NSVARS,AK,AM,CC,NDOFEL)
          ENDIF
      ELSEIF (NSVARS.eq.208.D0) THEN
          !此单元为CPE4
          IF(KSTEP.EQ.1.AND.KINC.EQ.1) THEN
              call MutexInit( 1 )      ! initialize Mutex #1
              call MutexLock( 1 )      ! lock Mutex #1
              CALL FORMMATRXCPE4(AK,AM,CC,NDOFEL,NNODE,NSVARS,MCRD,
     1    MLVARX,MDLOAD,NPREDF,JELEM,COORDS,PROPS)
              call MutexUnlock( 1 )   ! unlock Mutex #1
              CALL STOREMATRICESCPE4(SVARS,NSVARS,AK,AM,CC,NDOFEL)
          ELSE
              CALL READMATRICESCPE4(SVARS,NSVARS,AK,AM,CC,NDOFEL)
          ENDIF
      ENDIF
      CALL OUTPUTVARIABLE(RHS,AMATRX,SVARS,PROPS,ENERGY,U,V,A,
     1 LFLAGS,DTIME,NDOFEL,NRHS,NSVARS,MLVARX,JELEM,PARAMS,AK,AM,CC)
      RETURN
      END
C
      SUBROUTINE FORMMATRX(AK,AM,CC,NDOFEL,NNODE,NSVARS,MCRD,
     1 MLVARX,MDLOAD,NPREDF,JELEM,COORDS,PROPS)
C
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION RHS(20,1),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(6),COORDS(MCRD,NNODE),
     2     U(NDOFEL),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),
     4     PREDEF(2,NPREDF,NNODE),
     5     SRESID(20),AK(NDOFEL,NDOFEL),AM(NDOFEL,NDOFEL),
     6     CC(NDOFEL,NDOFEL),dNdrs(NNODE,2),
     7     TaJacob(2,2),dNdxy(NNODE,2),Bm(3,8),
     8     Kmm(8,3),Kmmm(8,8),Kmmmm(8,8),GT(NDOFEL),AN(2,8),
     9     TN(8,2),AMM(NDOFEL,NDOFEL)
      DIMENSION D(3,3),TD(3,3),Coo(2,4),Gauss(2,4),AMu(8,8),ACu(8,8),
     1 AKu(8,8),ACuS(8,12),AKuS(8,12),AKs(12,12),AMs(12,12),ACs(12,12),
     2 Zer2(2,1),AJacob(2,2),g(2,4),Bbb(2,4),ShapeFunM(1,4),
     3 Bb1(2,8),Bb2(2,8),ANN1(1,8),ANN2(1,8),AN12(2,8),
     4 ANN3(3,12),XB(2,1),YB(2,1),ANeta(1,2),X0(2,1),Y0(2,1),
     5 Adiff_eta(1,2),Aj11(1,1),Aj12(1,1),Aj21(1,1),Aj22(1,1),
     6 A_eta(2,2),B_eta(2,2),TEMP(1,1),Ea(2,2),Eb(2,2),EYE(2,2),
     7 Fa(2,2),Fb(2,2),A1(1,2),A2(1,2),B1(1,2),B2(1,2),TEMP_CuS1(8,1),
     8 TEMP_CuS2(8,1),TEMP_CuS3(8,2),AL1(3,2),AL2(3,2),TEMP_CuS4(8,1),
     9 TEMP_CuS5(8,1),TEMP_CuS6(8,2),T1(20,20),T2(20,20)
      PARAMETER ( ZERO = 0.D0, HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0,
     1 IM1 = 720,IM2 = 803)! IM1(IM2) is the number of the elements(nodes) in the SBPML domain.
      DIMENSION TEMP_KuS1(8,1),TEMP_KuS2(8,1),TEMP_KuS3(8,2),
     1 TEMP_KuS4(8,1),TEMP_KuS5(8,1),TEMP_KuS6(8,2),AM1(20,20),
     2 AK1(20,20),AC1(20,20),AM2(20,20),AK2(20,20),AC2(20,20),
     3 AM3(20,20),AK3(20,20),AC3(20,20),T3(20,20),T4(20,20),
     4 ACSu(12,8),AKSu(12,8),AN1212(8,8),
     5 AInf_NOC_Na(IM1,2),AInf_NOC_Nb(IM1,2),
     6 Aele_first(IM1,2),AXofN_pm(IM2,2),
     7 Amaterials(IM1,1)
      integer I3_index(4,1),I3,I31
      common /matrix/ AInf_NOC_Na,AInf_NOC_Nb,
     1 Aele_first,AXofN_pm,Amaterials
C
      
      DO J1 = 1,20
          DO J2=1,20
              AM1(J1,J2) = ZERO
              AK1(J1,J2) = ZERO
              AC1(J1,J2) = ZERO
          ENDDO
      ENDDO
C     文件读取
      Ans = 1.D0
      ALpml = 1.D0
      EYE(1,1) = ONE
      EYE(1,2) = ZERO
      EYE(2,1) = ZERO
      EYE(2,2) = ONE

      DO I4=1,8
          DO N4=1,8
              AMu(I4,N4) = ZERO
              ACu(I4,N4) = ZERO
              AKu(I4,N4) = ZERO
          ENDDO
C
          !ACu(I4) = ZERO
          !AKu(I4) = ZERO
          DO J4=1,12
              ACuS(I4,J4) = ZERO
              AKuS(I4,J4) = ZERO
          ENDDO
      ENDDO
      DO I5 = 1,12
          DO J5=1,12
              AKs(I5,J5) = ZERO
              AMs(I5,J5) = ZERO
              ACs(I5,J5) = ZERO
          ENDDO
      ENDDO 
      Zer2(1,1) = ZERO
      Zer2(2,1) = ZERO
      Gauss(1,1)=-SQRT(1/3.0D0)!高斯积分点
      Gauss(2,1)=-SQRT(1/3.0D0)
      Gauss(1,2)=SQRT(1/3.0D0)
      Gauss(2,2)=-SQRT(1/3.0D0)
      Gauss(1,3)=SQRT(1/3.0D0)
      Gauss(2,3)=SQRT(1/3.0D0)
      Gauss(1,4)=-SQRT(1/3.0D0)
      Gauss(2,4)=SQRT(1/3.0D0)
      N = JELEM!单元号？
      !WRITE (*,*) 'NSVARS'
      !WRITE (*,*) NSVARS
      !WRITE (*,*) 'N'
      !WRITE (*,*) N
      !WRITE (*,*) 'COORDS'
      !WRITE (*,*) COORDS
      !WRITE (*,*) 'N'
      !WRITE (*,*) N
      !WRITE (*,*) 'COORDS'
      !WRITE (*,*) COORDS
      !WRITE (*,*) COORDS(1,1)
      !WRITE (*,*) COORDS(2,1)
      !WRITE (*,*) COORDS(1,2)
      !WRITE (*,*) COORDS(2,2)
      !WRITE (*,*) COORDS(1,3)
      !WRITE (*,*) COORDS(2,3)
      !WRITE (*,*) COORDS(1,4)
      !WRITE (*,*) COORDS(2,4)
      EM = 2.D0
      NMAT = Amaterials(N,1)
      IF (NMAT.eq.1) THEN!The material parameters of SBPML elements
          C_ref= 200.D0!Shear wave velocity
          AGs = 80000000.D0!Shear modulus
          Av = 0.30D0!Poisson's ratio
          Ars = 2000.D0!density
      ELSEIF (NMAT.eq.2) THEN
          C_ref= 300.D0
          AGs = 180000000.D0
          Av = 0.30D0
          Ars = 2000.D0
      ENDIF
      
      !C_ref= PROPS(1)
      APHY_LPML=PROPS(2)
      AR=10**-6D0
      b_ele = PROPS(3)*10
      !Av = PROPS(4)!泊松比
      !AGs = PROPS(6)!剪切模量
      AEt = TWO*AGs*(ONE+Av)!杨氏模量
      Alamds = TWO*AGs*Av/(ONE-TWO*Av)
      !Ars = PROPS(5)!密度
      D(1,1) = AEt/((ONE+Av)*(ONE-TWO*Av))*(ONE-Av)
      D(1,2) = AEt/((ONE+Av)*(ONE-TWO*Av))*Av
      D(1,3) = ZERO
      D(2,1) = AEt/((ONE+Av)*(ONE-TWO*Av))*Av
      D(2,2) = AEt/((ONE+Av)*(ONE-TWO*Av))*(ONE-Av)
      D(2,3) = ZERO
      D(3,1) = ZERO
      D(3,2) = ZERO
      D(3,3) = AEt/((ONE+Av)*(ONE-TWO*Av))*(HALF-Av)
      DD = D(1,1)*D(2,2)*D(3,3)-D(1,1)*D(2,3)*D(3,2)
     1 -D(1,2)*D(2,1)*D(3,3)+D(1,2)*D(2,3)*D(3,1)
     2 +D(1,3)*D(2,1)*D(3,2)-D(1,3)*D(2,2)*D(3,1)
      TD(1,1) = (D(2,2)*D(3,3)-D(2,3)*D(3,2))/DD
      TD(1,2) = (D(1,3)*D(3,2)-D(1,2)*D(3,3))/DD
      TD(1,3) = (D(1,2)*D(2,3)-D(1,3)*D(2,2))/DD
      TD(2,1) = (D(2,3)*D(3,1)-D(2,1)*D(3,3))/DD
      TD(2,2) = (D(1,1)*D(3,3)-D(1,3)*D(3,1))/DD
      TD(2,3) = (D(1,3)*D(2,1)-D(1,1)*D(2,3))/DD
      TD(3,1) = (D(2,1)*D(3,2)-D(2,2)*D(3,1))/DD
      TD(3,2) = (D(1,2)*D(3,1)-D(1,1)*D(3,2))/DD
      TD(3,3) = (D(1,1)*D(2,2)-D(1,2)*D(2,1))/DD
      
      
      
      
      Axbsa11 = Aele_first(N,1)!第一变量
      Axbsa22 = Aele_first(N,2)!第二变量
      Coo(1,1) = Axbsa11
      Coo(1,2) = Axbsa11
      Coo(1,3) = Axbsa22
      Coo(1,4) = Axbsa22
      Coo(2,1) = ONE
      Coo(2,2) = -ONE
      Coo(2,3) = -ONE
      Coo(2,4) = ONE
      X1 = Coo(1,1)
      X2 = Coo(1,2)
      X3 = Coo(1,3)
      X4 = Coo(1,4)
      Y1 = Coo(2,1)
      Y2 = Coo(2,2)
      Y3 = Coo(2,3)
      Y4 = Coo(2,4)
      XB1 = AXofN_pm(AInf_NOC_Nb(N,1),1)
      XB2 = AXofN_pm(AInf_NOC_Nb(N,2),1)
      YB1 = AXofN_pm(AInf_NOC_Nb(N,1),2)
      YB2 = AXofN_pm(AInf_NOC_Nb(N,2),2)
      AFL1 = (-COORDS(1,1)+COORDS(1,2))*(-XB1+XB2)+(-COORDS(2,1)
     1 +COORDS(2,2))*(-YB1+YB2)
      !WRITE (*,*) 'AFL1'
      !WRITE (*,*) AFL1
      IF (AFL1.lt.ZERO) THEN
          XB(1,1) = XB2
          XB(2,1) = XB1
          YB(1,1) = YB2
          YB(2,1) = YB1
      ELSE
          XB(1,1) = XB1
          XB(2,1) = XB2
          YB(1,1) = YB1
          YB(2,1) = YB2
      ENDIF
      !WRITE (*,*) 'XB'
      !WRITE (*,*) XB
      !WRITE (*,*) 'YB'
      !WRITE (*,*) YB
      X01 = AXofN_pm(AInf_NOC_Na(N,1),1)
      X02 = AXofN_pm(AInf_NOC_Na(N,2),1)
      Y01 = AXofN_pm(AInf_NOC_Na(N,1),2)
      Y02 = AXofN_pm(AInf_NOC_Na(N,2),2)
      AFL2 = (-COORDS(1,1)+COORDS(1,2))*(-X01+X02)+(-COORDS(2,1)
     1 +COORDS(2,2))*(-Y01+Y02)
      !WRITE (*,*) 'AFL2'
      !WRITE (*,*) AFL2
      IF (AFL2.lt.ZERO) THEN
          X0(1,1) = X02
          X0(2,1) = X01
          Y0(1,1) = Y02
          Y0(2,1) = Y01
      ELSE
          X0(1,1) = X01
          X0(2,1) = X02
          Y0(1,1) = Y01
          Y0(2,1) = Y02
      ENDIF
      !WRITE (*,*) 'X0'
      !WRITE (*,*) X0
      !WRITE (*,*) 'Y0'
      !WRITE (*,*) Y0
      !WRITE (*,*) 'Y0'
      !WRITE (*,*) Y0
      !beta0 = ((EM+1)*C_ref/(2.D0*APHY_LPML)*log(1/ABS(AR)))
      !alfa0 = (EM+1)*b_ele/(2.D0*APHY_LPML)*log(1/ABS(AR))
      beta0 = 310.D0!Attenuation coefficient of waves
      alfa0 = 5.D0!Space tensile coefficient
      !WRITE (*,*) 'beta0'
      !WRITE (*,*) beta0
      !WRITE (*,*) 'alfa0'
      !WRITE (*,*) alfa0
      DO I6=1,4
          XI = Gauss(1,I6)
          ETA = Gauss(2,I6)
          TJ11 = ((1-ETA)*(X2-X1)+(1+ETA)*(X3-X4))/4.D0
          TJ12 = ((1-ETA)*(Y2-Y1)+(1+ETA)*(Y3-Y4))/4.D0
          TJ21 = ((1-XI)*(X4-X1)+(1+XI)*(X3-X2))/4.D0
          TJ22 = ((1-XI)*(Y4-Y1)+(1+XI)*(Y3-Y2))/4.D0
          DJ = TJ11*TJ22-TJ12*TJ21
          !DJ = ABS(DJ)
          !WRITE (*,*) 'TJ11'
          !WRITE (*,*) TJ11
          !WRITE (*,*) 'TJ12'
          !WRITE (*,*) TJ12
          !WRITE (*,*) 'TJ21'
          !WRITE (*,*) TJ21
          !WRITE (*,*) 'TJ22'
          !WRITE (*,*) TJ22
          !WRITE (*,*) 'DJ'
          !WRITE (*,*) DJ
C     %  --- Inversion of JACOBIAN
          AJacob(1,1) = TJ22/DJ
          AJacob(1,2) = -TJ12/DJ
          AJacob(2,1) = -TJ21/DJ
          AJacob(2,2) = TJ11/DJ
C     %  --- Local derivatives of shape function
          g(1, 1) = -(1 - ETA) / 4.D0
          g(2, 1) = -(1 - XI) / 4.D0
          g(1, 2) = (1 - ETA) / 4.D0
          g(2, 2) = -(1 + XI) / 4.D0
          g(1, 3) = (1 + ETA) / 4.D0
          g(2, 3) = (1 + XI) / 4.D0
          g(1, 4) = -(1 + ETA) / 4.D0
          g(2, 4) = (1 - XI) / 4.D0
          Bbb = matmul(AJacob,g)
          !WRITE (*,*) 'AJacob'
          !WRITE (*,*) AJacob
          !WRITE (*,*) 'g'
          !WRITE (*,*) g
C
          ShapeFunM(1,1)=0.25D0*(1-XI)*(1-ETA)
          ShapeFunM(1,2)=0.25D0*(1+XI)*(1-ETA)
          ShapeFunM(1,3)=0.25D0*(1+XI)*(1+ETA)
          ShapeFunM(1,4)=0.25D0*(1-XI)*(1+ETA)
C
          Bb1(1,1) = Bbb(1,1)
          Bb1(1,2) = ZERO
          Bb1(1,3) = Bbb(1,2)
          Bb1(1,4) = ZERO
          Bb1(1,5) = Bbb(1,3)
          Bb1(1,6) = ZERO
          Bb1(1,7) = Bbb(1,4)
          Bb1(1,8) = ZERO
          Bb1(2,1) = Bbb(2,1)
          Bb1(2,2) = ZERO
          Bb1(2,3) = Bbb(2,2)
          Bb1(2,4) = ZERO
          Bb1(2,5) = Bbb(2,3)
          Bb1(2,6) = ZERO
          Bb1(2,7) = Bbb(2,4)
          Bb1(2,8) = ZERO
C
          Bb2(1,1) = ZERO
          Bb2(1,2) = Bbb(1,1)
          Bb2(1,3) = ZERO
          Bb2(1,4) = Bbb(1,2)
          Bb2(1,5) = ZERO
          Bb2(1,6) = Bbb(1,3)
          Bb2(1,7) = ZERO
          Bb2(1,8) = Bbb(1,4)
          Bb2(2,1) = ZERO
          Bb2(2,2) = Bbb(2,1)
          Bb2(2,3) = ZERO
          Bb2(2,4) = Bbb(2,2)
          Bb2(2,5) = ZERO
          Bb2(2,6) = Bbb(2,3)
          Bb2(2,7) = ZERO
          Bb2(2,8) = Bbb(2,4)
C
          ANN1(1,1) = ShapeFunM(1,1)
          ANN1(1,2) = ZERO
          ANN1(1,3) = ShapeFunM(1,2)
          ANN1(1,4) = ZERO
          ANN1(1,5) = ShapeFunM(1,3)
          ANN1(1,6) = ZERO
          ANN1(1,7) = ShapeFunM(1,4)
          ANN1(1,8) = ZERO
C
          ANN2(1,1) = ZERO
          ANN2(1,2) = ShapeFunM(1,1)
          ANN2(1,3) = ZERO
          ANN2(1,4) = ShapeFunM(1,2)
          ANN2(1,5) = ZERO
          ANN2(1,6) = ShapeFunM(1,3)
          ANN2(1,7) = ZERO
          ANN2(1,8) = ShapeFunM(1,4)
          AN12(1,:) = ANN1(1,:)
          AN12(2,:) = ANN2(1,:)
C
          ANN3(1,1) = ShapeFunM(1,1)
          ANN3(1,2) = ZERO
          ANN3(1,3) = ZERO
          ANN3(1,4) = ShapeFunM(1,2)
          ANN3(1,5) = ZERO
          ANN3(1,6) = ZERO
          ANN3(1,7) = ShapeFunM(1,3)
          ANN3(1,8) = ZERO
          ANN3(1,9) = ZERO
          ANN3(1,10) = ShapeFunM(1,4)
          ANN3(1,11) =ZERO
          ANN3(1,12) =ZERO
          ANN3(2,1) = ZERO
          ANN3(2,2) = ShapeFunM(1,1)
          ANN3(2,3) = ZERO
          ANN3(2,4) = ZERO
          ANN3(2,5) = ShapeFunM(1,2)
          ANN3(2,6) = ZERO
          ANN3(2,7) = ZERO
          ANN3(2,8) = ShapeFunM(1,3)
          ANN3(2,9) = ZERO
          ANN3(2,10) = ZERO
          ANN3(2,11) = ShapeFunM(1,4)
          ANN3(2,12) =ZERO
          ANN3(3,1) = ZERO
          ANN3(3,2) = ZERO
          ANN3(3,3) = ShapeFunM(1,1)
          ANN3(3,4) = ZERO
          ANN3(3,5) = ZERO
          ANN3(3,6) = ShapeFunM(1,2)
          ANN3(3,7) = ZERO
          ANN3(3,8) = ZERO
          ANN3(3,9) = ShapeFunM(1,3)
          ANN3(3,10) = ZERO
          ANN3(3,11) = ZERO
          ANN3(3,12) = ShapeFunM(1,4)
C
          Aeta = ZERO
          Axbsa = ZERO
          DO I7 = 1,4
              Aeta = Aeta + ShapeFunM(1,I7)*Coo(2,I7)
              Axbsa = Axbsa + ShapeFunM(1,I7)*Coo(1,I7)
          ENDDO
C
          AN1 = HALF*(ONE+Aeta)
          AN2 = HALF*(ONE-Aeta)
          ANeta(1,1) = AN1
          ANeta(1,2) = AN2
          Adiff_eta(1,1) = HALF
          Adiff_eta(1,2) = -HALF
          Aj11 = matmul(ANeta,(XB-X0))*matmul(Adiff_eta,Y0)
          Aj12 = matmul(ANeta,(XB-X0))*matmul(Adiff_eta,YB)
          Aj21 = matmul(ANeta,(YB-Y0))*matmul(Adiff_eta,X0)
          Aj22 = matmul(ANeta,(YB-Y0))*matmul(Adiff_eta,XB)
          Aa_eta = Aj11(1,1)-Aj21(1,1)
          Ab_eta = -Aa_eta+Aj12(1,1)-Aj22(1,1)
          TEMP = matmul(Adiff_eta,Y0)
          A_eta(1,1) = TEMP(1,1)
          TEMP = matmul(-ANeta,(YB-Y0))
          A_eta(1,2) = TEMP(1,1)
          TEMP = matmul(-Adiff_eta,X0)
          A_eta(2,1) = TEMP(1,1)
          TEMP = matmul(ANeta,XB-X0)
          A_eta(2,2) = TEMP(1,1)
          TEMP = matmul(Adiff_eta,YB-Y0)
          B_eta(1,1) = TEMP(1,1)
          B_eta(1,2) = ZERO
          TEMP = matmul(-Adiff_eta,XB-X0)
          B_eta(2,1) = TEMP(1,1)
          B_eta(2,2) = ZERO
          XX0 = 0.D0
          alfa_xbsa = ONE+alfa0*((Axbsa-XX0)*Ans/ALpml)**EM
          beta_xbsa = beta0*((Axbsa-XX0)*Ans/ALpml)**EM
          alfa_xb11 = Axbsa+alfa0*(Ans/ALpml)**EM*(Axbsa-XX0)
     1    **(EM+1)/(EM+1)
          beta_xb11 = beta0*(Ans/ALpml)**EM*(Axbsa-XX0)**(EM+1)
     2    /(EM+1)
          Ea(1,1) = 1
          Ea(1,2) = ZERO
          Ea(2,1) = ZERO
          Ea(2,2) = alfa_xbsa
          Eb(1,1) = ZERO
          Eb(2,1) = ZERO
          Eb(1,2) = ZERO
          Eb(2,2) = beta_xbsa
          Fa = alfa_xb11*EYE
          Fb = beta_xb11*EYE
C
          a11 = Aa_eta*alfa_xbsa+Ab_eta*alfa_xb11*alfa_xbsa
          a22 = Aa_eta*beta_xbsa+alfa_xb11*beta_xbsa*Ab_eta+
     1    Ab_eta*beta_xb11*alfa_xbsa
          a33 = Ab_eta*beta_xb11*beta_xbsa
C
          A1(1,:) = A_eta(1,:)
          A2(1,:) = A_eta(2,:)
          B1(1,:) = B_eta(1,:)
          B2(1,:) = B_eta(2,:)
C
          TEMP_CuS1 = matmul(transpose(Bb1),transpose(matmul(A1,Ea)+
     1     matmul(matmul(B1,Fa),Ea)))
          TEMP_CuS2 = matmul(transpose(Bb2),transpose(matmul(A1,Ea)+
     1     matmul(matmul(B1,Fa),Ea)))
          !WRITE (*,*) 'Bb1'
          !WRITE (*,*) Bb1
          !WRITE (*,*) 'A1'
          !WRITE (*,*) A1
          !WRITE (*,*) 'A2'
          !WRITE (*,*) A2
          !WRITE (*,*) 'A_eta'
          !WRITE (*,*) A_eta
          !WRITE (*,*) 'Y0'
          !WRITE (*,*) Y0
          !WRITE (*,*) 'YB'
          !WRITE (*,*) YB
          !WRITE (*,*) 'Ea'
          !WRITE (*,*) Ea
          !WRITE (*,*) 'B1'
          !WRITE (*,*) B1
          !WRITE (*,*) 'Fa'
          !WRITE (*,*) Fa
          !WRITE (*,*) 'Ea'
          !WRITE (*,*) Ea
          !WRITE (*,*) 'Bb2'
          !WRITE (*,*) Bb2
          TEMP_CuS3(:,1) = TEMP_CuS1(:,1)
          TEMP_CuS3(:,2) = TEMP_CuS2(:,1)
          TEMP_CuS4 = matmul(transpose(Bb1),transpose(matmul(A2,Ea)+
     1     matmul(matmul(B2,Fa),Ea)))
          TEMP_CuS5 = matmul(transpose(Bb2),transpose(matmul(A2,Ea)+
     1     matmul(matmul(B2,Fa),Ea)))
          TEMP_CuS6(:,1) = TEMP_CuS4(:,1)
          TEMP_CuS6(:,2) = TEMP_CuS5(:,1)
          DO I8=1,3
              DO J8=1,2
                  AL1(I8,J8) = ZERO
                  AL2(I8,J8) = ZERO
              ENDDO
          ENDDO
          AL1(1,1) = ONE
          AL1(3,2) = ONE
          AL2(2,2) = ONE
          AL2(3,1) = ONE
          ACuS = ACuS+matmul(matmul(TEMP_CuS3,transpose(AL1)),
     1     ANN3)*DJ+matmul(matmul(TEMP_CuS6,transpose(AL2)),
     2     ANN3)*DJ
          !WRITE (*,*) 'TEMP_CuS3'
          !WRITE (*,*) TEMP_CuS3
          !WRITE (*,*) 'AL1'
          !WRITE (*,*) AL1
          !WRITE (*,*) 'ANN3'
          !WRITE (*,*) ANN3
          !WRITE (*,*) 'DJ'
          !WRITE (*,*) DJ
          !WRITE (*,*) 'TEMP_CuS6'
          !WRITE (*,*) TEMP_CuS6
          !WRITE (*,*) 'AL2'
          !WRITE (*,*) AL2
          !WRITE (*,*) 'ACuS'
          !WRITE (*,*) ACuS
C
          TEMP_KuS1 = matmul(transpose(Bb1),transpose(matmul(A1,Eb)+
     1     matmul(matmul(B1,Fa),Eb)+matmul(matmul(B1,Fb),Ea)))
          TEMP_KuS2 = matmul(transpose(Bb2),transpose(matmul(A1,Eb)+
     1     matmul(matmul(B1,Fa),Eb)+matmul(matmul(B1,Fb),Ea)))
          TEMP_KuS3(:,1) = TEMP_KuS1(:,1)
          TEMP_KuS3(:,2) = TEMP_KuS2(:,1)
          TEMP_KuS4 = matmul(transpose(Bb1),transpose(matmul(A2,Eb)+
     1     matmul(matmul(B2,Fa),Eb)+matmul(matmul(B2,Fb),Ea)))
          TEMP_KuS5 = matmul(transpose(Bb2),transpose(matmul(A2,Eb)+
     1     matmul(matmul(B2,Fa),Eb)+matmul(matmul(B2,Fb),Ea)))
          TEMP_KuS6(:,1) = TEMP_KuS4(:,1)
          TEMP_KuS6(:,2) = TEMP_KuS5(:,1)
          AKuS = AKuS+matmul(matmul(TEMP_KuS3,transpose(AL1)),
     1     ANN3)*DJ+matmul(matmul(TEMP_KuS6,transpose(AL2)),
     2     ANN3)*DJ
          !WRITE (*,*) 'AKuS'
          !WRITE (*,*) AKuS
C
          AMu = AMu+Ars*a11*DJ*matmul(transpose(AN12),AN12)
          ACu = ACu+Ars*a22*DJ*matmul(transpose(AN12),AN12)
          AKu = AKu+Ars*a33*DJ*matmul(transpose(AN12),AN12)
C
          AMs = AMs+a11*DJ*matmul(matmul(transpose(ANN3),TD),ANN3)
          ACs = ACs+a22*DJ*matmul(matmul(transpose(ANN3),TD),ANN3)
          AKs = AKs+a33*DJ*matmul(matmul(transpose(ANN3),TD),ANN3)
          !WRITE (*,*) 'AMu'
          !WRITE (*,*) AMu 
          !WRITE (*,*) 'ACu'
          !WRITE (*,*) ACu 
          !WRITE (*,*) 'AKu'
          !WRITE (*,*) AKu 
          !WRITE (*,*) 'AMs'
          !WRITE (*,*) AMs 
          !WRITE (*,*) 'ACs'
          !WRITE (*,*) ACs 
          !WRITE (*,*) 'AKs'
          !WRITE (*,*) AKs
          !WRITE (*,*) 'AN12'
          !WRITE (*,*) AN12
          !WRITE (*,*) 'ANN3'
          !WRITE (*,*) ANN3
          !WRITE (*,*) 'TD'
          !WRITE (*,*) TD
          !WRITE (*,*) 'a11'
          !WRITE (*,*) a11
          !WRITE (*,*) 'DJ'
          !WRITE (*,*) DJ
          AN1212 = matmul(transpose(AN12),AN12)
          !WRITE (*,*) 'AN1212'
          !WRITE (*,*) AN1212
      ENDDO
      ACSu = transpose(ACuS)
      AKSu = transpose(AKuS)
C     检查主对线元素是否全为负
      FLAG1 = ZERO
      DO I17=1,8
          IF (AMu(I17,I17).ge.ZERO) THEN
              FLAG1=FLAG1+ONE
          ENDIF
      ENDDO
      IF (FLAG1.eq.ZERO) THEN
          AMu=-AMu
          ACu=-ACu
          AKu=-AKu
          ACuS=-ACuS
          AKuS=-AKuS
      ENDIF
      FLAG2=ZERO
      DO I18=1,8
          IF (AMs(I18,I18).ge.ZERO) THEN
              FLAG2=FLAG2+ONE
          ENDIF
      ENDDO
      IF (FLAG2.eq.ZERO) THEN
          AMs=-AMs
          ACs=-ACs
          AKs=-AKs
          ACSu=-ACSu
          AKSu=-AKSu
      ENDIF
      !WRITE (*,*) 'AMu'
      !WRITE (*,*) AMu 
      !WRITE (*,*) 'ACu'
      !WRITE (*,*) ACu 
      !WRITE (*,*) 'AKu'
      !WRITE (*,*) AKu 
      !WRITE (*,*) 'AMs'
      !WRITE (*,*) AMs 
      !WRITE (*,*) 'ACs'
      !WRITE (*,*) ACs 
      !WRITE (*,*) 'AKs'
      !WRITE (*,*) AKs
      !WRITE (*,*) 'ACuS'
      !WRITE (*,*) ACuS
      !WRITE (*,*) 'ACSu'
      !WRITE (*,*) ACSu
C
      DO I9=1,8
          DO J9=1,8
              AM1(I9,J9) = AMu(I9,J9)
              AK1(I9,J9) = AKu(I9,J9)
              AC1(I9,J9) = ACu(I9,J9)
          ENDDO
      ENDDO
      
      DO I10=1,12
          DO J10=1,12
              AM1(I10+8,J10+8) = AMs(I10,J10)
              AK1(I10+8,J10+8) = AKs(I10,J10)
              AC1(I10+8,J10+8) = ACs(I10,J10)
          ENDDO
      ENDDO
      !WRITE (*,*) 'AM1'
      !WRITE (*,*) AM1
C
      DO I12=1,8
          DO J12=1,12
              AK1(I12,J12+8) = AKuS(I12,J12)
              AK1(J12+8,I12) = -AKSu(J12,I12)
              AC1(I12,J12+8) = ACuS(I12,J12)
              AC1(J12+8,I12) = -ACSu(J12,I12)
          ENDDO
      ENDDO
C
      DO I13=1,20
          DO J13=1,20
              T1(I13,J13)=ZERO
              T2(I13,J13)=ZERO
              T3(I13,J13)=ZERO
              T4(I13,J13)=ZERO
          ENDDO
      ENDDO
      T1(1,1)=ONE
      T1(2,2)=ONE
      T1(3,9)=ONE
      T1(4,10)=ONE
      T1(5,11)=ONE
      T1(6,3)=ONE
      T1(7,4)=ONE
      T1(8,12)=ONE
      T1(9,13)=ONE
      T1(10,14)=ONE
      T1(11,5)=ONE
      T1(12,6)=ONE
      T1(13,15)=ONE
      T1(14,16)=ONE
      T1(15,17)=ONE
      T1(16,7)=ONE
      T1(17,8)=ONE
      T1(18,18)=ONE
      T1(19,19)=ONE
      T1(20,20)=ONE
C
      T2=transpose(T1)
      AM2 = matmul(matmul(T1,AM1),T2)
      AK2 = matmul(matmul(T1,AK1),T2)
      AC2 = matmul(matmul(T1,AC1),T2)
      !WRITE (*,*) 'AM2'
      !WRITE (*,*) AM2
C
      T3(1,1)=ONE
      T3(2,2)=ONE
      T3(3,3)=ONE
      T3(4,4)=ONE
      T3(5,5)=ONE
      T3(6,16)=ONE
      T3(7,17)=ONE
      T3(8,18)=ONE
      T3(9,19)=ONE
      T3(10,20)=ONE
      T3(11,11)=ONE
      T3(12,12)=ONE
      T3(13,13)=ONE
      T3(14,14)=ONE
      T3(15,15)=ONE
      T3(16,6)=ONE
      T3(17,7)=ONE
      T3(18,8)=ONE
      T3(19,9)=ONE
      T3(20,10)=ONE
C
      T4=transpose(T3)
      !AM3 = matmul(matmul(T3,AM2),T4)
      !AK3 = matmul(matmul(T3,AK2),T4)
      !AC3 = matmul(matmul(T3,AC2),T4)
      
      !WRITE (*,*) 'TIME'
      !WRITE (*,*) TIME(1)
      !AT   = PROPS(1)!厚度
      !E    = PROPS(2)!弹性模量
      !Av  = PROPS(3)!泊松比
      !RHO  = PROPS(4)!密度
      !WRITE (*,*) 'COORDS'
      !WRITE (*,*) COORDS
      !Av = ZERO!泊松比
      !AGs = HALF!剪切模量
      
C
      DO K1 = 1, NDOFEL  
          !AN(1,K1) = ZERO!形函数
          !AN(2,K1) = ZERO
          !TN(K1,1) = ZERO!形函数转置
          !TN(K1,2) = ZERO
          !SRESID(K1) = ZERO!外力
          !DO KRHS = 1, NRHS
          !RHS(K1,1) = ZERO!构造RHS零矩阵
          !ENDDO
          !SVARS(K1)=ZERO
          !SVARS(K1+8)=ZERO
          DO K2 = 1, NDOFEL
              AMATRX(K2,K1) = ZERO!构造M、C、K零矩阵
              AK(K2,K1) = ZERO
              AM(K2,K1) = ZERO
              !AM1(K2,K1) = ZERO
              !AMM(K2,K1) = ZERO
              !Km(K2,K1) = ZERO
              CC(K2,K1) = ZERO
              !Kmmmm(K2,K1) = ZERO
              !Kmmm(K2,K1) = ZERO
          ENDDO
      ENDDO
      AM = AM2
      CC = AC2
      AK = AK2
      !WRITE (*,*) 'AM'
      !WRITE (*,*) AM
      !WRITE (*,*) 'AK'
      !WRITE (*,*) AK
      !WRITE (*,*) 'CC'
      !WRITE (*,*) CC
      RETURN
      END

!     Store stiff/mass matrices
      SUBROUTINE STOREMATRICES(SVARS,NSVARS,AK,AM,CC,NDOFEL)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION SVARS(NSVARS), AK(NDOFEL,NDOFEL), AM(NDOFEL,NDOFEL),
     1 CC(NDOFEL,NDOFEL)
      DO I = 1,NDOFEL
          SVARS((I+1)*NDOFEL+1:(I+2)*NDOFEL) = AK(I,1:NDOFEL)      
          SVARS((NDOFEL+I+1)*NDOFEL+1:(NDOFEL+I+2)*NDOFEL) = 
     1     AM(I,1:NDOFEL)
          SVARS((NDOFEL*2+I+1)*NDOFEL+1:(NDOFEL*2+I+2)*NDOFEL) = 
     1     CC(I,1:NDOFEL)
      ENDDO
      RETURN
      END
C
      !     Read stiff/mass matrices      
      SUBROUTINE READMATRICES(SVARS,NSVARS,AK,AM,CC,NDOFEL)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION SVARS(NSVARS), AK(NDOFEL,NDOFEL), AM(NDOFEL,NDOFEL),
     1 CC(NDOFEL,NDOFEL)
      
      DO I = 1,NDOFEL
      AK(I,1:NDOFEL) = SVARS((I+1)*NDOFEL+1:(I+2)*NDOFEL)     
      AM(I,1:NDOFEL) = SVARS((NDOFEL+I+1)*NDOFEL+1:(NDOFEL+I+2)*NDOFEL)
      CC(I,1:NDOFEL) = SVARS((NDOFEL*2+I+1)
     1 *NDOFEL+1:(NDOFEL*2+I+2)*NDOFEL)
      ENDDO
      
      
      RETURN
      
      END
C
      SUBROUTINE OUTPUTVARIABLE(RHS,AMATRX,SVARS,PROPS,ENERGY,U,V,A,
     1 LFLAGS,DTIME,NDOFEL,NRHS,NSVARS,MLVARX,JELEM,PARAMS,AK,AM,CC)
      
      INCLUDE 'ABA_PARAM.INC'
      
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL)
      DIMENSION SVARS(NSVARS),LFLAGS(*)  
      DIMENSION U(NDOFEL),V(NDOFEL),A(NDOFEL)
      DIMENSION PARAMS(3),ENERGY(8),PROPS(*),SRESID(NDOFEL),GT(NDOFEL),
     1 CC(NDOFEL,NDOFEL), AK(NDOFEL,NDOFEL), AM(NDOFEL,NDOFEL)
      PARAMETER ( ZERO = 0.D0, HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0)
      DO K1=1,NDOFEL
       DO K2 = 1,NRHS
            RHS(K1,K2) = 0.D0
       ENDDO
      ENDDO
      
      IF (LFLAGS(3).EQ.1) THEN     
C       Normal Incrementation
        
        IF (LFLAGS(1).EQ.1 .OR. LFLAGS(1).EQ.2) THEN
C         *STATIC                                        
          AMATRX = AK 
          SRESID = MATMUL(AK,U)   
	    RHS(1:NDOFEL,1)=RHS(1:NDOFEL,1)-SRESID                          
 

        ELSE IF (LFLAGS(1).EQ.11 .OR. LFLAGS(1).EQ.12) THEN
C         *DYNAMIC
          ALPHA = PARAMS(1)
          BETA  = PARAMS(2)
          GAMMA = PARAMS(3) 
          !WRITE (*,*) 'ALPHA'
          !WRITE (*,*) ALPHA
          !WRITE (*,*) 'BETA'
          !WRITE (*,*) BETA
          !WRITE (*,*) 'GAMMA'
          !WRITE (*,*) GAMMA
          !WRITE (*,*) '112'
          !WRITE (*,*) 'DTIME'
          !WRITE (*,*) DTIME
          DADU = 1.D0/(BETA*DTIME**2)
          DVDU = GAMMA/(BETA*DTIME)
          !WRITE (*,*) 'DADU'
          !WRITE (*,*) DADU
          !WRITE (*,*) 'DVDU'
          !WRITE (*,*) DVDU
          !WRITE (*,*) 'AM'
          !WRITE (*,*) AM
          !WRITE (*,*) 'CC'
          !WRITE (*,*) CC
          !WRITE (*,*) 'AK'
          !WRITE (*,*) AK
          !WRITE (*,*) 'GT'
          !WRITE (*,*) GT
          AMATRX = DADU*AM+(1.D0+ALPHA)*DVDU*CC+(1.D0+ALPHA)*AK
          !WRITE (*,*) 'AMATRX'
          !WRITE (*,*) AMATRX
          GT =  SVARS(1:NDOFEL)    
          SRESID = MATMUL(AM,A) + (1.D0+ALPHA)*MATMUL(AK,U)
     *              +(1.D0+ALPHA)*MATMUL(CC,V) - ALPHA*GT
          !WRITE (*,*) 'A'
          !WRITE (*,*) A
          !WRITE (*,*) 'V'
          !WRITE (*,*) V
          !WRITE (*,*) 'U'
          !WRITE (*,*) U
          !WRITE (*,*) 'SRESID'
          !WRITE (*,*) SRESID
	    RHS(1:NDOFEL,1)=RHS(1:NDOFEL,1)-SRESID
          !WRITE (*,*) 'RHS'
          !WRITE (*,*) RHS
          SRESID = MATMUL(AK,U) + MATMUL(CC,V)
          SVARS(NDOFEL+1:2*NDOFEL) = SVARS(1:NDOFEL)
          SVARS(1:NDOFEL) = SRESID

          
          ENERGY(1) = 0.0D0
          
          DO I=1,NDOFEL
              DO J=1,NDOFEL
                 ENERGY(1) = ENERGY(1) + 0.5D0*V(I)*AM(I,J)*V(J)
              ENDDO    
          ENDDO
                   
          DO I=1,NDOFEL
              DO J=1,NDOFEL
                 ENERGY(2) = ENERGY(2) + 0.5D0*U(I)*AK(I,J)*U(J)
              ENDDO    
          ENDDO
        !WRITE (*,*) 'AMATRX'
        !WRITE (*,*) AMATRX 
        END IF

      ELSE IF (LFLAGS(3).EQ.2) THEN 
C      Stiffness matrix 
      AMATRX=AK
      !WRITE (*,*) '32'        
      ELSE IF (LFLAGS(3).EQ.3) THEN 
C      Damping matrix     
       AMATRX=CC
     	!WRITE (*,*) '33'			   
      ELSE IF (LFLAGS(3).EQ.4) THEN 
C     Mass matrix 
       AMATRX=AM
       !WRITE (*,*) '34'
      ELSE IF (LFLAGS(3).EQ.5) THEN 
C       Half-step residual calculation
       ALPHA = PARAMS(1) 
      !WRITE (*,*) '35'
       GT =  0.5D0*(SVARS(1:NDOFEL)+SVARS(NDOFEL+1:2*NDOFEL))

       
       SRESID = MATMUL(AM,A)+(1.D0+ALPHA)*MATMUL(AK,U)
     *          +(1.D0+ALPHA)*MATMUL(CC,V)- ALPHA*GT
          
	 RHS(1:NDOFEL,1) = RHS(1:NDOFEL,1)-SRESID
 
      ELSE IF (LFLAGS(3).EQ.6) THEN
C       Initial acceleration calculation
       AMATRX=AM 
       !WRITE (*,*) '36'
      SRESID = MATMUL(AK,U)+MATMUL(CC,V)    
	RHS(1:NDOFEL,1) = RHS(1:NDOFEL,1)-SRESID
      SVARS(1:NDOFEL) = SRESID

       
       DO I=1,NDOFEL
          DO J=1,NDOFEL
             ENERGY(1) = ENERGY(1) + 0.5D0*V(I)*AM(I,J)*V(J)
          ENDDO    
       ENDDO
              
       DO I=1,NDOFEL
          DO J=1,NDOFEL
             ENERGY(2) = ENERGY(2) + 0.5D0*U(I)*AK(I,J)*U(J)
          ENDDO    
       ENDDO

      ELSE
	  !WRITE(17,*) "**** LFLAGS(1)=,",LFLAGS(1)
	  !WRITE(17,*) "**** LFLAGS(2)=,",LFLAGS(2)
	  !WRITE(17,*) "**** LFLAGS(3)=,",LFLAGS(3)
	  !WRITE(17,*) "**** LFLAGS(4)=,",LFLAGS(4)
	  !WRITE(17,*) "**** LFLAGS(5)=,",LFLAGS(5)
        !WRITE(17,*) "Error: Temporary unsupported function"

      END IF
      
      RETURN
      END
C
      SUBROUTINE FORMMATRXCPE4(AK,AM,CC,NDOFEL,NNODE,NSVARS,MCRD,
     1 MLVARX,MDLOAD,NPREDF,JELEM,COORDS,PROPS)
C
      INCLUDE 'ABA_PARAM.INC'
C
      PARAMETER ( ZERO = 0.D0, HALF = 0.5D0, ONE = 1.D0, TWO = 2.D0)
      DIMENSION RHS(MLVARX,2),AMATRX(NDOFEL,NDOFEL),
     1     SVARS(NSVARS),ENERGY(8),PROPS(4),COORDS(MCRD,NNODE),
     2     U(NDOFEL),V(NDOFEL),A(NDOFEL),TIME(2),
     3     PARAMS(3),
     4     PREDEF(2,NPREDF,NNODE),
     5     SRESID(8),AK(NDOFEL,NDOFEL),AM(NDOFEL,NDOFEL),
     6     CC(NDOFEL,NDOFEL),Gauss(4,2),D(3,3),dNdrs(NNODE,2),
     7     aJacob(2,2),TaJacob(2,2),dNdxy(NNODE,2),Bm(3,8),
     8     AKmm(8,3),AKmmm(8,8),AKmmmm(8,8),GT(NDOFEL),AN(2,8),
     9     TN(8,2),AMM(NDOFEL,NDOFEL),BmT(8,3),AM1(8,8)
C     这个CPE4采用了4点高斯积分和集中质量阵
      Aheight = (COORDS(2,1)+COORDS(2,2)+COORDS(2,3)+COORDS(2,4))/4.D0
      IF (Aheight.gt.70.D0) THEN
          E = 208000000.D0
      ELSE IF (Aheight.lt.70.D0 .and. Aheight.gt.40.D0) THEN
          E = 468000000.D0
      ELSE IF (Aheight.lt.40.D0) THEN
          E = 832000000.D0
      ENDIF
      AT   = PROPS(1)!厚度
      !E    = PROPS(2)!弹性模量
      ANU  = PROPS(3)!泊松比
      RHO  = PROPS(4)!密度
      DO K1 = 1, NDOFEL  
          AN(1,K1) = ZERO!形函数
          AN(2,K1) = ZERO
          TN(K1,1) = ZERO!形函数转置
          TN(K1,2) = ZERO
          SRESID(K1) = ZERO!外力
          DO KRHS = 1, NRHS
              RHS(K1,KRHS) = ZERO!构造RHS零矩阵
          ENDDO
          SVARS(K1)=ZERO
          DO K2 = 1, NDOFEL
              AMATRX(K2,K1) = ZERO!构造M、C、K零矩阵
              AK(K2,K1) = ZERO
              AM(K2,K1) = ZERO
              AM1(K2,K1) = ZERO
              AMM(K2,K1) = ZERO
              !Km(K2,K1) = ZERO
              CC(K2,K1) = ZERO
              AKmmmm(K2,K1) = ZERO
              AKmmm(K2,K1) = ZERO
          ENDDO
      END DO
      Gauss(1,1)=-SQRT(1/3.0D0)!按照第3->4->1->2象限的顺序
      Gauss(1,2)=-SQRT(1/3.0D0)
      Gauss(2,1)=SQRT(1/3.0D0)
      Gauss(2,2)=-SQRT(1/3.0D0)
      Gauss(3,1)=SQRT(1/3.0D0)
      Gauss(3,2)=SQRT(1/3.0D0)
      Gauss(4,1)=-SQRT(1/3.0D0)
      Gauss(4,2)=SQRT(1/3.0D0)
      !WRITE (*,*) 'Gauss'
      !WRITE (*,*) Gauss(1,1:2)
      !WRITE (*,*) Gauss(2,1:2)
      !WRITE (*,*) Gauss(3,1:2)
      !WRITE (*,*) Gauss(4,1:2)
      D(1,1)=(ONE-ANU)*E*(ONE/((ONE+ANU)*(ONE-2.D0*ANU)))
      D(2,2)=(ONE-ANU)*E*(ONE/((ONE+ANU)*(ONE-2.D0*ANU)))
      D(1,2)=ANU*E*(ONE/((ONE+ANU)*(ONE-2.D0*ANU)))
      D(2,1)=ANU*E*(ONE/((ONE+ANU)*(ONE-2.D0*ANU)))
      D(1,3)=ZERO
      D(2,3)=ZERO
      D(3,1)=ZERO
      D(3,2)=ZERO
      D(3,3)=(ONE-2.D0*ANU)/TWO*E*(ONE/((ONE+ANU)*(ONE-2.D0*ANU)))
      !WRITE (*,*) 'D'
      !WRITE (*,*) D
      !WRITE (*,*) 'COORDS'
      !WRITE (*,*) COORDS
C
      DO I = 1,4
          dNdrs(1,1) = -(ONE-Gauss(I,2))*0.25D0
          dNdrs(2,1) = (ONE-Gauss(I,2))*0.25D0
          dNdrs(3,1) = (ONE+Gauss(I,2))*0.25D0
          dNdrs(4,1) = -(ONE+Gauss(I,2))*0.25D0
          dNdrs(1,2) = -(ONE-Gauss(I,1))*0.25D0
          dNdrs(2,2) = -(ONE+Gauss(I,1))*0.25D0
          dNdrs(3,2) = (ONE+Gauss(I,1))*0.25D0
          dNdrs(4,2) = (ONE-Gauss(I,1))*0.25D0
          !WRITE (*,*) 'dNdrs'
          !WRITE (*,*) dNdrs
          aJacob=matmul(COORDS(1:2,1:4),dNdrs)      !雅各比矩阵
          !WRITE (*,*) 'aJacob'
          !WRITE (*,*) aJacob
C雅克比矩阵求逆
          TaJacob(1,1)=aJacob(2,2)/(aJacob(1,1)*aJacob(2,2)
     1                -aJacob(1,2)*aJacob(2,1))
          TaJacob(1,2)=-aJacob(1,2)/(aJacob(1,1)*aJacob(2,2)
     1                -aJacob(1,2)*aJacob(2,1))
          TaJacob(2,1)=-aJacob(2,1)/(aJacob(1,1)*aJacob(2,2)
     1                -aJacob(1,2)*aJacob(2,1))
          TaJacob(2,2)=aJacob(1,1)/(aJacob(1,1)*aJacob(2,2)
     1                -aJacob(1,2)*aJacob(2,1))
          !WRITE (*,*) 'TaJacob'
          !WRITE (*,*) TaJacob
          dNdxy=matmul(dNdrs,TaJacob)
          !WRITE (*,*) 'dNdxy'
          !WRITE (*,*) dNdxy
          DO I3=1,3
              DO J3=1,8
                  Bm(I3,J3)=ZERO
              ENDDO
          ENDDO
C
          DO J = 1,4
              Bm(1, (J-1)*2+1) = dNdxy(J, 1)
              Bm(2, (J-1)*2+2) = dNdxy(J, 2)
              Bm(3, (J-1)*2+1) = dNdxy(J, 2)
              Bm(3, (J-1)*2+2) = dNdxy(J, 1)
          END DO
          !WRITE (*,*) 'Bm'
          !WRITE (*,*) Bm
          wt=1.D0                                 !积分点权重
          !BmT = transpose(Bm)
          DO I1 = 1,8
              DO J1 = 1,3
                  BmT(I1,J1) = Bm(J1,I1)
              ENDDO
          ENDDO
          !WRITE (*,*) 'BmT1'
          !WRITE (*,*) BmT(1,1:3)
          !WRITE (*,*) 'BmT2'
          !WRITE (*,*) BmT(2,1:3)
          !WRITE (*,*) 'BmT3'
          !WRITE (*,*) BmT(3,1:3)
          AKmm=matmul(BmT,D)
          AKmmm=matmul(AKmm,Bm)
          aJacob_det=aJacob(1,1)*aJacob(2,2)
     1        -aJacob(1,2)*aJacob(2,1)
          !WRITE (*,*) 'matmul(BmT,D)'
          !WRITE (*,*) matmul(BmT,D)
          !WRITE (*,*) 'aJacob_det'
          !WRITE (*,*) aJacob_det
          !WRITE (*,*) 'AT'
          !WRITE (*,*) AT
          !WRITE (*,*) 'wt'
          !WRITE (*,*) wt
          !WRITE (*,*) 'AKmmm'
          !WRITE (*,*) AKmmm
          !WRITE (*,*) 'AKmm'
          !WRITE (*,*) AKmm
          AKmmmm=AKmmm*wt*AT*aJacob_det
          AK = AK+AKmmmm
C
          AN1 = (ONE-Gauss(I,1))*(ONE-Gauss(I,2))*0.25D0
          AN2 = (ONE+Gauss(I,1))*(ONE-Gauss(I,2))*0.25D0
          AN3 = (ONE+Gauss(I,1))*(ONE+Gauss(I,2))*0.25D0
          AN4 = (ONE-Gauss(I,1))*(ONE+Gauss(I,2))*0.25D0
          AN(1,1) = AN1
          AN(2,2) = AN1
          AN(1,3) = AN2
          AN(2,4) = AN2
          AN(1,5) = AN3
          AN(2,6) = AN3
          AN(1,7) = AN4
          AN(2,8) = AN4
          !WRITE (*,*) 'AN'
          !WRITE (*,*) AN
          DO I2 = 1,8
              DO J2 = 1,2
                  TN(I2,J2) = AN(J2,I2)
              ENDDO
          ENDDO
          !WRITE (*,*) 'TN'
          !WRITE (*,*) TN
          AMM = RHO*matmul(TN,AN)*AT*aJacob_det*wt
          AM1 = AM1 + AMM
      END DO
C     集中质量阵
      DO I3 = 1,8
          DO J3 = 1,8
              AM(I3,I3) = AM(I3,I3)+AM1(I3,J3)
          ENDDO
      ENDDO
      RETURN
      END
C
      SUBROUTINE STOREMATRICESCPE4(SVARS,NSVARS,AK,AM,CC,NDOFEL)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION SVARS(NSVARS), AK(NDOFEL,NDOFEL), AM(NDOFEL,NDOFEL),
     1 CC(NDOFEL,NDOFEL)
      DO I = 1,NDOFEL
          SVARS((I+1)*NDOFEL+1:(I+2)*NDOFEL) = AK(I,1:NDOFEL)      
          SVARS((NDOFEL+I+1)*NDOFEL+1:(NDOFEL+I+2)*NDOFEL) = 
     1     AM(I,1:NDOFEL)
          SVARS((NDOFEL*2+I+1)*NDOFEL+1:(NDOFEL*2+I+2)*NDOFEL) = 
     1     CC(I,1:NDOFEL)
      ENDDO
      RETURN
      END
C
      SUBROUTINE READMATRICESCPE4(SVARS,NSVARS,AK,AM,CC,NDOFEL)
      
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION SVARS(NSVARS), AK(NDOFEL,NDOFEL), AM(NDOFEL,NDOFEL),
     1 CC(NDOFEL,NDOFEL)
      
      DO I = 1,NDOFEL
      AK(I,1:NDOFEL) = SVARS((I+1)*NDOFEL+1:(I+2)*NDOFEL)     
      AM(I,1:NDOFEL) = SVARS((NDOFEL+I+1)*NDOFEL+1:(NDOFEL+I+2)*NDOFEL)
      CC(I,1:NDOFEL) = SVARS((NDOFEL*2+I+1)
     1 *NDOFEL+1:(NDOFEL*2+I+2)*NDOFEL)
      ENDDO
      RETURN
      END