      SUBROUTINE UMAT(STRESS, STATEV, DDSDDE, SSE, SPD, SCD, RPL,
     1 DDSDDT, DRPLDE, DRPLDT, STRAN, DSTRAN, TIME, DTIME, TEMP, DTEMP,
     2 PREDEF, DPRED, CMNAME, NDI, NSHR, NTENS, NSTATV, PROPS, NPROPS,
     3 COORDS, DROT, PNEWDT, CELENT, DFGRD0, DFGRD1, NOEL, NPT, LAYER,
     4 KSPT, KSTEP, KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*8 CMNAME
C
      DIMENSION STRESS(NTENS), STATEV(NSTATV), DDSDDE(NTENS, NTENS),
     1 DDSDDT(NTENS), DRPLDE(NTENS), STRAN(NTENS), DSTRAN(NTENS),
     2 PREDEF(1), DPRED(1), PROPS(NPROPS), COORDS(3), DROT(3, 3),
     3 DFGRD0(3, 3),DFGRD1(3,3),IDENTITY(NTENS,NTENS),XIDENTITY(NTENS),
     4 StressTER(NTENS),STRILALDEV(NTENS),STRESSDEV(NTENS),
     5 TERM1(NTENS,NTENS),TERM2(NTENS,NTENS),TERM3(NTENS,NTENS),
     6 TERM4(NTENS,NTENS),DELTASTRAINP(NTENS),DELTASTRAINE(NTENS)
C
      PARAMETER(ZERO=0.D0, ONE=1.D0, TWO=2.D0, THREE=3.D0, SIX=6.D0,
     1 ENUMAX=.4999D0, TOLER=1.0D-6, FOUR=4.D0, THO=(2/3)**0.5,
     2 THREETWO=(3/2)**0.5,TOLERANCE=0.001)
C
C ----------------------------------------------------------------
C UMAT FOR ISOTROPIC ELASTICITY AND MISES PLASTICITY
C WITH NON-LINEAR KINEMATIC HARDENING - CANNOT BE USED FOR PLANE STRESS
C           By H. Badnava
C ----------------------------------------------------------------
C PROPS(1) - E
C PROPS(2) - NU
C PROPS(3) - SYIELD0
C PROPS(4) - Q
C PROPS(5) - b
C ----------------------------------------------------------------
C
C ELASTIC PROPERTIES
C
      EMOD=PROPS(1)
      ENU=MIN(PROPS(2), ENUMAX)
      EBULK3=EMOD/(ONE-TWO*ENU)
      EBULK = EBULK3/THREE
      EG2=EMOD/(ONE+ENU)
      EG=EG2/TWO
      EG3=THREE*EG
      ELAM=(EBULK3-EG2)/THREE
C
C ELASTIC STIFFNESS + DDSDDE FOR ELASTIC CONDITION
      DO K1=1, NDI
      DO K2=1, NDI
      DDSDDE(K2, K1)=ELAM
      END DO
      DDSDDE(K1, K1)=EG2+ELAM
      END DO
      DO K1=NDI+1, NTENS
      DDSDDE(K1, K1)=EG
      END DO
C
      DO K1=1, NTENS
      DO K2=1, NTENS      
      EE2= DDSDDE(K2, K1)
      END DO
      END DO
C      
C     RECOVER STATE VARIABLES
      P = STATEV(1)
      r = STATEV(2) 
C
C     CALCULATE PREDICTOR STRESS 
            DO K1=1, NTENS
            DO K2=1, NTENS
            STRESS(K2)=STRESS(K2)+DDSDDE(K2, K1)*DSTRAN(K1)
            END DO
            END DO
C            
C     CALCULATE DEVIATORIC TERIAL STRESS
            SHYDROSTATIC = 0.0D0
            DO K1=1, NDI
            SHYDROSTATIC=SHYDROSTATIC+STRESS(K1)
            END DO
            SHYDROSTATIC=SHYDROSTATIC/THREE
            DO K1=1,NDI
            STRILALDEV(K1)=STRESS(K1)-SHYDROSTATIC
            END DO
            DO K1=(NDI+1),NTENS
            STRILALDEV(K1)=STRESS(K1)
            END DO
C     
C     CALCULATE EQUIVALENT VON MISES TERIAL STRESS
      SMISESTER =0.D0
      DO K1=1,NDI
      SMISESTER = SMISESTER+STRILALDEV(K1)**2
      END DO
      DO K1=(NDI+1),NTENS
      SMISESTER=SMISESTER+2*(STRILALDEV(K1)**2)
      END DO
      SMISESTER=SQRT(SMISESTER*(THREE/TWO))
C
C     GET YIELD STRESS AND HARDENING MODULUS
C
      SYIELD0 = PROPS(3)
      Q=PROPS(4)
      b=PROPS(5)
C     SET EQUIVALENT STRAIN = ZERO      
      DP = ZERO
C
C     DETERMINE IF ACTIVELY YIELDING
C
      IF((SMISESTER-r).GT.(ONE+TOLER)*SYIELD0) THEN
C     ACTIVELY YIELDING
C
C     CALCULATE DELTAP -- NEWTON-RAPHSON METHOD
      r0=r
      DO 10 KTH=1,40      
      HARD = Q*b*EXP(-b*DP)
C     r = r+dr AND dr(p)= b(Q-r)dp            
      r = r0 + b*(Q-r)*DP
      DDP = ((SMISESTER-3*EG*DP-r-SYIELD0)/(3*EG+HARD))
      DP=DP+DDP
C     CHECK CONVERGENCE
      YF = SMISESTER-THREE*EG*DP-r-SYIELD0
            IF(DABS(YF).LT.TOLERANCE) GO TO 10
10    CONTINUE            
            P = P + DP
C
C     CALCULATE DELTA STRAIN PELASTIC
            DO K1=1,NTENS
      DELTASTRAINP(K1) = (THREE/TWO)*DP*(STRILALDEV(K1)/SMISESTER)
            END DO
C     CALCULATE ELASTIC STRAIN
            DO K1=1,NTENS
            DELTASTRAINE(K1) = DSTRAN(K1)- DELTASTRAINP(K1)
            END DO
C     UPDATE STRESS            
            DO K1=1, NTENS            
            STRESS(K1)=STRESS(K1)-TWO*EG*DELTASTRAINP(K1)            
            END DO
C
C     DETERMINE DEVIATORIC STRESS
            SHYDRO = ZERO
            DO K1=1,NDI
            SHYDRO=SHYDRO+(STRESS(K1)/THREE)
            END DO
            DO K1=1,NDI
            STRESSDEV(K1)=STRESS(K1)-SHYDRO
            END DO
            DO K1=(NDI+1),NTENS
            STRESSDEV(K1)=STRESS(K1)
            END DO            
C     
C     CALCULATE EFFECTIVE STRESS
            EFFSTRESS=0.D0
            DO K1=1,NDI
            EFFSTRESS = EFFSTRESS+STRESSDEV(K1)**2
            END DO
            DO K1=(NDI+1),NTENS
            EFFSTRESS=EFFSTRESS+2*(STRESSDEV(K1)**2)
            END DO
            EFFSTRESS=SQRT(EFFSTRESS*(THREE/TWO))
C      
C ***********************************************************
C ************  DETERMINE CONSISTENT TANGENT MODULUS ********
C ***********************************************************
      CC = EFFSTRESS/SMISESTER          
      BB = -(THREE/(TWO*(SMISESTER**2)))*(CC-(ONE/((3*EG/HARD)+1)))
C
C     DEFINE IDENTITY TENSOR
            DO K1=1,NTENS
            DO K2=1,NTENS
                IF (K1.EQ.K2)THEN
                IDENTITY(K1,K2)=ONE
                ELSE
                IDENTITY(K1,K2)=ZERO
                END IF
            END DO
            END DO
C
C     DEFINE IDENTITY VECTOR
            DO K1=1,NDI
            XIDENTITY(K1) = ONE
            END DO
            DO K1=(NDI+1),NTENS
            XIDENTITY(K1) = ZERO
            END DO
C
C     
            DO K1=1,NTENS
            DO K2=1,NTENS
            TERM1(K2,K1)=TWO*EG*IDENTITY(K2,K1)
            TERM2(K2,K1)=EBULK*XIDENTITY(K2)*XIDENTITY(K1)
            TERM3(K2,K1)=(TWO/THREE)*EG*CC*XIDENTITY(K2)*XIDENTITY(K1)
            TERM4(K2,K1)=EG2*BB* STRILALDEV(K2)*STRILALDEV(K1)
            END DO
            END DO
C
            DO K1=1,NTENS
            DO K2=1,NTENS
       DDSDDE(K2,K1)=TERM1(K2,K1)+TERM2(K2,K1)-TERM3(K2,K1)+TERM4(K2,K1)
            END DO
            END DO            
      END IF
C     SAVE STATE VARIABLES       
      STATEV(1)=P
      STATEV(2)=r
      RETURN
      END