C     4X4 JACOBIAN THEN CONDENSATION WITH REGULARIZATION 
c     (energy due to viscous regularization is calculated)
      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1     RPL,DDSDDT,DRPLDE,DRPLDT,
     2     STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3     NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4     CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C     
      INCLUDE 'ABA_PARAM.INC'
C     
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1     DDSDDE(NTENS,NTENS),
     2     DDSDDT(NTENS),DRPLDE(NTENS),
     3     STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     4     PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)
      
      DIMENSION STRANT(6),TSTRANT(4)
      DIMENSION CFULL(6,6),CDFULL(6,6)
      DIMENSION DDFDE(6), DDMDE(6), DCDDF(6,6), DCDDM(6,6)
      DIMENSION ATEMP1(6), ATEMP2(6), TDDSDDE(6,6)
      DIMENSION OLD_STRESS(6)
      DIMENSION DOLD_STRESS(6),D_STRESS(6)
      PARAMETER (ZERO = 0.D0,ONE = 1.D0,TWO = 2.D0, HALF = 0.5D0)
C****************************
C     STRANT..... STRAIN AT THE END OF THE INCREMENT
C     TSTRANT.....TEMPORARY ARRAY TO HOLD THE STRAIN FOR PLANE STRESS PROBLEM
C     CFULL.......FULL 6X6 ELASTICITY MATRIX
C     CDFULL......FULL 6X6 DAMAGED ELASTICITY MATRIX
C     DDFDE....... D DF/D E
C     DDMDE....... D DM/D E
C     DCDDF....... D C/ D DF THE DERIVATIVE OF THE FULL MATRIX OVER DF
C     DCDDM........D C/ D DM THE DERIVATIVE OF THE FULL MATRIX OVER DM
C     ATEMP1,ATEMP2...TEMPORARY ARRAY USED IN JACOBIAN CALCULATION
C     TDDSDDE.....UNCONDENSED JACOBIAN MATRIX FOR PLANE STRESS PROBLEM
C     OLD_STRESS...STRESS AT THE BEGINNING OF THE INCREMENT, SAVED FOR THE ENERGY
C                  COMPUTATION
C     DOLD_STRESS...STRESS AT THE BEGINNING OF THE INCREMENT, 
C                  IF THERE'S NO VISCOUS REGULARIZATION
C     D_STRESS...STRESS IF THERE'S NO VISCOUS REGULARIZATION, THE ABOVE IS CALCULATED
C                TO CALCULATE THE SCD, ENERGY CAUSED BY VISCOUS REGULARIZATION
C     STATEV(1)   damage variable df
C     STATEV(2)   damage variable dm
C     STATEV(3)   regularized damage variable dfv
C     STATEV(4)   regularizaed damage variable dmv
C     STATEV(5:10) TEMPORARY ARRAYS TO SAVE DOLD_STRESS
C************
C
C     GET THE MATERIAL PROPERTIES---ENGINEERING CONSTANTS
C
      TENL = PROPS(1)           !YOUNG'S MODULUS IN DIRECTION 1 (L)
      TENT = PROPS(2)           !YOUNG'S MODULUS IN DIRECTION 2 (T)
      SHRLT = PROPS(3)          !SHEAR MODULUS IN 12 PLANE
      SHRTT = PROPS(4)          !SHEAR MODULUS IN 23 PLANE
      XNULT = PROPS(5)          !POISON'S RATIO POI_12
      XNUTT = PROPS(6)          !POISON'S RATIO POI_23
      XNUTL = XNULT / TENL * TENT !POI_21
C     
C     GET THE FAILURE PROPERTIES
C
      SIGTL = PROPS(7)          !FAILURE STRESS IN 1 DIRECTION IN TENSION
      SIGCL = PROPS(8)          !FAILURE STRESS IN 1 DIRECTION IN COMPRESSION
      SIGTT = PROPS(9)          !FAILURE STRESS IN 2 DIRECTION IN TENSION
      SIGCT = PROPS(10)          !FAILURE STRESS IN 2 DIRECTION IN COMPRESSION
      SIGSLT = PROPS(11)        !FAILURE STRESS IN SHEAR IN 1-2 PLANE
      GFMAT = PROPS(12)         !FRACTURE ENERGY IN MATRIX
      GFFIB = PROPS(13)         !FRACTURE ENERGY IN FIBER
      ETA = PROPS(14)           ! VISCOSITY FOR REGULARIZATION
C     
C     CALCULATE THE STRAIN AT THE END OF THE INCREMENT
C     
      DO I = 1, NTENS
         STRANT(I) = STRAN(I) + DSTRAN(I)
      END DO
c     
C     FILL THE 6X6 FULL STIFFNESS MATRIX
      DO I = 1, 6
         DO J = 1, 6
            CFULL(I,J)=ZERO
         END DO
      END DO
      ATEMP = ONE - TWO * XNULT * XNUTL - XNUTT ** TWO
     1     - TWO * XNULT * XNUTL * XNUTT
      CFULL(1,1) = TENL * (ONE - XNUTT ** TWO) / ATEMP
      CFULL(2,2) = TENT * (ONE - XNULT * XNUTL) / ATEMP
      CFULL(3,3) = CFULL(2,2)
      CFULL(1,2) = TENT * (XNULT + XNULT * XNUTT) / ATEMP
      CFULL(1,3) = CFULL(1,2)
      CFULL(2,3) = TENT * (XNUTT + XNULT * XNUTL) / ATEMP
      CFULL(4,4) = SHRLT
      CFULL(5,5) = SHRLT
      CFULL(6,6) = SHRTT
      DO I = 2, 6
         DO J = 1, I-1
            CFULL(I,J) = CFULL(J,I)
         END DO
      END DO
c calculate the failure strain by failure stress
      EPITL = SIGTL / cfull(1,1) !FAILURE STRAIN 1 DIRECTION IN TENSION
      EPICL = SIGCL / cfull(1,1) !FAILURE STRAIN 1 DIRECTION IN COMPRESSION
      EPITT = SIGTT / cfull(2,2) !TENSILE FAILURE STRAIN 2 DIRECTION
      EPICT = SIGCT / cfull(2,2) !COMPRESSIVE FAILURE STRAIN 2 DIRECTION
      EPISLT = SIGSLT/ SHRLT    ! FAILURE SHEAR STRAIN ...ENGINEERING STRAIN
C     
C     CHECK THE FAILURE INITIATION CONDITION
c     
      DFOLD = STATEV(1)
      DMOLD = STATEV(2)
      DFVOLD = STATEV(3)
      DMVOLD = STATEV(4)
      CALL CheckFailureIni(EPITL,EPICL,EPITT,EPICT,EPISLT,STRANT,
     1     GFMAT,GFFIB, CELENT, CFULL, DF, DM, DDFDE, DDMDE, NTENS,
     2  DFOLD, DMOLD,NDI)
C     
C     ! USE VISCOUS REGULARIZATION
C     
      DFV = ETA / (ETA + DTIME) * DFVOLD + DTIME / (ETA + DTIME) * DF
      DMV = ETA / (ETA + DTIME) * DMVOLD + DTIME / (ETA + DTIME) * DM
C     SAVE THE OLD STRESS TO OLD_STRESS
      DO I = 1, NTENS
         OLD_STRESS(I) = STRESS(I)
      END DO

C     CALL ROUTINE TO CALCULATE THE STRESS
C     CALCULATE THE STRESS IF THERE'S NO VISCOUS REGULARIZATION
      CALL GetStress(CFULL,CDFULL,DF,DM,D_STRESS,STRANT,NDI,NTENS)
C     CALCULATE THE STRESS IF THERE'S VISCOUS REGULARIZATION
      CALL GetStress(CFULL,CDFULL,DFV,DMV,STRESS,STRANT,NDI,NTENS)
C     GET THE OLD STRESS IF THERE'S NO VISCOUS REGULARIZATION
      DO I=1,NTENS
         DOLD_STRESS(I)=STATEV(I+4)
      END DO
C     SAVE THE CURRENT STRESS IF THERE'S NO VISCOUS REGULARIZATION
      DO I=1,NTENS
         STATEV(I+4)=D_STRESS(I)
      END DO
C     
C     CALCULATE THE DERIVATIVE MATRIX DC/DDM, DC/DDF OF THE DAMAGED MATRIX
C     
      CALL ElasticDerivative(CFULL,DMV,DFV, DCDDM,DCDDF)
C     
C     UPDATE THE JACOBIAN
C     
C     FULL 3D CASE
      IF (NDI .EQ. 3) THEN
         DO I = 1, NTENS
            ATEMP1(I) = ZERO
            DO J = 1, NTENS
               ATEMP1(I) = ATEMP1(I) + DCDDM(I,J) * STRANT(J)
            END DO
         END DO
         
         DO I = 1, NTENS
            ATEMP2(I) = ZERO
            DO J = 1, NTENS
               ATEMP2(I) = ATEMP2(I) + DCDDF(I,J) * STRANT(J)
            END DO
         END DO
         
         DO I = 1, NTENS
            DO J = 1, NTENS
               DDSDDE(I,J)=CDFULL(I,J) + ( ATEMP1(I) * DDMDE(J)
     1              + ATEMP2(I) * DDFDE(J) ) * DTIME / (DTIME + ETA)
            END DO
         END DO
C     
C     ! PLANE STRESS CASE
C     
      ELSE IF (NDI .EQ.2) THEN
         TSTRANT(1) = STRANT(1)
         TSTRANT(2) = STRANT(2)
         TSTRANT(3) = -CDFULL(1,3) / CDFULL(3,3) * STRANT(1)
     1        - CDFULL(2,3) / CDFULL(3,3) * STRANT(2)
         TSTRANT(4) = STRANT(3)
         DO I = 1, 4
            ATEMP1(I) = ZERO
            DO J = 1, 4
               ATEMP1(I) = ATEMP1(I) + DCDDM(I,J) * TSTRANT(J)
            END DO
         END DO
         
         DO I = 1, 4
            ATEMP2(I) = ZERO
            DO J = 1, 4
               ATEMP2(I) = ATEMP2(I) + DCDDF(I,J) * TSTRANT(J)
            END DO
         END DO
         DO I = 1,6
            DO J = 1,6
            TDDSDDE(I,J) = ZERO
            END DO
         END DO
C     TO GET THE UNCONDENSED JACOBIAN FOR PLANE STRESS CASE
         DO I = 1, NTENS
            DO J = 1, NTENS
               DDSDDE(I,J) = ZERO
            END DO
         END DO
         DO I = 1, 4
            DO J = 1, 4
               TDDSDDE(I,J)=CDFULL(I,J) + ( ATEMP1(I) * DDMDE(J)
     1              + ATEMP2(I) * DDFDE(J) ) * DTIME / (DTIME + ETA)
            END DO
         END DO
C     
C     CONDENSE THE JACOBIAN MATRIX FOR PLANE STRESS PROBLEM
C     
         CALL MatrixCondense(TDDSDDE,DDSDDE)
      END IF 
C     
C     TO UPDATE THE STATE VARIABLE
C     
      STATEV(1) = DF
      STATEV(2) = DM
      STATEV(3) = DFV
      STATEV(4) = DMV
C     
C     TO COMPUTE THE ENERGY
C     
      DO I = 1, NDI
         SSE = SSE + HALF * (STRESS(I) + OLD_STRESS(I)) * DSTRAN(I)
      END DO
      DO I = NDI+1, NTENS
         SSE = SSE + (STRESS(I) + OLD_STRESS(I)) * DSTRAN(I)
      END DO
C     TO COMPUTE THE INTERNAL ENERGY WITHOUT VISCOUS REGULARIZATION
      DO I = 1, NDI
         SCD = SCD + HALF * (STRESS(I) + OLD_STRESS(I)
     1        -D_STRESS(I)-DOLD_STRESS(I)) * DSTRAN(I)
      END DO
      DO I = NDI+1, NTENS
         SCD = SCD + (STRESS(I) + OLD_STRESS(I)
     1        -D_STRESS(I)-DOLD_STRESS(I)) * DSTRAN(I)
      END DO
      
      RETURN
      END
      
C******************************************************************************
C CALCULATE THE STRESS BASED ON THE DAMAGE VARAIBLES***************************
C******************************************************************************
      SUBROUTINE GetStress(CFULL,CDFULL,DFV,DMV,STRESS,STRANT,NDI,NTENS)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION CFULL(6,6),CDFULL(6,6),STRESS(NTENS),
     1     STRANT(6),CDTHREE(3,3)
      PARAMETER (ZERO=0.D0, ONE=1.D0)
C     CDTHREE.....DAMAGED CONDENSED-ELASTICITY MATRIX FOR PLANE STRESS PROBLEM
      DO I = 1, 6
         DO J = 1, 6
            CDFULL(I,J)=CFULL(I,J)
         END DO
      END DO
      IF ( (DFV .NE. ZERO) .OR. (DMV .NE. ZERO)) THEN
         CDFULL(1,1) = (ONE - DFV) * CFULL(1,1)
         CDFULL(1,2) = (ONE - DFV) * (ONE - DMV) * CFULL(1,2)
         CDFULL(2,1) = CDFULL(1,2)
         CDFULL(2,2) = (ONE - DMV) * CFULL(2,2)
         CDFULL(1,3) = (ONE - DFV) * CFULL(1,3)
         CDFULL(3,1) = CDFULL(1,3)
         CDFULL(2,3) = (ONE- DMV) * CFULL(2,3)
         CDFULL(3,2) = CDFULL(2,3)
         CDFULL(4,4) = (ONE - DMV) * (ONE - DFV) * CFULL(4,4)
      END IF
C   UPDATE THE STRESS STATE IF 3D CASE
C
      IF (NDI .EQ. 3) THEN
         DO I = 1, NTENS
            STRESS(I)=ZERO
            DO J = 1, NTENS
               STRESS(I)=STRESS(I)+CDFULL(I,J) * STRANT(J)
            END DO
         END DO
         
C     
C     INITIALIZE THE 3X3 CONDENSED STIFFNESS MATRIX IF PLANE STRESS CASE
C     
      ELSE IF ( NDI .EQ. 2) THEN
         DO I = 1, NTENS
            DO J = 1, NTENS
               CDTHREE(I,J)=ZERO
            END DO
         END DO
C     
C     
C     CONDENSE THE UNDAMAGED STIFFNESS MATRIX
C     
         CALL MatrixCondense(CDFULL,CDTHREE)
C     
C     UPDATE THE STRESS
C     
         DO I = 1, NTENS
            STRESS(I)=ZERO
            DO J = 1, NTENS
               STRESS(I)=STRESS(I)+CDTHREE(I,J) * STRANT(J)
            END DO
         END DO
      END IF 
      RETURN
      END
C******************************************************************************
C     TO CHECK THE FAILURE INITIATION AND THE CORRESPONDING DERIVATIVE*********
C******************************************************************************
      SUBROUTINE CheckFailureIni(EPITL,EPICL,EPITT,EPICT,EPISLT,STRANT,
     1     GFMAT,GFFIB, CELENT, CFULL, DF, DM, DDFDE, DDMDE, NTENS,
     2     DFOLD,DMOLD, NDI )
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION DDFDE(6), DDMDE(6), STRANT(6), CFULL(6,6)
      DIMENSION DFMNDE(6), DFFNDE(6)
      PARAMETER (ZERO = 0.D0, ONE = 1.D0, TWO = 2.D0, HALF = 0.5D0)
C     
C     CHECK THE INITIATION CONDITION FOR MATRIX
C     FMN=FM/EPITT > 1 THEN EVALUATE THE DAMAGE VARIABLE AND DERIVATIVE
C     
      TERM1 = STRANT(2)**TWO / EPICT / EPITT
      TERM2 = (EPICT - EPITT) / EPICT / EPITT * STRANT(2)
      IF (NDI .EQ. 3) THEN
         TERM3 = (STRANT(4))**TWO / EPISLT**TWO
      ELSE IF (NDI .EQ. 2) THEN
         TERM3 = (STRANT(3))**TWO / EPISLT**TWO
      END IF
      TERM = TERM1 + TERM2 + TERM3
      IF (TERM .GT. ZERO) THEN
         FMN = SQRT(TERM)
      ELSE
         FMN = ZERO
      END IF
C
C     INITIALIZE THE ARRAY AND VARIABLE
C
      DM = ZERO
      DDMDFMN = ZERO
      DO I = 1, 6
         DFMNDE(I) = ZERO
         DDMDE(I) = ZERO
      END DO
      IF (FMN .GT. ONE) THEN
C     CALCULATE DM, DDMDFMN
         CALL DamageEvaluation( CFULL(2,2), FMN, GFMAT, CELENT,
     1        EPITT, DM, DDMDFMN)
C     CALCULATE DFMNDE
         IF (DM .GT. DMOLD) THEN
            DFMNDE(2) = HALF / FMN * (TWO * STRANT(2) + EPICT - EPITT)
     1           / EPICT / EPITT
            IF (NDI .EQ. 3) THEN
               DFMNDE(4) = ONE / FMN * STRANT(4) / EPISLT**TWO
            ELSE IF (NDI .EQ. 2) THEN
               DFMNDE(4) = ONE / FMN * STRANT(3) / EPISLT**TWO
            END IF
            DO I = 1, 6
               DDMDE(I) = DFMNDE(I) * DDMDFMN
            END DO
         END IF
      END IF
      DM = MAX (DM, DMOLD)
C     
C     CHECK THE INITIATION CONDITION FOR FIBER
C     FFN=FF/EPITL>1 THEN CALCULATE THE DAMAGE VARIABLE AND DERIVATIVE
C     
	TERM1 = STRANT(1)**TWO / EPICL / EPITL
	TERM2 = (EPICL - EPITL) / EPICL / EPITL *STRANT(1)
        TERM = TERM1 + TERM2
        IF (TERM .GT. ZERO) THEN
           FFN = SQRT(TERM)
        ELSE
           FFN = ZERO
        END IF
	DF = ZERO
	DDFDFFN = ZERO
	DO I = 1, 6
           DFFNDE(I) = ZERO
           DDFDE(I) = ZERO
	END DO
	IF (FFN .GT. ONE) THEN
C     CALCULATE DF, DDFDFFN
           CALL DamageEvaluation( CFULL(1,1), FFN, GFFIB, CELENT,
     1          EPITL, DF, DDFDFFN)
C     CALCULATE DFFNDE
           IF (DF .GT. DFOLD) THEN
              DFFNDE(1) = HALF / FFN * (TWO * STRANT(1) + EPICL - EPITL)
     1             /  EPICL / EPITL
              DDFDE(1) = DFFNDE(1) * DDFDFFN
           END IF
	END IF
        DF = MAX (DF, DFOLD)
	RETURN
	END
C******************************************************************************
C     SUBROUTINE TO EVALUATE THE DAMAGE AND THE
c     DERIVATIVE************************
C******************************************************************************
      SUBROUTINE DamageEvaluation(STIFF, FN, GF, CELENT, EPIT, D,
     1     DDDFN)
C     CALCULATE DAMAGE VARIABLE
      INCLUDE 'ABA_PARAM.INC'
      PARAMETER (ONE = 1.D0,tol=1d-3, zero = 0.d0)
      TERM1 = STIFF * EPIT**2 * CELENT / GF
      TERM2 = (ONE - FN) * TERM1
      D = ONE - EXP(TERM2) / FN
C     CALCULATE THE DERIVATIVE OF DAMAGE VARIABLE WITH RESPECT TO FAILURE
C     RITERION
      DDDFN = (ONE / FN + TERM1) * (ONE - D)
      RETURN
      END
C******************************************************************************
C     SUBROUTINE TO CONDENSE THE 4X4 MATRIX INTO 3X3 MATRIX********************
C******************************************************************************
      SUBROUTINE MatrixCondense(CFULL,CTHREE)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION CFULL(6,6),CTHREE(3,3)
C     
      CTHREE(1,1) = CFULL(1,1) - CFULL(1,3) * CFULL(3,1) / CFULL(3,3)
      CTHREE(1,2) = CFULL(1,2) - CFULL(1,3) * CFULL(3,2) / CFULL(3,3)
      CTHREE(2,1) = CFULL(2,1) - CFULL(2,3) * CFULL(3,1) / CFULL(3,3)
      CTHREE(2,2) = CFULL(2,2) - CFULL(2,3) * CFULL(3,2) / CFULL(3,3)
      CTHREE(3,3) = CFULL(4,4)
      RETURN
      END
C*******************************************************************************
C     SUBROUTINE TO GET THE DERIVATIVE MATRIX OF CONDENSE DAMAGED MATRIX OVER
C**** THE DAMAGE VARIABLE******************************************************
C*******************************************************************************
      SUBROUTINE ElasticDerivative(CFULL,DMV,DFV, DCDDM,DCDDF)
      INCLUDE 'ABA_PARAM.INC'
      DIMENSION CFULL(6,6), DCDDM(6,6),
     1     DCDDF(6,6)
      PARAMETER (ZERO = 0.D0, ONE = 1.D0, HALF = 0.5D0)
C     initialize the data to zero
      DO I = 1, 6
         DO J = 1, 6
            DCDDM(I,J) = ZERO
            DCDDF(I,J) = ZERO
         END DO
      END DO
C     
C     CALCULATE DC/DDF
C     
      DCDDF(1,1) = - CFULL(1,1)
      DCDDF(1,2) = - (ONE - DMV) * CFULL(1,2)
      DCDDF(2,1) = DCDDF(1,2)
      DCDDF(1,3) = -CFULL(1,3)
      DCDDF(3,1) = DCDDF(1,3)
      DCDDF(4,4) = -(ONE - DMV) * CFULL(4,4)
C     
C     CALCULATE DC/DDM
C     
      DCDDM(1,2) = - (ONE - DFV) * CFULL(1,2)
      DCDDM(2,1) = DCDDM(1,2)
      DCDDM(2,2) = -CFULL(2,2)
      DCDDM(2,3) = -CFULL(2,3)
      DCDDM(3,2) = DCDDM(2,3)
      DCDDM(4,4) = -(ONE - DFV) * CFULL(4,4)
      RETURN
      END
      
