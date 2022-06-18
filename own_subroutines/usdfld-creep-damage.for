      SUBROUTINE USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,
     3 LACCFLA)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME,ORNAME
      CHARACTER*3  FLGRAY(15)
      DIMENSION FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      DIMENSION ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),
     1 COORD(*)
C------------------------------------------------------------
C STATEV(1)--PEEQ    | STATEV(2)--Damage flag
C STATEV(3)--DamageT | STATEV(4)--PEEQ0(damage initiation)
C STATEV(5)--PEEQ-PEEQ0
C------------------------------------------------------------
C Get PEEQ and CEEQ
C
      CALL GETVRM('PE',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)
      FIELD(1) = ARRAY(7)
      CALL GETVRM('CE',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)
      FIELD(1) = FIELD(1)+ARRAY(7)
C
C Store the PEEQ as a solution dependent state variable
C
      STATEV(1) = FIELD(1)
      STATEV(3) = STATEV(1)/0.15
C
C Before damage initiation, store the PEEQ
C
      IF(STATEV(3).LT.1) THEN
        STATEV(4)=STATEV(1)
C
C Damage initiation
C
      ELSE
        STATEV(3)=1
      END IF
C
C Damage evolution
C
      IF(STATEV(3).NE.1) THEN
        GOTO 10
      END IF
      STATEV(5)=STATEV(1)-STATEV(4)
      IF(STATEV(5).GT.0.02) THEN
        IF(STATEV(2).NE.0) THEN
          STATEV(2)=0
          WRITE(16,*) 'ELEMENT',NOEL,'DELETED AFTER INCREMENT'
     1                 ,KINC,'.'
        END IF
      END IF
C
C If error, write comment to .DAT file:
C
10    IF(JRCD.NE.0)THEN
       WRITE(6,*) 'REQUEST ERROR IN USDFLD FOR ELEMENT NUMBER ',
     1     NOEL,'INTEGRATION POINT NUMBER ',NPT
      ENDIF
C
      RETURN
      END

      SUBROUTINE CREEP(DECRA,DESWA,STATEV,SERD,EC,ESW,P,QTILD,
     1 TEMP,DTEMP,PREDEF,DPRED,TIME,DTIME,CMNAME,LEXIMP,LEND,
     2 COORDS,NSTATV,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
C
      INCLUDE 'ABA_PARAM.INC'
C
      CHARACTER*80 CMNAME
C
      DIMENSION DECRA(5),DESWA(5),STATEV(*),PREDEF(*),DPRED(*),
     1 TIME(3),COORDS(*),EC(2),ESW(2)
C
C DEFINE CONSTANTS
C
C      A=2.9e-29
      A =2.9E-22
      AN=8
      AM=0
C
      DECRA(1) = A*(QTILD**AN)*(TIME(1)**AM)*DTIME
      IF(LEXIMP.EQ.1) THEN
       DECRA(5) = AN*A*(QTILD**(AN-1.))*DTIME*
     1             (TIME(1)**AM)
      END IF
C
      RETURN
      END

