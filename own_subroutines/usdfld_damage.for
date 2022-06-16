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
C Get PEEQ
C
      CALL GETVRM('PE',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)
      FIELD(1) = ARRAY(7)
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
      IF(STATEV(5).GT.0.04) THEN
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