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
C
C Absolute value of current strain:
      CALL GETVRM('E',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)
      EPS = ABS( ARRAY(1) )
C Maximum value of strain up to this point in time:
      CALL GETVRM('SDV',ARRAY,JARRAY,FLGRAY,JRCD,JMAC,JMATYP,
     1 MATLAYO,LACCFLA)
      EPSMAX = ARRAY(1)
C Use the maximum strain as a field variable
      FIELD(1) = MAX( EPS , EPSMAX )
C Store the maximum strain as a solution dependent state 
C variable
      STATEV(1) = FIELD(1)
C If error, write comment to .DAT file:
      IF(JRCD.NE.0)THEN
       WRITE(6,*) 'REQUEST ERROR IN USDFLD FOR ELEMENT NUMBER ',
     1     NOEL,'INTEGRATION POINT NUMBER ',NPT
      ENDIF
C
      RETURN
      END
