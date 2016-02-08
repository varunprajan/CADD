C
C
C     User Defined Field Subroutine (USDFLD)
C     Written by Varun P. Rajan, July 2015
C
      subroutine USDFLD(FIELD,STATEV,PNEWDT,DIRECT,T,CELENT,
     1 TIME,DTIME,CMNAME,ORNAME,NFIELD,NSTATV,NOEL,NPT,LAYER,
     2 KSPT,KSTEP,KINC,NDI,NSHR,COORD,JMAC,JMATYP,MATLAYO,LACCFLA)
C
C     include 'ABA_PARAM.INC'
C      
      character*80 CMNAME,ORNAME
      character*3 FLGRAY(15)
      dimension FIELD(NFIELD),STATEV(NSTATV),DIRECT(3,3),
     1 T(3,3),TIME(2)
      dimension ARRAY(15),JARRAY(15),JMAC(*),JMATYP(*),COORD(*)
C
      real x, y, z
      real pos1, pos2
C      
      x = COORD(1)
      y = COORD(2)
      z = COORD(3)
C     
C     Elastic modulus is assumed to vary linearly with the (1st) field
C     variable. Field variable ranges from 0 (Emin) to 1 (Emax).
C     Here, I have implemented the variation for a material with a linear
C     modulus gradient in the y-direction where E = Emin at y = pos1
C     and E = Emax at y = pos2, and all y-values are between pos1 and pos2
C      
      FIELD(1) = (y - pos1)/(pos2 - pos1)
C      
      return
      end