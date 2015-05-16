
*FILE GETPISEN.FOR
C ----------------------------------------------------------------
        SUBROUTINE GETPISEN(E,MM,RHO,MB,BB,P,D,T,U0,IFOK)
C ----------------------------------------------------------------
*  Description:
*     Iterative solver which is used to find density pressure
*     and total 4-velocity from energy, rest-mass, total
*     momentum density , and magnetic field B.
*    ( isentropic equation of state ) 
C ----------------------------------------------------------------
*  Input:
        REAL*8 E,MM,RHO,MB,BB
*         E  - energy density
*         MM - momentum density m[i]m[i]; i=1,3
*         RHO - mass density
*         MB  - scalar product m[i]B[i]
*         BB  - scalar product B[i]B[i]
C ----------------------------------------------------------------
*   Input/Output:
        REAL*8 P,D,T,U0
        LOGICAL IFOK
*         P   - pressure
*         D   - rest-mass density
*         U0   - Lorentz factor
*         It is assumed that when the subroutine is called U0 gives the 
*         initial guess (solution from the previous step). 
C ----------------------------------------------------------------
*   Commons:
C        INCLUDE 'FLOOR.CMN'
        INCLUDE 'STATE.CMN'
C ----------------------------------------------------------------
*   Constants:
        REAL*8 ALFA
        PARAMETER (ALFA=1.E-10) ! precision parameter
        INTEGER CLIM
        PARAMETER (CLIM=30)  ! maximum number of iterations
       
C ----------------------------------------------------------------
*   Variables:
        INTEGER I
        REAL*8 LOR1,LOR2,F,DFDX
C ----------------------------------------------------------------
                 
        IFOK =.TRUE.        

        LOR1=U0
        DO I=1,CLIM+1    ! Newton iterations
             IF(I.EQ.CLIM)THEN 
                 IFOK=.FALSE.
                 RETURN
             ENDIF
             CALL FDFDX(MM,RHO,MB,BB,LOR1,F,DFDX)
             LOR2=LOR1-F/DFDX
C             print *,i,lor1,F,DFDX
             IF(ABS(LOR2-LOR1).LT.ALFA) GOTO 10
             LOR1=LOR2
         ENDDO

10       U0=LOR2
         IF(U0.LT.1.)IFOK=.FALSE.
         D=RHO/U0
         P=KAPPA_I*D**GAM
         T=P/D

         END 

* -------------------------------------------------------------
         SUBROUTINE FDFDX(MM,RHO,MB,BB,LOR,FX,DFDX)
* -------------------------------------------------------------
*   Input:
      REAL*8 MM,RHO,MB,BB,LOR
*   Output:
      REAL*8 FX,DFDX
*   Commons:
      INCLUDE 'STATE.CMN'
*   Variables: 
      REAL*8 A,DF,DA,V2,F,F2,F3,F4,P,E,B4,MB2 

      B4=BB**2
      MB2=MB**2
      V2=(LOR**2-1.D0)/LOR**2
      P=KAPPA_I*RHO**GAM
      E=(GAM/(GAM-1.D0))*P
      F=RHO*LOR + E*LOR**(2.D0-GAM)
      F2=F*F
      F3=F2*F
      F4=F2*F2
      A=V2*(F4+B4*F2+2.D0*BB*F3)
      FX=BB*MB2 + MM*F2+2.D0*MB2*F - A
      DF=RHO+E*(2.D0-GAM)*LOR**(1.D0-GAM)
      DA = (F4+B4*F2+2.D0*BB*F3)*2.D0/LOR**3
      DA = DA+(4.D0*F3+2.D0*B4*F+6.D0*BB*F2)*V2*DF 
      DFDX = (2.D0*MM*F+2.D0*MB2)*DF-DA

      END

! STATE.CMN
!        INTEGER IEOS
!        COMMON /EOS_TYPE/ IEOS
!        REAL*8 GAM, KAPPA_I
!        COMMON /ISEN_POLY/ GAM,KAPPA_I
