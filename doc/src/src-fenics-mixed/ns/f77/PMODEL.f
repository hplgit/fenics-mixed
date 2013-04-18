      SUBROUTINE PMODEL(P, P_1, R_D, Q, C, N, T)
C     Integrate P in N steps from 0 to T, given start value P_1
      INTEGER N
      REAL*8 P(0:N), P_1, R_D, Q, C, T
      REAL*8 DT
      INTEGER I
Cf2py intent(in) P0, R, Q, C, N
Cf2py intent(out) P

      DT = T/N
      P(0) = P_1
      DO I = 0, N-1
        P(i+1) = P(i) + DT*(Q - P(i)/R_D)/C
      END DO
      END SUBROUTINE PMODEL
