import bcmodelf77

class BCModel:
    def __init__(self, C, R_d, N, T):
        self.C, self.R, self.N, self.T = C, R_d, N, T

    def __call__(self, P, Q):
        P_ = bcmodelf77.pmodel(
             P, self.R, Q, self.C, self.N, self.T)
        return P_

import nose.tools as nt

def test_BCModel():
    C = 0.127; R_d = 5.43; N = 2; T = 0.01
    pmodel = BCModel(C, R_d, N, T)
    P_1 = 16000
    Q = 1000
    P_ = pmodel(P_1, Q)
    # Manual formula:
    P1_ = P_1 + T/2*(Q - P_1/R_d)/C
    P1_ = P1_ + T/2*(Q - P1_/R_d)/C
    nt.assert_almost_equal(
        P_[-1], P1_, places=10, msg='F77: %g, manual coding: %s' %
        (P_[-1], P1_))

if __name__ == "__main__":
  test_BCModel()






