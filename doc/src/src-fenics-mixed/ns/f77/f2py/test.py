import nose.tools as nt

def test_bcmodelf77():
    import bcmodelf77
    C = 0.127
    R_d = 5.43
    N = 2
    P_1 = 16000
    Q = 1000
    T = 0.01
    P_ = bcmodelf77.pmodel(P_1, R_d, Q, C, N, T)

    # Manual formula:
    P1_ = P_1 + T/2*(Q - P_1/R_d)/C
    P1_ = P1_ + T/2*(Q - P1_/R_d)/C
    nt.assert_almost_equal(
        P_[-1], P1_, places=10, msg='F77: %g, manual coding: %s' %
        (P_[-1], P1_))

if __name__ == '__main__':
    test_bcmodelf77()

