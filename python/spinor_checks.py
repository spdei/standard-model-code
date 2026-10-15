"""Other representations: the spinor of so(14) minus its highest weight gives d_Z = 3, and the 27
of E6 is free of cubic anomalies and gives d_Z = 2 whichever weight is removed."""
import numpy as np

from sm16 import spinor, z_distance


def cubic_free(W):
    A = np.array(W)
    return not np.einsum('li,lj,lk->ijk', A, A, A).any()


def main():
    W7 = spinor(7)
    assert cubic_free(W7)
    s7 = W7.index((1,) * 7)
    for p in (3, 5, 7):
        assert z_distance(W7, (s7,), p) == 3
    print('spinor of so(14): no cubic anomaly; removing the highest weight gives [[63,1,3]]_p '
          'for p = 3, 5, 7')
    # 27 of E6 = 16_1 + 10_{-2} + 1_4 under so(10) x u(1), in integral coordinates
    W27 = [tuple(q) + (1,) for q in spinor(5)] \
        + [tuple(2 * sg if j == i else 0 for j in range(5)) + (-2,) for i in range(5) for sg in (1, -1)] \
        + [(0, 0, 0, 0, 0, 4)]
    assert len(set(W27)) == 27 and cubic_free(W27)
    for p in (3, 5, 7, 11):
        for s in range(27):
            assert z_distance(W27, (s,), p) == 2
    print('27 of E6: no cubic anomaly; removing any weight gives d_Z = 2 for p = 3, 5, 7, 11')


if __name__ == '__main__':
    main()
