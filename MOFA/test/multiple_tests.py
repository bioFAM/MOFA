from test_funct import test
import numpy as np
import sys


def multiple_tests(res_dir, to_test='K', val=None,  it_ix=0):
    if val is None:
        if to_test =='K':
            #define default values
            pass
        elif to_test =='D':
            pass
        elif to_test =='M':
            pass
        elif to_test =='N':
            pass

    elif type(val) == float or type(val) == int:
        val = np.array([val])

    # varry K
    if to_test == 'K':
        time_K = np.zeros(len(val))
        for i in range(len(val)):
            time_K[i] = test(K=val[i])
        import pdb; pdb.set_trace()
        res = np.concatenate((val[:,None], time_K[:,None]), axis=1)

    # varry D
    elif to_test == 'D':
        time_D = np.zeros(len(val))
        for i in range(len(val)):
            time_D[i] = test(D=val[i])
        res = np.concatenate((val[:,None], time_D[:,None]), axis=1)

    # varry M
    elif to_test == 'M':
        time_M = np.zeros(len(val))
        for i in range(len(val)):
            time_M[i] = test(M=val[i])
        res = np.concatenate((val[:,None], time_M[:,None]), axis=1)

    # varry N
    elif to_test == 'N':
        time_N = np.zeros(len(val))
        for i in range(len(val)):
            time_N[i] = test(N=val[i])
        res = np.concatenate((val[:,None], time_N[:,None]), axis=1)

    file_name = res_dir + '/' + to_test + '_' + str(it_ix)
    with open(file_name, 'w') as f:
        np.savetxt(f,
                   res,
                   delimiter=' ')




if __name__ == '__main__':
    to_test = sys.argv[1]
    res_dir = sys.argv[2]
    val = int(sys.argv[3])
    it = int(sys.argv[4])

    multiple_tests(res_dir=res_dir, to_test=to_test, val=val, it_ix=it)
