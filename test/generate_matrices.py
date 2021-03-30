import numpy as np
import scipy.linalg as la

width = int(input())

while (True):
    M_1 = np.random.uniform(low=10.0, high=126.0, size=(width,width))
    M_2 = np.random.uniform(low=52.0, high=303.5, size=(width,width))
    print("Calculating product matrix...");
    P = np.matmul(M_1, M_2)
    Pinv = np.zeros((width, width));
    try:
        Pinv = la.inv(P)
    except LinAlgError:
        print("Ended up with a singular matrix. Trying again")
        continue
    except ValueError:
        print("The matrix is not a square 2D matrix");
        break;
    print("Calculating condition number...");
    condition_number = la.norm(P, ord='fro') * la.norm(Pinv, ord='fro');
    print(condition_number);
    # it turns out that matrices generated using
    # purely uniform distributions tend to be
    # poorly conditioned
    if (condition_number < 1.0e9):
        break;
    # else continue

np.savetxt('m_1.dat', M_1, delimiter=' ')
np.savetxt('m_2.dat', M_2, delimiter=' ')
np.savetxt('matmul.dat', P, delimiter=' ')
