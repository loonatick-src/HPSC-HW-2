import numpy as np
import scipy.linalg as la

widths = [100, 1000, 5000, 10000]

condition_number_threshold = 1.0e7;
for width in widths:
    while (True):
        print(f"Generating random matrices for width {width}")
        M_1 = np.random.uniform(low=10.0, high=126.0, size=(width,width))
        M_2 = np.random.uniform(low=52.0, high=303.5, size=(width,width))
        print("Calculating product matrix...");
        P = np.matmul(M_1, M_2)
        Pinv = np.zeros((width, width));
        try:
            Pinv = la.inv(P)
        except la.LinAlgError:
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
        if (condition_number < condition_number_threshold):
            break;
        else:
            print(f"Condition number larger than {condition_number_threshold}, trying again");
        # else continue

    condition_number_threshold *= 50.0;
    mfile_1 = "m_1_" + f"{width}" + ".dat"
    np.savetxt(mfile_1, M_1, delimiter=' ', fmt="%lf");
    mfile_2 = "m_2_" + f"{width}" + ".dat"
    np.savetxt(mfile_2, M_2, delimiter=' ', fmt="%lf");
    pfile = "matmul_" + f"{width}" + ".dat"
    np.savetxt(pfile, P, delimiter=' ', fmt="%lf");
    mfile_2t = "m_2_" + f"{width}t" + ".dat"
    np.savetxt(mfile_2t, M_2.T, delimiter=' ', fmt="%lf");
    pfilet = "matmul_" + f"{width}t" + ".dat"
    np.savetxt(pfilet, P.T, delimiter=' ', fmt="%lf");
