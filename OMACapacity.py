import cvxpy as cvx
import matplotlib.pyplot as plot
import numpy as np
import spectrum as sp

def opt_allocatin(w, h, P_max, B_max):
    N = len(w)
    P = cvx.Variable(N)
    B = cvx.Variable(N)
    w = np.array(w)
    h = np.array(h)
    r = -cvx.multiply(w, cvx.kl_div(B, B + cvx.multiply(P, h)) - cvx.multiply(P, h));
    obj = cvx.Maximize(cvx.sum(r))
    con = [P >= 0.0, B >= 0.0, cvx.sum(P) <= P_max, cvx.sum(B) <= B_max]    
    prob = cvx.Problem(obj, con)
    prob.solve()
    return prob.status, prob.value, P.value, B.value

def uniform_allocation(h, P_max, B_max):
    N = len(h);
    R = B_max / N * np.log2(1 + P_max * h);
    return np.sum(R)


def main1():
    K = 5  # Cell number
    T = 100  #
    radii = np.array([2, 2, 2, 2, 2])  # Cell radii
    lambs = np.array([0.5, 0.5, 0.5, 0.5, 0.5])  # Cell user achieve rate
    P_max = np.ones(K)
    B_max = np.ones(K)
    means = lambs * np.pi * radii ** 2
    dis = []  # Distance of j-th user in i-th cell dis[i-1][j-1] 
    cha = []  # CSI of j-th user in i-th cell dis[i-1][j-1] 
    gain = []  # Channel gain of j-th user in i-th cell dis[i-1][j-1]
    P_all = []
    B_all = [] 
    R = np.zeros(K)
    path_loss = 2
    N_k = np.random.poisson(means)
    for i in range(K):
        tmpdis = np.random.uniform(0, radii[i], N_k[i]);
        tmpcha = (2.0) ** 0.5 / 2.0 * (np.random.rand(N_k[i]) + 1j * np.random.rand(N_k[i]))
        tmpgain = tmpdis ** path_loss * np.linalg.norm(tmpcha) ** 2
        dis.append(tmpdis)
        cha.append(tmpcha)
        gain.append(tmpgain)
        w = np.ones(N_k[i]);
        status, R[i], P, B = opt_allocatin(w, tmpgain, P_max[i], B_max[i])
        P_all.append(P)
        B_all.append(B)
    

def main2():
    N = 10;
    R = 5;
    path_loss = 2
    r = np.random.uniform(0, R, N);
    cha = (2.0) ** 0.5 / 2.0 * (np.random.rand(N) + 1j * np.random.rand(N))
    gain = r ** path_loss * np.linalg.norm(cha) ** 2
    SNR_all = np.arange(10)
    P_max_all = sp.db2pow(SNR_all);
    # P_max_all = np.arange(10)*100
    B_max = 1;
    R = np.zeros(len(P_max_all))
    R2 = np.zeros(len(P_max_all))
    for i in range(len(R)):
        status, R[i], P, B = opt_allocatin(np.ones(N), gain, P_max_all[i], B_max)
        R2[i] = uniform_allocation(gain, P_max_all[i], B_max)
    R = R / np.log(2)
    plot.rc('text', usetex=True)
    plot.figure
    opt_curve, = plot.plot(P_max_all, R, '-o')
    uni_curve, = plot.plot(P_max_all, R2, '-o')
    plot.xlabel('P')
    plot.ylabel('R')
    plot.legend([opt_curve, uni_curve], ['Opt', 'Uni'])
    plot.show()
    
        
if __name__ == '__main__':
    main2()