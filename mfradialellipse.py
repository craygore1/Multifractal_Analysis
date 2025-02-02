import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import solve
from scipy.stats import linregress

def mfradialellipse(bwimage, objectcolor, qvals, plots):
    sz = bwimage.shape
    
    if objectcolor:
        mylog = bwimage == 1
    else:
        mylog = bwimage == 0
    
    centerx, centery = sz[1] // 2, sz[0] // 2
    radiusx, radiusy = centerx, centery
    fullradius = radiusx + radiusy
    
    imcols, imrows = np.meshgrid(np.arange(1, sz[1] + 1), np.arange(1, sz[0] + 1))
    myellipse = ((imrows - centery)**2 / radiusy**2 + (imcols - centerx)**2 / radiusx**2) <= 1
    mylog = mylog & myellipse
    
    xvec = np.round(np.linspace(-sz[1] // 2, sz[1] // 2, sz[1]))
    yvec = np.round(np.linspace(-sz[0] // 2, sz[0] // 2, sz[0]))
    Xim, Yim = np.meshgrid(xvec, yvec)
    
    theta, rho = np.arctan2(Yim, Xim), np.hypot(Xim, Yim)
    
    myind = np.where(mylog)
    thetaind, rhoind = theta[myind], rho[myind]
    Xind, Yind = Xim[myind], Yim[myind]
    
    totalarea = np.pi * max(xvec) * max(yvec)
    aspect = max(max(xvec), max(yvec)) / min(max(xvec), max(yvec))
    
    maxboxes = 5
    M = np.zeros((4**maxboxes, maxboxes + 1))
    
    for i in range(maxboxes + 1):
        d1 = np.pi * np.ones(2**i)
        d2 = -np.pi * np.ones(2**i - 1)
        A = np.diag(d1) + np.diag(d2, -1)
        
        areas = (totalarea / (2**i)) * np.ones(2**i)
        ab = solve(A, areas)
        minorrange = np.sqrt(ab / aspect)
        majorrange = ab / minorrange
        
        minorrange = np.insert(minorrange, 0, 0)
        majorrange = np.insert(majorrange, 0, 0)
        
        thetarange = np.linspace(-np.pi, np.pi, (2**i) + 1)
        
        counter = 0
        for j in range(len(thetarange) - 1):
            for k in range(len(majorrange) - 1):
                temp1 = (thetarange[j] <= thetaind) & (thetaind < thetarange[j+1])
                temp2 = (1 <= (Xind**2 / majorrange[k]**2 + Yind**2 / minorrange[k]**2)) & \
                        ((Xind**2 / majorrange[k+1]**2 + Yind**2 / minorrange[k+1]**2) < 1)
                temp3 = temp1 & temp2
                M[counter, i] = np.sum(temp3)
                counter += 1
    
    prbM = M / (np.sum(M, axis=0) + 1e-10)
    truesz = (2 * np.pi) / (2 ** np.arange(maxboxes + 1))
    X = np.log2(truesz)
    
    yD, yalph, yf = np.zeros((len(qvals), len(prbM[0]))), np.zeros((len(qvals), len(prbM[0]))), np.zeros((len(qvals), len(prbM[0])))
    
    for idx, k in enumerate(qvals):
        for a in range(len(prbM[0])):
            nonzero_prbM = prbM[:, a][prbM[:, a] > 0]
            if k == 1:
                yD[idx, a] = np.sum(nonzero_prbM * np.log2(nonzero_prbM))
            else:
                yD[idx, a] = np.log2(np.maximum(np.sum(nonzero_prbM**k), 1e-10))
                mu = (nonzero_prbM**k) / np.sum(nonzero_prbM**k)
                yalph[idx, a] = np.sum(mu * np.log2(nonzero_prbM))
                yf[idx, a] = np.sum(mu * np.log2(mu))
    
    Dq, tauq, myalpha, falpha = np.zeros(len(qvals)), np.zeros(len(qvals)), np.zeros(len(qvals)), np.zeros(len(qvals))
    
    for currq in range(len(qvals)):
        slope, _, _, _, _ = linregress(X, yD[currq])
        if qvals[currq] == 1:
            Dq[currq] = abs(slope)
        else:
            tauq[currq] = slope
            Dq[currq] = tauq[currq] / (qvals[currq] - 1)
            myalpha[currq] = linregress(X, yalph[currq])[0]
            falpha[currq] = linregress(X, yf[currq])[0]
    
    h = qvals[-1] - qvals[-2]
    alphaleg = np.gradient(tauq, h)
    fleg = qvals * alphaleg - tauq
    
    myalpha, falpha = alphaleg, fleg
    
    if plots:
        plt.figure()
        plt.plot(qvals, Dq, color='#0F6FC6', linewidth=1.25)
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xlabel('$q$', fontsize=16)
        plt.ylabel('$D(q)$', fontsize=16)
        plt.show()
        
        plt.figure()
        plt.scatter(myalpha, falpha, color='#0F6FC6', marker='.')
        plt.grid(True, which='both', linestyle='--', linewidth=0.5)
        plt.xlabel('$\alpha$', fontsize=16)
        plt.ylabel('$f(\alpha)$', fontsize=16)
        plt.show()
    
    return Dq, myalpha, falpha