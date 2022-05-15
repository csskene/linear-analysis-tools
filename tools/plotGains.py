# Plot the obtained singular values with cubic splines
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import matplotlib.colors as mcolors


def cubeFit(x1,x2,fx1,fx2,fx1d,fx2d):
    # Finds coefficients f(x) = Ax^3+Bx^2+Cx+D such that
    # f(x1)  = fx1
    # f(x2)  = fx2
    # f'(x1) = fx1d
    # f'(x2) = fx2d

    RHS = np.array([fx1,fx2,fx1d,fx2d]);
    row1 = np.array([x1**3, x1**2, x1, 1])
    row2 = np.array([x2**3, x2**2, x2, 1])
    row3 = np.array([3*x1**2, 2*x1, 1, 0])
    row4 = np.array([3*x2**2, 2*x2, 1, 0])
    LHS = np.vstack((row1,row2,row3,row4))
    ANS = np.linalg.solve(LHS,RHS)
    A = ANS[0]
    B = ANS[1]
    C = ANS[2]
    D = ANS[3]
    return A,B,C,D


if __name__ == '__main__':

    # Sort out the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-gainsFile',type=str,help='Folder where the outputted gains file is.')
    parser.add_argument('-mode',type=int,help='Plot only this mode.')
    args = parser.parse_args()
    gainsFile = args.gainsFile
    mode      = args.mode

    if(not args.mode):
        mode = 0
        print('Plotting only the leading singular value (default)')
    elif(args.mode and mode ==-1):
        print('Plotting the gains for all modes.')
    else:
        print('Plotting the gains for mode %d' % mode)

    # Plot settings
    plt.rc('text', usetex=True)
    plt.rc('font', family='times',size=14)

    # Read in the data as a pandas data frame
    df = pd.read_csv(gainsFile,sep="\s+",header=0)

    colourList = list(mcolors.TABLEAU_COLORS.values())
    modes = []
    if(mode != -1):
        modes.append(mode)
    else:
        modes = list(set(df['i']))

    # For the legend
    lines = []
    labels = []
    for modeNum in modes:
        # Get only the mode of interest
        interMode = df['i'] == modeNum
        dfmode1 = df[interMode].reset_index()
        dfmode1 = dfmode1.sort_values('Omega').reset_index()

        # Fit cubic splines between consecutive singular values and plot
        for iter in range(len(dfmode1['Omega'])-1):
             omL    = dfmode1.loc[iter,'Omega']
             omR    = dfmode1.loc[iter+1,'Omega']
             gainL  = dfmode1.loc[iter,'sigma']
             gainR  = dfmode1.loc[iter+1,'sigma']
             dgainL = dfmode1.loc[iter,'dsigma/domega']
             dgainR = dfmode1.loc[iter+1,'dsigma/domega']
             A,B,C,D = cubeFit(omL,omR,gainL,gainR,dgainL,dgainR)
             x = np.linspace(omL,omR,100)
             f = A*x**3+B*x**2+C*x+D
             if(iter==0):
                 line, = plt.plot(x,f,colourList[modeNum])
                 lines.append(line)
                 labels.append('Mode = '+str(modeNum))
             else:
                 plt.plot(x,f,colourList[modeNum])

    plt.xlabel(r'$\omega$')
    plt.ylabel(r'$\sigma$')
    plt.legend(lines, labels)
    if(mode!=1):
        fileName = 'gainsMode'+str(mode)
    else:
        fileName = 'gainsAllModes'
    plt.savefig(fileName,dpi=300,bbox_inches='tight')
