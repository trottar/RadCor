# import matplotlib
# matplotlib.use('Agg') 
import matplotlib.pyplot as plt
import math as m
import numpy as np

energyStr = "10380"
Ebeam = 10.38
theta_values = [11, 18, 30]
theta = theta_values[0]
thetaStr = str(theta)
deg2rad = 0.0174533
Mp = 0.938
degree_sign = u'\N{DEGREE SIGN}'
unpolLabel = "_unpol"
unpolLabel2 = " Unpolarized"
longLabel = "_long"
longLabel2 = " Longitudinal"
transLabel = "_trans"
transLabel2 = " Transverse"
AparLabel = "_Apar"
AparLabel2 = r'$\mathrm{A}_{\parallel}^{^{3}He}$'
AperpLabel = "_Aperp"
AperpLabel2 = r'$\mathrm{A}_{\perp}^{^{3}He}$'

def q2_calc(Ep, theta):
    return 2*Ebeam*Ep*(1 - m.cos(theta*deg2rad))

def x_calc(Ep, theta):
    Q2 = q2_calc(Ep, theta)
    nu = Ebeam - Ep
    denom = 2*Mp*nu
    return Q2/denom

def main():

    # Theta value loop
    for theta_value in theta_values:

        theta = theta_value
        thetaStr = str(theta)

        print("\nCreating Plots for theta = " + thetaStr)
        for i in range(30): print("-", end='')
        print()

        # Get the input filename from the user
        inUnpolFileName = "CAnalyzer-master/example/output/radiated_model_" + thetaStr + "deg" + unpolLabel + ".dat"
        inLongFileName = "CAnalyzer-master/example/output/radiated_model_" + thetaStr + "deg" + longLabel + ".dat"
        inTransFileName = "CAnalyzer-master/example/output/radiated_model_" + thetaStr + "deg" + transLabel + ".dat"
        outUnpolFileName = "Plots/" + thetaStr + "deg" + unpolLabel + ".png"
        outLongFileName = "Plots/" + thetaStr + "deg" + longLabel + ".png"
        outTransFileName = "Plots/" + thetaStr + "deg" + transLabel + ".png"
        outAparFileName = "Plots/" + thetaStr + "deg" + AparLabel + ".png"
        outAperpFileName = "Plots/" + thetaStr + "deg" + AperpLabel + ".png"
        outApardiffFileName = "Plots/" + thetaStr + "deg" + AparLabel + "diff.png"
        outAperpdiffFileName = "Plots/" + thetaStr + "deg" + AperpLabel + "diff.png"

        # Read data from the input file
        try:
            with open(inUnpolFileName, 'r') as infile:
                next(infile)
                dataUnpol = []
                for line in infile:
                    line.replace("\n","")
                    vals = line.strip().split(" ")
                    while "" in vals:
                        vals.remove("")
                    for i in range(len(vals)):
                        vals[i] = float(vals[i])    # Ebeam, nu, xs_rad, xs_born, stat, syst
                    if (vals[0] != 10380 or vals[1] > (Ebeam-1)*1000): continue
                    dataUnpol.append(vals)
        except IOError:
            print("Unable to open file:", inUnpolFileName)
            return

        try:
            with open(inLongFileName, 'r') as infile:
                next(infile)
                dataLong = []
                for line in infile:
                    line.replace("\n","")
                    vals = line.strip().split(" ")
                    while "" in vals:
                        vals.remove("")
                    for i in range(len(vals)):
                        vals[i] = float(vals[i])    # Ebeam, nu, xs_rad, xs_born, stat, syst
                    if (vals[0] != 10380 or vals[1] > (Ebeam-1)*1000): continue
                    dataLong.append(vals)
        except IOError:
            print("Unable to open file:", inLongFileName)
            return

        try:
            with open(inTransFileName, 'r') as infile:
                next(infile)
                dataTrans = []
                for line in infile:
                    line.replace("\n","")
                    vals = line.strip().split(" ")
                    while "" in vals:
                        vals.remove("")
                    for i in range(len(vals)):
                        vals[i] = float(vals[i])    # Ebeam, nu, xs_rad, xs_born, stat, syst
                    if (vals[0] != 10380 or vals[1] > (Ebeam-1)*1000): continue
                    dataTrans.append(vals)
        except IOError:
            print("Unable to open file:", inTransFileName)
            return

        if (theta == 11):
            xvals = np.array([0.394,0.504,0.663,0.912])
            aparvals = np.array([-0.0044,-0.0000,0.0042,-0.001])
            aparerrs = np.array([0.0016,0.0019,0.0028,0.006])
            aperpvals = np.array([0.0007,-0.0001,-0.0023,-0.0073])
            aperperrs = np.array([0.0006,0.0007,0.0011,0.0023])
        elif (theta == 18):
            xvals = np.array([0.542,0.637,0.747,0.885])
            aparvals = np.array([-0.002,-0.002,0.003,-0.008])
            aparerrs = np.array([0.005,0.006,0.010,0.020])
            aperpvals = np.array([-0.0004,-0.0014,0.0050,-0.018])
            aperperrs = np.array([0.0016,0.0020,0.0032,0.006])
        else:
            xvals = np.array([0.423,0.440,0.465,0.492,0.517,0.544,0.571,0.599,0.629,0.647,0.676,0.708,0.737,0.767,0.799,0.834,0.867,0.901])
            aparvals = np.array([-0.00615,-0.00233,-0.00462,0.00448,0.00402,0.00348,0.00753,0.00873,0.00217,0.00302,0.00387,0.00429,0.00777,0.008,0.00812,0.00501,0.00717,0.01247])
            aparerrstat = np.array([0.00711,0.003,0.00326,0.00356,0.00242,0.00262,0.00285,0.00323,0.00189,0.002,0.0018,0.00206,0.00239,0.00282,0.00346,0.00584,0.00907,0.01179])
            aparerrsyst = np.array([0.00078,0.00125,0.00138,0.0015,0.00106,0.00097,0.00095,0.00114,0.0008,0.00078,0.00072,0.00031,0.0005,0.00037,0.00055,0.00072,0.00111,0.00141])
            aparerrs = np.sqrt(aparerrstat**2 + aparerrsyst**2)
            aperpvals = np.array([-0.02459,0.00831,0.00093,-0.01353,0.00866,0.01762,0.01071,0.0152,0.00809,0.005,0.00936,0.00627,0.00723,0.01034,0.01797,-0.0043,-0.00852,0.00376])
            aperperrstat = np.array([0.01754,0.00729,0.00787,0.00864,0.00626,0.00666,0.00746,0.00823,0.0048,0.00463,0.0039,0.00448,0.00518,0.00607,0.00759,0.01359,0.02221,0.02985])
            aperperrsyst = np.array([0.00128,0.00066,0.00085,0.00154,0.00279,0.00281,0.00332,0.00517,0.00299,0.00149,0.00134,0.00174,0.00254,0.00309,0.00378,0.00615,0.00948,0.01257])
            aperperrs = np.sqrt(aperperrstat**2 + aperperrsyst**2)

        xsradUnpol = [[],[]]
        xsbornUnpol = [[],[]]
        xsradLong = [[],[]]
        xsbornLong = [[],[]]
        xsradTrans = [[],[]]
        xsbornTrans = [[],[]]

        Apar_rad = [[],[]]
        Apar_born = [[],[]]
        Aperp_rad = [[],[]]
        Aperp_born = [[],[]]
        Apar_diff = [[],[]]
        Aperp_diff = [[],[]]

        for d in dataUnpol:
            xsradUnpol[0].append(Ebeam - d[1]/1000)
            xsradUnpol[1].append(d[2])
            xsbornUnpol[0].append(Ebeam - d[1]/1000)
            xsbornUnpol[1].append(d[3])
        for d in dataLong:
            xsradLong[0].append(Ebeam - d[1]/1000)
            xsradLong[1].append(d[2])
            xsbornLong[0].append(Ebeam - d[1]/1000)
            xsbornLong[1].append(d[3])
        for d in dataTrans:
            xsradTrans[0].append(Ebeam - d[1]/1000)
            xsradTrans[1].append(d[2])
            xsbornTrans[0].append(Ebeam - d[1]/1000)
            xsbornTrans[1].append(d[3])

        for i in range(len(dataLong)):
            if (dataLong[i][3] == 0): continue
            tempEp = Ebeam - (dataLong[i][1] + dataUnpol[i][1])/2000
            tempx = x_calc(tempEp, theta)
            Apar_rad[0].append(tempx)
            Apar_born[0].append(tempx)
            Apar_rad[1].append(dataLong[i][2]/(2.0*dataUnpol[i][2]))
            Apar_born[1].append(dataLong[i][3]/(2.0*dataUnpol[i][3]))
            Apar_diff[0].append(tempx)
            Apar_diff[1].append(dataLong[i][3]/(2.0*dataUnpol[i][3]) - dataLong[i][2]/(2.0*dataUnpol[i][2]))
            
        for i in range(len(dataTrans)):
            if (dataTrans[i][3] == 0): continue
            tempEp = Ebeam - (dataTrans[i][1] + dataUnpol[i][1])/2000
            tempx = x_calc(tempEp, theta)
            Aperp_rad[0].append(tempx)
            Aperp_born[0].append(tempx)
            Aperp_rad[1].append(dataTrans[i][2]/(2.0*dataUnpol[i][2]))
            Aperp_born[1].append(dataTrans[i][3]/(2.0*dataUnpol[i][3]))
            Aperp_diff[0].append(tempx)
            Aperp_diff[1].append(dataTrans[i][3]/(2.0*dataUnpol[i][3]) - dataTrans[i][2]/(2.0*dataUnpol[i][2]))
                
        horizLineUnpol = [[xsradUnpol[0][0],xsradUnpol[0][len(xsradUnpol[0])-1]],[0,0]]
        horizLineLong = [[xsradLong[0][0],xsradLong[0][len(xsradLong[0])-1]],[0,0]]
        horizLineTrans = [[xsradTrans[0][0],xsradTrans[0][len(xsradTrans[0])-1]],[0,0]]
        horizLineApar = [[Apar_born[0][0],Apar_born[0][len(Apar_rad[0])-1]],[0,0]]
        horizLineAperp = [[Aperp_rad[0][0],Aperp_rad[0][len(Aperp_rad[0])-1]],[0,0]]

        plt.figure(figsize=(10, 5))  # Adjust width and height as needed

        # Plot the data
        plt.plot(xsradUnpol[0], xsradUnpol[1], color="blue", label='XS rad')
        # plt.plot(xsbornUnpol[0], xsbornUnpol[1], color="blue", label='XS born', linestyle='--')
        plt.plot(xsbornUnpol[0], xsbornUnpol[1], color="green", label='XS born')        
        # plt.plot(xsraw[0], xsraw[1], color='black', label='XS_raw')
        plt.plot(horizLineUnpol[0], horizLineUnpol[1], color='black', linestyle='--')
        plt.xlabel('E\' [GeV]')
        plt.ylabel('XS')
        plt.title(unpolLabel2 + ' E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign)
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outUnpolFileName, format='png')
        print("Plot saved to " + outUnpolFileName)

        plt.clf()

        # Plot the data
        plt.plot(xsradLong[0], xsradLong[1], color="blue", label='XS rad')
        # plt.plot(xsbornLong[0], xsbornLong[1], color="blue", label='XS born', linestyle='--')
        plt.plot(xsbornLong[0], xsbornLong[1], color="green", label='XS born')
        # plt.plot(xsraw[0], xsraw[1], color='black', label='XS_raw')
        plt.plot(horizLineLong[0], horizLineLong[1], color='black', linestyle='--')
        plt.xlabel('E\' [GeV]')
        plt.ylabel('XS')
        plt.title(longLabel2 + ' E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign)
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outLongFileName, format='png')
        print("Plot saved to " + outLongFileName)

        plt.clf()

        # Plot the data
        plt.plot(xsradTrans[0], xsradTrans[1], color="blue", label='XS rad')
        # plt.plot(xsbornTrans[0], xsbornTrans[1], color="blue", label='XS born', linestyle='--')
        plt.plot(xsbornTrans[0], xsbornTrans[1], color="green", label='XS born')
        # plt.plot(xsraw[0], xsraw[1], color='black', label='XS_raw')
        plt.plot(horizLineTrans[0], horizLineTrans[1], color='black', linestyle='--')
        plt.xlabel('E\' [GeV]')
        plt.ylabel('XS')
        plt.title(transLabel2 + ' E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign)
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outTransFileName, format='png')
        print("Plot saved to " + outTransFileName)

        plt.clf()

        # Plot the data
        plt.plot(Apar_rad[0], Apar_rad[1], color="blue", label=AparLabel2+' rad')
        # plt.plot(Apar_born[0], Apar_born[1], color="blue", label=AparLabel2+' born', linestyle='--')
        plt.plot(Apar_born[0], Apar_born[1], color="green", label=AparLabel2+' born')
        # plt.plot(xsraw[0], xsraw[1], color='black', label='XS_raw')
        plt.plot(horizLineApar[0], horizLineApar[1], color='black', linestyle='--')
        plt.errorbar(xvals, aparvals, yerr=aparerrs, fmt='o', color='red', markersize=5, capsize=5, label=AparLabel2+' data')
        plt.xlabel('x')
        plt.ylabel(AparLabel2)
        plt.title(AparLabel2 + ' E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign)
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outAparFileName, format='png')
        print("Plot saved to " + outAparFileName)

        plt.clf()

        # Plot the data
        plt.plot(Aperp_rad[0], Aperp_rad[1], color="blue", label=AperpLabel2+' rad')
        # plt.plot(Aperp_born[0], Aperp_born[1], color="blue", label=AperpLabel2+' born', linestyle='--')
        plt.plot(Aperp_born[0], Aperp_born[1], color="green", label=AperpLabel2+' born')
        # plt.plot(xsraw[0], xsraw[1], color='black', label='XS_raw')
        plt.plot(horizLineAperp[0], horizLineAperp[1], color='black', linestyle='--')
        plt.errorbar(xvals, aperpvals, yerr=aperperrs, fmt='o', color='red', markersize=5, capsize=5, label=AperpLabel2+' data')
        plt.xlabel('x')
        plt.ylabel(AperpLabel2)
        plt.title(AperpLabel2 + ' E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign)
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outAperpFileName, format='png')
        print("Plot saved to " + outAperpFileName)

        plt.clf()

        # Plot the data
        plt.plot(Apar_diff[0], Apar_diff[1], color="blue", label=AparLabel2+'(born - rad)')
        # plt.plot(Aperp_born[0], Aperp_born[1], color="blue", label=AperpLabel2+' born', linestyle='--')
        # plt.plot(Aperp_born[0], Aperp_born[1], color="green", label=AperpLabel2+' born')
        # plt.plot(xsraw[0], xsraw[1], color='black', label='XS_raw')
        plt.plot(horizLineAperp[0], horizLineAperp[1], color='black', linestyle='--')
        # plt.errorbar(xvals, aperpvals, yerr=aperperrs, fmt='o', color='red', markersize=5, capsize=5, label=AperpLabel2+' data')
        plt.xlabel('x')
        plt.ylabel(AparLabel2 + '(born - rand)')
        plt.title(AparLabel2 + '(born - rad) E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign)
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outApardiffFileName, format='png')
        print("Plot saved to " + outApardiffFileName)

        plt.clf()

        # Plot the data
        plt.plot(Aperp_diff[0], Aperp_diff[1], color="blue", label=AperpLabel2+'(born - rad)')
        # plt.plot(Aperp_born[0], Aperp_born[1], color="blue", label=AperpLabel2+' born', linestyle='--')
        # plt.plot(Aperp_born[0], Aperp_born[1], color="green", label=AperpLabel2+' born')
        # plt.plot(xsraw[0], xsraw[1], color='black', label='XS_raw')
        plt.plot(horizLineAperp[0], horizLineAperp[1], color='black', linestyle='--')
        # plt.errorbar(xvals, aperpvals, yerr=aperperrs, fmt='o', color='red', markersize=5, capsize=5, label=AperpLabel2+' data')
        plt.xlabel('x')
        plt.ylabel(AperpLabel2 + '(born - rand)')
        plt.title(AperpLabel2 + '(born - rad) E = ' + energyStr + ' MeV 'r'$\theta$ = ' + thetaStr + degree_sign)
        plt.grid(False)
        plt.legend()
        # plt.show()
        plt.tight_layout()
        plt.savefig(outAperpdiffFileName, format='png')
        print("Plot saved to " + outAperpdiffFileName)

        plt.clf()
    
    for i in range(30): print("-", end='')
    print()




if __name__ == "__main__":
    main()