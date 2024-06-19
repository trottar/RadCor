#! /usr/bin/python

#
# Description:
# ================================================================
# Time-stamp: "2024-06-19 11:42:30 trottar"
# ================================================================
#
# Author:  Richard L. Trotta III <trottar.iii@gmail.com>
#
# Copyright (c) trottar
#

import csv
import numpy as np
from collections import defaultdict
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from matplotlib.backends.backend_pdf import PdfPages
import sys

class Point:
    def __init__(self, x=0, q2=0, F1=0, F2=0, g1=0, g2=0):
        self.x = x
        self.q2 = q2
        self.F1 = F1
        self.F2 = F2
        self.g1 = g1
        self.g2 = g2

    def set_vals(self, x, q2, F1, F2, g1, g2):
        self.x = x
        self.q2 = q2
        self.F1 = F1
        self.F2 = F2
        self.g1 = g1
        self.g2 = g2

    def set_equal(self, pt):
        self.set_vals(pt.x, pt.q2, pt.F1, pt.F2, pt.g1, pt.g2)

def get_data_points(file_name):
    data_points = defaultdict(list)
    q2_keys = []
    
    with open(file_name, 'r') as jam_file:
        reader = csv.DictReader(jam_file)
        for row in reader:
            pt = Point()

            #col_type = "Inelastic"            
            #col_type = "QE"
            col_type = "IpQE"
            
            if file_name == "table_3He_JAM_smeared_no_QE_ipol1IA14_SF23_AC11.csv":
                pt.q2 = float(row['Q2'])
                pt.x = float(row['X'])
                pt.F2 = float(row['F2'])
                pt.g1 = float(row['G1'])
                pt.g2 = float(row['G2'])
                pt.F1 = float(row['F1'])
            elif file_name == "XZ_table_3He_JAM_smeared_kpsv_onshell_ipol1_ipolres1_IA14_SF23_AC11.csv":
                pt.q2 = float(row['Q2'])
                pt.x = float(row['X'])
                pt.F2 = float(row[f'F2_{col_type}'])
                pt.F1 = float(row[f'F1_{col_type}'])
                pt.g1 = float(row[f'G1_{col_type}'])
                pt.g2 = float(row[f'G2_{col_type}'])
            else:
                print(f"ERROR: Not a valid file name, {file_name}")
                sys.exit(2)
                
            if pt.q2 not in data_points:
                data_points[pt.q2] = []
            data_points[pt.q2].append(pt)
            
            if pt.q2 not in q2_keys:
                q2_keys.append(pt.q2)
    
    q2_keys.sort()
    return data_points, q2_keys

def interpolate_2d_linear(x, q2, data_points, q2_keys):
    pt = Point()

    q2min, q2max = 0, 0
    for q2_key in q2_keys:
        if q2_key <= q2:
            q2min = q2_key
        if q2_key > q2:
            q2max = q2_key
            break

    def find_x_range(data, x):
        ixmin, ixmax = 0, 0
        for i in range(len(data)):
            if data[i].x <= x:
                ixmin = i
            if data[i].x > x:
                ixmax = i
                break
        return ixmin, ixmax

    ixmin1, ixmax1 = find_x_range(data_points[q2min], x)
    ixmin2, ixmax2 = find_x_range(data_points[q2max], x)

    xmin1 = data_points[q2min][ixmin1].x
    xmax1 = data_points[q2min][ixmax1].x
    xmin2 = data_points[q2max][ixmin2].x
    xmax2 = data_points[q2max][ixmax2].x

    def interpolate(val_min, val_max, xmin, xmax):
        return (xmax - x) / (xmax - xmin) * val_min + (x - xmin) / (xmax - xmin) * val_max

    F1 = interpolate(data_points[q2min][ixmin1].F1, data_points[q2min][ixmax1].F1, xmin2, xmax2) * (q2max - q2) / (q2max - q2min) + \
         interpolate(data_points[q2max][ixmin1].F1, data_points[q2max][ixmax1].F1, xmin1, xmax1) * (q2 - q2min) / (q2max - q2min)

    F2 = interpolate(data_points[q2min][ixmin1].F2, data_points[q2min][ixmax1].F2, xmin2, xmax2) * (q2max - q2) / (q2max - q2min) + \
         interpolate(data_points[q2max][ixmin1].F2, data_points[q2max][ixmax1].F2, xmin1, xmax1) * (q2 - q2min) / (q2max - q2min)

    g1 = interpolate(data_points[q2min][ixmin1].g1, data_points[q2min][ixmax1].g1, xmin2, xmax2) * (q2max - q2) / (q2max - q2min) + \
         interpolate(data_points[q2max][ixmin1].g1, data_points[q2max][ixmax1].g1, xmin1, xmax1) * (q2 - q2min) / (q2max - q2min)

    g2 = interpolate(data_points[q2min][ixmin1].g2, data_points[q2min][ixmax1].g2, xmin2, xmax2) * (q2max - q2) / (q2max - q2min) + \
         interpolate(data_points[q2max][ixmin1].g2, data_points[q2max][ixmax1].g2, xmin1, xmax1) * (q2 - q2min) / (q2max - q2min)

    pt.set_vals(x, q2, F1, F2, g1, g2)
    return pt

# Proton and neutron mass (GeV)
mp = 0.938 
mn = 0.939

# Set max value of x
xmax = 1.5

# Usage
#data_points, q2_keys = get_data_points("table_3He_JAM_smeared_no_QE_ipol1IA14_SF23_AC11.csv")
data_points, q2_keys = get_data_points("XZ_table_3He_JAM_smeared_kpsv_onshell_ipol1_ipolres1_IA14_SF23_AC11.csv")

# Interpolated values for different x and q2
q2_values_gev = [0.1, 0.3, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0]
x_values = np.linspace(0.000, xmax, 1500)

# Store results
interpolated_points = []

for x in x_values:
    for q2 in q2_values:
        pt = interpolate_2d_linear(x, q2, data_points, q2_keys)
        interpolated_points.append(pt)
        
# Define the deviation threshold
#threshold = 55.0 # table_3He_JAM_smeared_no_QE_ipol1IA14_SF23_AC11.csv
threshold = 100.0

deviation_points = []

# Initialize a PDF file
with PdfPages('strufunc_jam.pdf') as pdf:

    # Plot the results
    for q2 in q2_values:
        F1_vals = [pt.F1 for pt in interpolated_points if pt.q2 == q2]
        F2_vals = [pt.F2 for pt in interpolated_points if pt.q2 == q2]
        g1_vals = [pt.g1 for pt in interpolated_points if pt.q2 == q2]
        g2_vals = [pt.g2 for pt in interpolated_points if pt.q2 == q2]
        x_vals = [pt.x for pt in interpolated_points if pt.q2 == q2]

        plt.figure(figsize=(12, 8))
        #plt.plot(x_vals, F1_vals, label=f'F1')
        plt.plot(x_vals, F2_vals, label=f'F2')
        #plt.plot(x_vals, g1_vals, label=f'g1')
        #plt.plot(x_vals, g2_vals, label=f'g2')

        # Find where F1 deviates greatly from 0.0
        deviation_x = None
        for i in range(1, len(F1_vals)):  # Start from the second element to avoid index error
            percent_change = abs((F1_vals[i] - F1_vals[i-1]) / F1_vals[i-1]) * 100
            if percent_change > threshold and x_vals[i] > 0.01:
                deviation_x = x_vals[i-1]
                #deviation_x = x_vals[i-10]  # 10 bins before spike
                plt.axvline(x=deviation_x, color='r', linestyle='--', label=f'x={deviation_x:.3f}')
                break

        if deviation_x is not None:
            deviation_points.append((q2, deviation_x))
        else:
            deviation_points.append((q2, xmax))

        plt.xlabel('x')
        plt.ylabel('Values')
        plt.title(f'Interpolated $Q^2$={q2}')
        plt.legend()
        plt.grid(True)
        pdf.savefig()  # Save the figure to the PDF
        plt.close()

    # List of colors
    colors = [
        "blue", "green", "red", "cyan", "magenta", 
        "black", "gray", "purple", "brown", "orange", 
        "pink", "olive", "lime", "navy", "teal", 
        "maroon", "silver"
    ]
    
    # Assuming deviation_points is a list of tuples (q2, x)
    deviation_points = [(pt[0], pt[1]) for pt in deviation_points]

    try:
        # Extract q2 and x data from deviation_points
        q2_deviation = np.array([pt[0] for pt in deviation_points])
        x_deviation = np.array([pt[1] for pt in deviation_points])

        # Define your exponential function with fixed c parameter
        def exp_func_mp(x, a, b):
            c = mp**2  # Fixed to proton mass squared
            return a * np.exp(b * x) + c

        # Define your exponential function with fixed c parameter
        def exp_func_mn(x, a, b):
            c = mn**2  # Fixed to neutron mass squared
            return a * np.exp(b * x) + c

        # Define your exponential function
        def exp_func(x, a, b):
            return a * np.exp(b * x)

        # Define your exponential function
        def poly_func(x, a, b, c):
            #return a + b * x + c * x**2
            #return (1-x)*x+a*np.exp(b*x)
            #return (1-x)*a*np.exp(b/(1-x))
            return x*((1-a*x)**b)*np.exp(c*x)

        # Fit the exponential function to the deviation points
        popt_mp, pcov_mp = curve_fit(exp_func_mp, x_deviation, q2_deviation, maxfev=2000)
        popt_mn, pcov_mn = curve_fit(exp_func_mn, x_deviation, q2_deviation, maxfev=2000)
        popt, pcov = curve_fit(exp_func, x_deviation, q2_deviation, maxfev=2000)
        popt_poly, pcov_poly = curve_fit(poly_func, x_deviation, q2_deviation, maxfev=2000)

        # Generate points for the fitted function
        x_fit = np.linspace(min(x_deviation), max(x_deviation), 100)
        q2_fit_mp = exp_func_mp(x_fit, *popt_mp)
        q2_fit_mn = exp_func_mn(x_fit, *popt_mn)
        q2_fit = exp_func(x_fit, *popt)
        q2_fit_poly = poly_func(x_fit, *popt_poly)

        # Plot deviation points and the fitted function
        plt.figure(figsize=(12, 8))
        plt.plot(x_fit, q2_fit_mp, linestyle='-.', color='purple', label=f'Fit: {popt_mp[0]:.3e} * exp({popt_mp[1]:.3e} * x) + $m_p^2$')
        plt.plot(x_fit, q2_fit_mn, linestyle='--', color='red', label=f'Fit: {popt_mn[0]:.3e} * exp({popt_mn[1]:.3e} * x) + $m_n^2$')
        plt.plot(x_fit, q2_fit, linestyle='-', color='green', label=f'Fit: {popt[0]:.3e} * exp({popt[1]:.3e} * x)')
        #plt.plot(x_fit, q2_fit_poly, linestyle='-', color='orange', label=f'Fit: x * (1 - {popt_poly[0]:.3e} * x)^({popt_poly[1]:.3e}) * exp({popt_poly[2]:.3e} * x)')
        plt.plot(x_fit, q2_fit_poly, linestyle='-', color='orange', label=f'Fit: x * (1 - {popt_poly[0]:.3f} * x)^({popt_poly[1]:.3f}) * exp({popt_poly[2]:.3f} * x)')
        plt.scatter(x_deviation, q2_deviation, marker='o', linestyle='-', color='b', label='Deviation Points')

        plt.axvline(x=1.000, color='r', linestyle='--', label=f'x={1.000:.3f}')

        plt.xlabel('x')
        plt.ylabel('$Q^2$')
        plt.xscale('log')
        plt.yscale('log')
        plt.title('Exponential Fit of Deviation Points')
        plt.legend()
        plt.grid(True)
        pdf.savefig()  # Save the figure to the PDF
        plt.close()
    except TypeError:
        print("Fit failed...")
    except RuntimeError:
        print("Fit failed...")            

    i=1
    plt.figure(figsize=(12, 8))
    # Plot the results
    for q2 in q2_values:
        filtered_x = [pt[1] for pt in deviation_points if pt[0] == q2]
        F1_vals = [pt.F1 for pt in interpolated_points if pt.q2 == q2 and pt.x < filtered_x]
        F2_vals = [pt.F2 for pt in interpolated_points if pt.q2 == q2 and pt.x < filtered_x]
        g1_vals = [pt.g1 for pt in interpolated_points if pt.q2 == q2 and pt.x < filtered_x]
        g2_vals = [pt.g2 for pt in interpolated_points if pt.q2 == q2 and pt.x < filtered_x]
        x_vals = [pt.x for pt in interpolated_points if pt.q2 == q2 and pt.x < filtered_x]

        plt.plot(x_vals, F2_vals, label=f'$Q^2$={q2}', color=colors[i])    
        i+=1

    plt.xlabel('x')
    plt.ylabel('F2')
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    pdf.savefig()  # Save the figure to the PDF
    plt.close()

        
    x_lst = np.linspace(0.01, xmax, 15)

    i=1
    plt.figure(figsize=(12, 8))
    for x in x_lst:
        # Filter points for the current x
        filtered_points = [pt for pt in interpolated_points if x - 0.01 < pt.x < x + 0.01]

        # Group by Q^2 and compute average F2 per Q^2 bin
        if filtered_points:
            q2_bins = {}
            for pt in filtered_points:
                if pt.q2 not in q2_bins:
                    # Check the deviation condition
                    for ptt in deviation_points:
                        if ptt[0] == pt.q2 and x < ptt[1]:
                            q2_bins[pt.q2] = []
                            break  # Found the matching deviation condition, break out of the loop
                if pt.q2 in q2_bins:
                    q2_bins[pt.q2].append(pt.F2)

            q2_vals = []
            avg_F2_vals = []

            for q2, F2_list in q2_bins.items():
                q2_vals.append(q2)
                avg_F2_vals.append(np.mean(F2_list))

            # Plot the results
            plt.scatter(q2_vals, avg_F2_vals, label=f'x={x:.3f}', color=colors[i])
            plt.plot(q2_vals, avg_F2_vals, color=colors[i])
            i += 1

    # Set plot labels and options
    plt.xlabel('$Q^2$')
    plt.ylabel('F2')
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    pdf.savefig()  # Save the figure to the PDF    
    #plt.show()
