import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import pandas as pd
import math
import sys

def gaussian(x, A, x0, std):
    '''
        Function:
            Create a gaussian curve

        Parameters:
            x (np.array(float32)): x data
            A (float32): peak of curve
            x0 (float32): center of the curve
            std (float32): standard deviation

        Returns:
            gaussianCurve (np.array(float32)): Gaussian Curve
    '''

    return A * np.exp((-(x - x0) ** 2) / (2 * (std ** 2)))

def gaussianFit(x_range, y_data, x0, a, std, plot=False):
    '''
        Function:
            Fits and plots a curve to a gaussian wave

        Parameters:
            x_range (np.array(float32)): the range of x_values
            x0 (float): x_value (mass) at the peak,
            a (float): y_value (intensity) at the peak
            plot (bool): whether or not you want to plot the curve

        Returns:
            x, y (np.array(float32), np.array(float32)): Gaussian Curve
    '''
    try:
        params, _ = curve_fit(gaussian, x_range, y_data, p0=[a, x0, std])
        fitted_a, fitted_x0, fitted_sigma = params
        fitted_y_data = gaussian(x_range, fitted_a, fitted_x0, fitted_sigma)

        if plot:
            plt.plot(x_range, fitted_y_data)
            plt.show()
    except Exception as e:
        return None

    return x_range, fitted_y_data, params

def multiGaussian(x, *params):
    '''
        Function:
            Gaussian function for multiple peaks

        Parameters:
            x (np.array(np.float32)): the x data
            params (list(float)): the "initial guesses" for the parameters of each curve involved in the multifit
                -> [a(0), x0(0), std(0), a(1), x0(1), std(1)...a(n), x0(n), std(n)]

        Return:
            curve (list(float)): y-data of the gaussian curve
    '''
    length = len(params) // 3
    result = np.zeros_like(x)
    for i in range(length):
        result += gaussian(x, params[i * 3], params[i * 3 + 1], params[i * 3 + 2])
    return result

def multiPeakGaussianFit(x, y, initial_guesses):
    '''
        Function:
            Fit multiple peaks whose singular fit curves intersect

        Parameters:
            x (np.array(np.float32)): range of m/z values on the actual mass to intensity graph
            y (np.array(np.float32)): the y-data corresponding to the x values
            initial_guesses (list(float)): the "initial guesses" for the parameters of each curve involved in the multifit
                -> [a(0), x0(0), std(0), a(1), x0(1), std(1)...a(n), x0(n), std(n)]
        
        Return:
            x, fitted_y, x0 (np.array(np.float32), np.array(np.float32), np.array(np.float32)):
                full extent of the x range,
                fit intensity,
                all fit x0 values
        
    '''
    try:
        params, _ = curve_fit(multiGaussian, x, y, initial_guesses)
        fitted_y = multiGaussian(x, *params)

        return x, fitted_y, params[1::3]
    except Exception as e:
        return None

def find_peaks(mz_b, mz_av, output_path, plot=False):
    '''
        Function:
            Find peaks on an mz_base to mz_av function and writes it to a file
        
        Parameters:
            mz_b (str): path to m/z data
            mz_av (str): path to intensity data
            outputPath (str): path to the file to write the x0 of the peaks in
            plot (bool): whether or not to plot the (mz_base, mz_av) graph
    '''
    mz_b = np.array(pd.read_csv(mz_b)['mz_base'])
    mz_av = np.array(pd.read_csv(mz_av)["mz_av"])

    dy = np.diff(mz_av)

    peaks = []
    lows = []
    highs = []
    stds = []
    gaussian_curves = []
    curve_buffer = []
    final_x0 = []
    
    tolerance = 1000

    for i in range(len(dy) - 1):
        if mz_av[i] < 1500:
            continue
        
        if dy[i] >= 0 and dy[i + 1] < 0:
            peak = i + 1
            
            low = i
            while (low > 0 and not (dy[low] >= 0 and dy[low - 1] < 0)):
                low -= 1

            high = i + 2
            while (high < len(dy) and not (dy[high] >= 0 and dy[high - 1] < 0)):
                high += 1

            domain = mz_b[low : high + 1]
            y_data = mz_av[low : high + 1]
            std = (domain[-1] - domain[0]) / (2 * math.sqrt(2 * math.log(2)))

            if len(domain) < 3:
                continue

            
            if (lows and highs and mz_b[low] <= mz_b[highs[-1]]):
                if len(curve_buffer) == 0:
                    curve_buffer.append(len(lows) - 1)
                    gaussian_curves.pop()
                curve_buffer.append(len(lows))
            else:
                if (curve_buffer):
                    gc = []
                    initial_guesses = []
                    for c in curve_buffer:
                        d = mz_b[ lows[c] : highs[c] ]
                        y = mz_av[ lows[c] : highs[c] ]
                        s = (d[-1] - d[0]) / (2 * math.sqrt(2 * math.log(2)))
                        curve = gaussianFit(d, y, mz_b[peaks[c]], mz_av[peaks[c]], s)
                        if curve:
                            x, y, params = curve
                            gc.append([x, y])
                            for j in params:
                                initial_guesses.append(j)
                        else:
                            gc.append([mz_b[ lows[c] : highs[c] ], mz_av[ lows[c] : highs[c] ]])
                            initial_guesses.append(mz_av[peaks[c]])
                            initial_guesses.append(mz_b[peaks[c]])
                            initial_guesses.append(stds[c])

                    
                    x_multi = []
                    y_multi = []
                    for j in gc:
                        for k in range(len(j[0])):
                            x_multi.append(j[0][k])
                            y_multi.append(j[1][k])
                                
                    curve_buffer = []

                    multifit = multiPeakGaussianFit(x_multi, y_multi, initial_guesses)
                    if not multifit:
                        continue

                    gaussian_curves.append(multifit[:2])
                    for mz_val in multifit[2]:
                        final_x0.append(mz_val)
                    continue;
                    
            
                curve = gaussianFit(domain, y_data, mz_b[peak], mz_av[peak], std)
                if not curve:
                    continue
                gaussian_curves.append(curve[:2])
            
            peaks.append(peak)
            lows.append(low)
            highs.append(high)
            stds.append(std)

    with open(output_path, "w+") as fw:
        for i in final_x0:
            fw.write(str(i) + "\n")

        fw.close()
    
    if plot:
        fig, ax = plt.subplots(1, 1)
        ax.plot(mz_b, mz_av, color='yellow')
        ax.set_title('Mass Spectra', fontsize=14)
        ax.set_xlabel('m/z', fontsize=12)
        ax.set_ylabel('Intensity', fontsize=12)
        for c in gaussian_curves:
            plt.plot(c[0], c[1], color="red")
            
        plt.show()


find_peaks(sys.argv[1],
           sys.argv[2],
           sys.argv[3])
