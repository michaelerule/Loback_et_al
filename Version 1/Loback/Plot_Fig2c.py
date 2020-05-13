#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  1 14:31:58 2019

Description: Generates panel (d) of main Fig. 2, i.e. the 
             Gaussian Process extrapolation results for the
             static linear decoder, for the PPC manuscript.
Note:        The data file (referred to here as "filename") 
             was generated in Matlab, using the function
             MAEvsN_static.m
             See the /Adrianna PPC ms/Code/Paper_fns subfolder 
             located in the HFSP shared Dropbox for this Matlab function.

@author: adrianna
"""

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import ConstantKernel
from sklearn.gaussian_process.kernels import WhiteKernel, Matern
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from matplotlib.pyplot import figure
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Tahoma']

# Define color scheme
tableau20 = np.float32([\
    ( 31, 119, 180), (174, 199, 232), (255, 127,  14), 
    (255, 187, 120), ( 44, 160,  44), (152, 223, 138), 
    (214,  39,  40), (255, 152, 150), (148, 103, 189), 
    (197, 176, 213), (140,  86,  75), (196, 156, 148),    
    (227, 119, 194), (247, 182, 210), (127, 127, 127), 
    (199, 199, 199), (188, 189,  34), (219, 219, 141), 
    ( 23, 190, 207), (158, 218, 229)])/255.        

def plot_gp(fidx, mu, cov, X, X_train=None, Y_train=None, samples=[], 
            xtitle="", ytitle="", xlimits=range(0, 500, 200), 
            ylimits=range(0, 1)):
    
    X  = X.ravel()
    mu = mu.ravel()
    uncertainty = 1.96 * np.sqrt(np.diag(cov))
    
    ax = plt.subplot(3,3,((fidx-1)*3)+1) 
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False) 
    
    plt.fill_between(X, mu + uncertainty, mu - uncertainty, 
                     color=tableau20[fidx*5-1], alpha=0.2)
    plt.plot(X, mu, color="black", lw=2)
    
    for i, sample in enumerate(samples):
        plt.plot(X, sample, lw=1, ls='--', label=f'Sample {i+1}')
    
    if X_train is not None:
        plt.plot(X_train, Y_train, 'x', color=tableau20[6])
    
    # plt.legend()
    plt.yticks(ylimits, fontsize=14)  
    plt.xticks(xlimits, fontsize=14)
    plt.ylim(0, ylimits[-1])  
    plt.xlabel(xtitle,fontsize=15)
    plt.ylabel(ytitle,fontsize=15) 
    plt.savefig("Fig2c.svg")

def data_loader(filename):
    df = pd.read_csv(filename, ',', header=None)
    X_train  = df.get(0).values # N values 
    X_train  = X_train.reshape(-1, 1)
    return df, X_train

if __name__ == '__main__':
    figure(num=None, figsize=(7, 7), dpi=80)
    # Load data 
    filename    = "./Results/MAEvsN_static_m04_Nthresh200_5.csv"
    df, X_train = data_loader(filename)
    # Define extrapolation points 
    X_extrap = np.array(range(1, 401, 1))
    X_extrap = X_extrap.reshape(-1, 1)
    for i in range(0,3):
        Y_train = df.get(i+1).values #mean MAE values 
        Y_train = Y_train.reshape(-1, 1)
        # Define kernel for GP
        k = ConstantKernel(1.0) * Matern(length_scale=1.0) + WhiteKernel(1.0)
        gp = GaussianProcessRegressor(kernel=k)
        gp.fit(X_train, Y_train)
        # Compute posterior predictive mean and covariance
        mu_s, cov_s = gp.predict(X_extrap, return_cov=True)
        # Show optimized kernel parameters
        print(gp.kernel_)
        if i==0:
            xtitle  = " "
            ytitle  = "MAE (m)" 
            ylimits = [0, 0.5, 1.0] #range(0,2,1)
        elif i==1:
            xtitle  = " " 
            ytitle  = "MAE (m/s)"
            ylimits = [0, 0.2]
        else:
            xtitle  = "# of PPC Neurons"
            ytitle  = "MAE (deg)"
            ylimits = [0, 25]
        # Plot
        plot_gp(i+1, mu_s, cov_s, X_extrap, X_train=X_train, Y_train=Y_train, 
                xtitle=xtitle, ytitle=ytitle,xlimits=[1, 200, 400],
                ylimits=ylimits)
