#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 22 22:05:35 2019

Description: Generates panel (e) of main Fig. 2, i.e. the 
             Gaussian Process extrapolation results for the
             static linear decoder, for the PPC manuscript.
             Unlike Fig. 2d, for this analysis, the "best"
             subset of size K neurons was chosen based on 
             ranking based on explained variance. 
Note:        The data file (referred to here as "filename") 
             was generated in Matlab, using the function
             MAEvsN_static_ranked.m, which calls on the 
             function rankedNeurons_explainedVar.m
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

    X = X.ravel()
    mu = mu.ravel()
    uncertainty = 1.96 * np.sqrt(np.diag(cov))
    
    ax = plt.subplot(3,3,((fidx-1)*3)+1) 
    # ax = plt.subplot(2,2,fidx) 
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
    plt.ylim(ylimits[0], ylimits[-1])  
    plt.xlabel(xtitle,fontsize=15)
    plt.ylabel(ytitle,fontsize=15) 

def data_loader(filename):
    df = pd.read_csv(filename, ',', header=None)
    X_train  = df.get(0).values # N values 
    X_train  = X_train.reshape(-1, 1)
    return df, X_train

if __name__ == '__main__':
    figure(num=None, figsize=(7, 7), dpi=80)
    
    # Load data 
    filename    = "Results/MAEvsN_static_ranked_m04_5.csv"
    df, X_train = data_loader(filename)
    
    # pdb.set_trace()
    
    # Define extrapolation points 
    X_extrap = np.array(range(1, 401, 1))
    X_extrap = X_extrap.reshape(-1, 1)
    
    for i in range(0,3):
        
        Y_train = df.get(i+1).values # mean MAE values 
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
            ylimits = [0.2, 0.9] #range(0,2,1)
        elif i==1:
            xtitle  = " " 
            ytitle  = "MAE (m/s)"
            ylimits = [0, 0.2]
        else:
            xtitle  = "# of PPC Neurons"
            ytitle  = "MAE (deg)"
            ylimits = [7, 14]
        
        # Plot
        plot_gp(i+1, mu_s, cov_s, X_extrap, X_train=X_train, Y_train=Y_train, 
                xtitle=xtitle, ytitle=ytitle,xlimits=[1, 200, 400],
                ylimits=ylimits)
        
        # Hack from MRULE to get axes to match
        if i==0:
            plt.ylim(0,1)
        if i==1:
            plt.ylim(0,0.15)
        if i==2:
            plt.ylim(0,25)
        plt.savefig("Fig2d.svg")
