#Useful functions in python for soil mechanics:
#Created by: Luis Zambrano-Cruzatty
#contact: luis.zambrano@maine.edu

import numpy as np 
import pandas as np

def GetWaterContent(WtPWc, WsPWc,Wc):
    #WtPWc=         Total weight + weight of container
    #WsPWc=         Weight of solids + weight of container
    #Wc=            Weight of container

    Ww=WtPWc-WsPWc
    Ws=WsPWc-Wc

    return (Ww/Ws)*100 #water content in percentage

def GetTotalUWeight(Gs, S, e, gamma_w=9.81):
    #Gs=            specific gravity
    #S=             saturation dimensionless
    #e=             void ratio
    #gamma_w        U. weight of water SI default

    return gamma_w*(Gs+S*e)/(1+e) #returns UW in same units as gamma_w

def Process_GSD(GSD_data):
    #GSD        is the GSD as a pandas dataframe
    Wtd=GSD_data[GSD_data.columns[1]].iloc[-1] #second column last row
    GSD_data['Weight retained']=[0]*len(GSD_data[GSD_data.columns[1]]) # add new column for retained weight
    GSD_data['Retained']=[0]*len(GSD_data[GSD_data.columns[1]])#adds a new column named 'Retained'
    GSD_data['Passing']=[0]*len(GSD_data[GSD_data.columns[1]])#adds a new column named 'Retained'
    W=0 #accumulates weight
    for i in range(len(GSD_data[GSD_data.columns[1]])):
        W += GSD_data[GSD_data.columns[1]][i] #accumulates weight
        GSD_data['Weight retained'][i]=W*100/Wtd
        GSD_data['Retained'][i]=W*100/Wtd
        GSD_data['Passing'][i]=100-GSD_data['Retained'][i]
        
