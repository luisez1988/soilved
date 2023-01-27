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
