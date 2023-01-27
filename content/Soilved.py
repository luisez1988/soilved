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
