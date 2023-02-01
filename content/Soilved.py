#Useful functions in python for soil mechanics:
#Created by: Luis Zambrano-Cruzatty
#contact: luis.zambrano@maine.edu

import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

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

def Process_GSD(GSD_data, plot_fig=True):
    #GSD        is the GSD as a pandas dataframe
    Wtd=GSD_data[GSD_data.columns[1]].iloc[-1] #second column last row
    GSD_data['Weight retained']=[0.000000]*len(GSD_data[GSD_data.columns[1]]) # add new column for retained weight
    GSD_data['Retained']=[0.0000000]*len(GSD_data[GSD_data.columns[1]])#adds a new column named 'Retained'
    GSD_data['Passing']=[0.000000]*len(GSD_data[GSD_data.columns[1]])#adds a new column named 'Passing'
    W=0 #accumulates weight
    for i in range(len(GSD_data[GSD_data.columns[1]])):
        W += GSD_data[GSD_data.columns[1]][i] #accumulates weight
        GSD_data['Weight retained'][i]=W #writes accumulated weight to dataframe
        GSD_data['Retained'][i]=W*100/Wtd #calculates Percent coarser
        GSD_data['Passing'][i]=100-GSD_data['Retained'][i] #calculates percent finer

    #Plot details####
    if (plot_fig):
        PlotGSD(GSD_data[:-2])

def PlotGSD(GSD_data):
        plt.plot(GSD_data[GSD_data.columns[0]], GSD_data[GSD_data.columns[4]], '*-') #GSD plot
        plt.xscale('log') #d scale in log scale        
        plt.ylabel(r'Percent finer [$\%$]') # adds y label
        plt.xlabel(r'Particle size [mm]') # adds x label
        plt.xlim(0.0001,100)
        plt.ylim(0,100)
        plt.gca().invert_xaxis() #inverts d axis
        x_line=4.75
        y_line=0
        plt.axvline(x=4.75, color='red', linestyle='--')# gravel-san limit
        plt.annotate(f'Sieve No 4', xy=(x_line, y_line), xytext=(x_line + 10, y_line + 10),
            arrowprops=dict(facecolor='black', shrink=0.001))

        x_line=0.075
        y_line=90
        plt.axvline(x=0.075, color='red', linestyle='--')# sand-fines limit
        plt.annotate(f'Sieve No 200', xy=(x_line, y_line), xytext=(x_line , y_line ),
            arrowprops=dict(facecolor='black', shrink=0.001))
        plt.grid(b=True, which='major')# show major gridlines
        plt.grid(b=True, which='minor')# show minor grid lines

        ax2= plt.twinx() # second axis
        ax2.set_ylim(100,0)
        ax2.set_ylabel(r'Percent coarser [$\%$]')

def log_interp(x, Ds, P): #interpolates in GSD
    logd = np.log10(Ds) #transforms diameters to log values  
    interp_func=interp1d(P, logd) 
    return np.power(10.0, interp_func(x)) #finds solution

def GetClayFraction(D, Passing):
    d_clay=np.log10(0.002) #2 microns expressed in mm
    logD=np.log10(D) #diameter in log space
    iterp_func=interp1d(logD, Passing)
    return iterp_func(d_clay)


def GetSrinkageLimit(V_dry, Ws, gamma_w=9.81, Gs=2.67):
    #Vdry=      Dry volume in consistent units
    #Ws=        Dry weight in consistent units
    #gamma_w=   Unit weight of water in consistent units default 9.81 kN/m^3
    #Gs=        Specific gravity default 2.67

    return ((V_dry*gamma_w/Ws)-(1/Gs))*100 


def GetUSCS(GF, SF, FF, Cu=1000, Cc=10000, is_organic=False,PI=0, LL=0):
    #GF=        Gravel fraction in %
    #SF=        Sand fraction in %
    #FF=        Fines fraction in %
    #Cu=        Coefficient of uniformity default is high meaning well graded
    #Cc=        Coefficient of curvature default for well graded mat.
    #PI=        Plasticity index in % default 0
    #LL=        Liquid limit, default 0
    #Is_organic Flag to indicate material is organic
    prefix=""
    sufix=""
    if (FF>= 50): #material is clay, silt, or organic
        if (LL<50): #low plasticity
            if (is_organic): #material is organic
                Group="OL"
                Name="Low-plastic organic soil"
            else: #material is not organic
                PIA=0.73*(LL-20) #Plasticity index on the A-line                
                if ((PI>7) and (PI>=PIA)):#CL
                    Group="CL"
                    prefix, sufix=GetDetailedName_fines(FF, SF, GF)
                    if (prefix==""):
                        Name="Lean clay %s" %sufix
                    else:
                        Name="%s lean clay %s" %(prefix, sufix)
                elif((PI>=4) and (PI<=7) and (PI>=PIA)): #CL-ML
                    Group="CL-ML"
                    prefix, sufix=GetDetailedName_fines(FF, SF, GF)
                    if (prefix==""):
                        Name="Silty clay %s" %sufix
                    else:
                        Name="%s silty clay %s" %(prefix, sufix)
                else: #ML
                    Group="ML"
                    prefix, sufix=GetDetailedName_fines(FF, SF, GF)
                    if (prefix==""):
                        Name="Silt %s" %sufix
                    else:
                        Name="%s silt %s" %(prefix, sufix)


        else: #high plasticity
            if (is_organic): #material is organic
                Group="OH"
                Name="High-plastic organic soil"
            else: #material is not organic
                PIA=0.73*(LL-20) #Plasticity index on the A-line                
                if (PI>=PIA):#CH
                    Group="CL"
                    prefix, sufix=GetDetailedName_fines(FF, SF, GF)
                    if (prefix==""):
                        Name="Fat clay %s" %sufix
                    else:
                        Name="%s fat clay %s" %(prefix, sufix)
                else: #MH
                    Group="MH"
                    prefix, sufix=GetDetailedName_fines(FF, SF, GF)
                    if (prefix==""):
                        Name="Elastic silt %s" %sufix
                    else:
                        Name="%s elastic silt %s" %(prefix, sufix)

    else: #material is coarse (sand or gravel)
        if (GF>=SF): #Gravel
            if (FF<5):
                if ((Cu>=4) and (Cc>=1) and (Cc<=3)): #GW
                    Group="GW"
                    if (SF<15):
                        Name="Well-graded gravel"
                    else:
                        Name="Well-graded gravel with sand"
                else: #GP
                    Group="GP"
                    if (SF<15):
                        Name="Poorly-graded gravel"
                    else:
                        Name="Poorly-graded gravel with sand"
            elif (FF>12):
                Group, _, prefix=GetMinorFraction_coarse(PI, LL, "G")
                if (SF>=15):
                    Name="%s gravel with sand" %prefix
                else:
                    Name="%s gravel" %prefix
            else:
                if ((Cu>=4) and (Cc>=1) and (Cc<=3)): #GW
                    Group="GW"
                    Name="Well-graded gravel"
                else: #GP
                    Group="GP"
                    Name="Poorly-graded gravel"
                second_group, sufix, _=GetMinorFraction_coarse(PI, LL, "G")
                Group= "%s-%s" %(Group, second_group)
                Name= "%s %s" %(Name, sufix)
                if (SF>=15):
                    Name="%s and sand" %Name
        else: #sand
            if (FF<5):
                if ((Cu>=4) and (Cc>=1) and (Cc<=3)): #SW
                    Group="SW"
                    if (GF<15):
                        Name="Well-graded sand"
                    else:
                        Name="Well-graded sand with gravel"
                else: #SP
                    Group="SP"
                    if (GF<15):
                        Name="Poorly-graded sand"
                    else:
                        Name="Poorly-graded sand with gravel"
            elif (FF>12):
                Group, _, prefix=GetMinorFraction_coarse(PI, LL, "S")
                if (GF>=15):
                    Name="%s sand with gravel" %prefix
                else:
                    Name="%s sand" %prefix
            else:
                if ((Cu>=4) and (Cc>=1) and (Cc<=3)): #SW
                    Group="SW"
                    Name="Well-graded sand"
                else: #sP
                    Group="SP"
                    Name="Poorly-graded sand"
                second_group, sufix, _=GetMinorFraction_coarse(PI, LL, "S")
                Group= "%s-%s" %(Group, second_group)
                Name= "%s %s" %(Name, sufix)
                if (GF>=15):
                    Name="%s and gravel" %Name
    return Group, Name

def GetDetailedName_fines(FF, SF, GF):
    Percent_coarser=100-FF #percent coarser
    sufix=""
    prefix=""
    if (Percent_coarser<30):
        if (Percent_coarser<15):
            sufix=""
        elif(SF>=GF):
            sufix="with sand"
        else:
            sufix="with gravel"
    else:
        if (SF>GF):
            if (GF<15):
                prefix="Sandy"
            else:
                prefix="Sandy"
                sufix="with gravel"
        else:
            if (SF<15):
                prefix="Gravelly"
            else:
                prefix="Gravelly"
                sufix="with sand"
    return prefix, sufix

def GetMinorFraction_coarse(PI, LL, FL):
    PIA=0.73*(LL-20) #Plasticity index on the A-line 
    sufix=""
    prefix=""
    if (LL<50): #low plasticity                      
        if ((PI>7) and (PI>=PIA)):#C
            Group="%sC" %FL
            sufix="with clay"
            prefix="Clayey"
        elif((PI>=4) and (PI<=7) and (PI>=PIA)): #C-M
            Group="%sC-%sM" %(FL, FL)
            sufix="with silty clay"
            prefix="Silty, clayey"
        else: #M
            Group="%sM" %FL
            sufix="with silt"
            prefix="Silty"
    else: #high plasticity              
        if (PI>=PIA):#C
            Group="%sC" %FL
            sufix="with clay"
            prefix="Clayey"
        else: #M
            Group="%sM" %FL
            sufix="with silt"
            prefix="Silty"
    return Group, sufix, prefix







        
