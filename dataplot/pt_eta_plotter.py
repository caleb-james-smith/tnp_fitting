import matplotlib.pyplot as plt
import argparse
import os
import shutil
import numpy as np
import pandas as pd
from tabulate import tabulate

#parse arguments
parser = argparse.ArgumentParser(description='Plot Maker')
parser.add_argument("-f", "--filename", type=str, help="file name for reading")
parser.add_argument("-g", "--foldername", type=str, help="folder name for storage")
args=parser.parse_args()

print('Found argument list: ',args)
path = args.foldername

if os.path.exists(path):
    shutil.rmtree(path)
os.mkdir(path) #math folder

os.mkdir(path + '/pt')
os.mkdir(path + '/eta')

#read file
df = pd.read_csv(args.filename)

data_name = df.head()

#create dictionaries

d = {} #data efficiency arrays
MC = {} #MC efficiency arrays
errd = {} #data error arrays
errMC = {} #MC error arrays
SF = {} #scale factor arrays
errSF = {} #SF error arrays

scalefactor = []


#make some lists
pt = [] #pt bins


for series_name, series in df.items():

    if series_name.startswith('bin'):
        eta = series
    
    elif series_name.startswith('effdata'):
        d["{0}".format(series_name)] = np.array(series)
        pt.append(series_name[7:])

    elif series_name.startswith('errdata'):
        errd["{0}".format(series_name)] = np.array(series)

    elif series_name.startswith('effMC'):
        MC["{0}".format(series_name)] = np.array(series)

    elif series_name.startswith('errMC'):
        errMC["{0}".format(series_name)] = np.array(series)

    elif series_name.startswith('SF'):
        SF["{0}".format(series_name)] = np.array(series)

    elif series_name.startswith('errSF'):
        errSF["{0}".format(series_name)] = np.array(series)

pt = [i for n, i in enumerate(pt) if i not in pt[:n]]


print()


#graphing

for i in pt:

    fig1, ax1 = plt.subplots()
    plt.errorbar(eta, d[f'effdata{i}'], yerr = errd[f'errdata{i}'], fmt='o', color='blue', ecolor='blue', capsize=2, label='Data')
    plt.errorbar(eta, MC[f'effMC{i}'], yerr = errMC[f'errMC{i}'], fmt='s', color='green', ecolor='green', capsize=2, label='MC')

    plt.title(f'Efficiencies, pt: {i}')
    plt.xlabel('eta bins')
    plt.ylabel('efficiencies')
    plt.legend()
    plt.grid()
    plt.savefig(f'{path}/eta/{i}eff_eta.png')

    plt.close()

    fig2, ax2 = plt.subplots()
    plt.errorbar(eta, SF[f'SF{i}'], yerr = errSF[f'errSF{i}'], fmt='o', capsize=2)

    plt.title(f'Scale Factors, pt: {i}')
    plt.xlabel('eta bins')
    plt.ylabel('scale factors')
    plt.grid()
    plt.savefig(f'{path}/eta/{i}SF_eta.png')

    plt.close()

for i in range(len(eta)):

    efflistData =  [ d[f'effdata5-10'][i], d[f'effdata10-20'][i], d[f'effdata20-30'][i], d[f'effdata30-40'][i], d[f'effdata40-60'][i], d[f'effdata60-100'][i] ]
    errlistData = [ errd[f'errdata5-10'][i], errd[f'errdata10-20'][i], errd[f'errdata20-30'][i], errd[f'errdata30-40'][i], errd[f'errdata40-60'][i], errd[f'errdata60-100'][i] ]
    
    
    efflistMC = [ MC[f'effMC5-10'][i], MC[f'effMC10-20'][i], MC[f'effMC20-30'][i], MC[f'effMC30-40'][i], MC[f'effMC40-60'][i], MC[f'effMC60-100'][i] ]
    errlistMC = [ errMC[f'errMC5-10'][i], errMC[f'errMC10-20'][i], errMC[f'errMC20-30'][i], errMC[f'errMC30-40'][i], errMC[f'errMC40-60'][i], errMC[f'errMC60-100'][i] ]

    SFlist = [ SF[f'SF5-10'][i], SF[f'SF10-20'][i], SF[f'SF20-30'][i], SF[f'SF30-40'][i], SF[f'SF40-60'][i], SF[f'SF60-100'][i] ]
    errSFlist = [ errSF[f'errSF5-10'][i], errSF[f'errSF10-20'][i], errSF[f'errSF20-30'][i], errSF[f'errSF30-40'][i], errSF[f'errSF40-60'][i], errSF[f'errSF60-100'][i] ]

    for j in pt:
        scalefactor.append(SF[f'SF{j}'][i])
    
    fig3, ax3 = plt.subplots()
    plt.errorbar(pt, efflistData, yerr = errlistData, fmt='o', color='blue', ecolor='blue', capsize=2, label='Data')
    plt.errorbar(pt, efflistMC, yerr = errlistMC, fmt='s', color='green', ecolor='green', capsize=2, label='MC')

    plt.title(f'Efficiencies, eta: {eta[i]}')
    plt.xlabel('pt bins')
    plt.ylabel('efficiencies')
    plt.legend()
    plt.grid()
    plt.savefig(f'{path}/pt/{eta[i]}eff_pt.png')

    plt.close()

    fig4, ax4 = plt.subplots()
    plt.errorbar(pt, SFlist, yerr = errSFlist, fmt='o', capsize=2)

    plt.title(f'Scale Factors, eta: {eta[i]}')
    plt.xlabel('pt bins')
    plt.ylabel('scale factors')
    plt.grid()
    plt.savefig(f'{path}/pt/{eta[i]}SF_pt.png')

    plt.close()


#make visualization table
    
mydata = [
    [pt[5], f"{SF['SF60-100'][0]} +/- {errSF['errSF60-100'][0]}", f"{SF['SF60-100'][1]} +/- {errSF['errSF60-100'][1]}", f"{SF['SF60-100'][2]} +/- {errSF['errSF60-100'][2]}", f"{SF['SF60-100'][3]} +/- {errSF['errSF60-100'][3]}", f"{SF['SF60-100'][4]} +/- {errSF['errSF60-100'][4]}"],
    [pt[4],f"{SF['SF40-60'][0]} +/- {errSF['errSF40-60'][0]}", f"{SF['SF40-60'][1]} +/- {errSF['errSF40-60'][1]}", f"{SF['SF40-60'][2]} +/- {errSF['errSF40-60'][2]}", f"{SF['SF40-60'][3]} +/- {errSF['errSF40-60'][3]}", f"{SF['SF40-60'][4]} +/- {errSF['errSF40-60'][4]}"],
    [pt[3], f"{SF['SF30-40'][0]} +/- {errSF['errSF30-40'][0]}", f"{SF['SF30-40'][1]} +/- {errSF['errSF30-40'][1]}", f"{SF['SF30-40'][2]} +/- {errSF['errSF30-40'][2]}", f"{SF['SF30-40'][3]} +/- {errSF['errSF30-40'][3]}", f"{SF['SF30-40'][4]} +/- {errSF['errSF30-40'][4]}"],
    [pt[2], f"{SF['SF20-30'][0]} +/- {errSF['errSF20-30'][0]}", f"{SF['SF20-30'][1]} +/- {errSF['errSF20-30'][1]}", f"{SF['SF20-30'][2]} +/- {errSF['errSF20-30'][2]}", f"{SF['SF20-30'][3]} +/- {errSF['errSF20-30'][3]}", f"{SF['SF20-30'][4]} +/- {errSF['errSF20-30'][4]}"],
    [pt[1], f"{SF['SF10-20'][0]} +/- {errSF['errSF10-20'][0]}", f"{SF['SF10-20'][1]} +/- {errSF['errSF10-20'][1]}", f"{SF['SF10-20'][2]} +/- {errSF['errSF10-20'][2]}", f"{SF['SF10-20'][3]} +/- {errSF['errSF10-20'][3]}", f"{SF['SF10-20'][4]} +/- {errSF['errSF10-20'][4]}"],
    [pt[0], f"{SF['SF5-10'][0]} +/- {errSF['errSF5-10'][0]}", f"{SF['SF5-10'][1]} +/- {errSF['errSF5-10'][1]}", f"{SF['SF5-10'][2]} +/- {errSF['errSF5-10'][2]}", f"{SF['SF5-10'][3]} +/- {errSF['errSF5-10'][3]}", f"{SF['SF5-10'][4]} +/- {errSF['errSF5-10'][4]}"],
    ['X', eta[0], eta[1], eta[2], eta[3], eta[4]]
]

#print(len(scalefactor))
scalefactor.sort()

scalefactor = np.array(scalefactor)

onelist = []
for i in range(len(scalefactor)):
    onelist.append(1)
    
onelist = np.array(onelist)

SFdiffs = scalefactor - onelist

absSFdiffs = np.abs(SFdiffs)

absSFdiffs.sort()

diffmax = max(absSFdiffs)
diffmin = min(absSFdiffs)
diffrange = diffmax - diffmin
diffstep = diffrange/4

steplist = [ diffmin, diffmin + diffstep, diffmin + 2*diffstep, diffmin + 3*diffstep]#, diffmin + 4*diffstep, diffmin + 5*diffstep ]

#print data breakdown file

with open(f'{path}/results.txt', 'w') as f:

    for i in pt:

        print('pt:', i, file = f)
        print('------------', file = f)
        print(file = f)

        for j in range(5):

            print('eta:', eta[j], file =f)
            print(file = f)
            print('Data Efficiency:', d[f'effdata{i}'][j], '+/-', errd[f'errdata{i}'][j], file = f)
            print(file = f)
            print('MC Efficiency:', MC[f'effMC{i}'][j], '+/-', errMC[f'errMC{i}'][j], file = f)
            print(file = f)
            print('Scale Factor:', SF[f'SF{i}'][j], '+/-', errSF[f'errSF{i}'][j], file = f)
            print(file = f)
            print(file = f)

    print(tabulate(mydata, tablefmt='grid'), file = f)

    print(file = f)
    
    print(scalefactor, file = f)
    print(SFdiffs, file = f)
    print(absSFdiffs, file = f)

    print(file = f)

    print(steplist, file = f)
    print(diffmax, file = f)


