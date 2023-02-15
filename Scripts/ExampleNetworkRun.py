# -*- coding: utf-8 -*-
"""
Created on Sun Jan 15 11:27:23 2023

@author: dgreen1
"""

#Example script for using the DepressionalNetwork class to analyze a network
#This script analyzes and reconstructs the example network plot given in Figure X of Green and Crumpton (2023)

import DepressionalNetworkv2142023 as dn
import numpy as np
import pandas as pd
import sys
import time
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
sys.setrecursionlimit(4155)

##Specify the locations of the input and output folders
inpath = r'E:\Backups\Backup\DepressionalCascade\Github\Inputs'
outpath = r'E:\Backups\Backup\DepressionalCascade\Github\Outputs'

##Read the depressional attribute and network files into Pandas dataframes
##The depression data file should have the form DepressionID|MaximumArea|MaximumVolume|MinimumContributingArea|CurveNumber
##MaximumDepth is not used in the algorithm, but can be estimated as the quotient of the maximum volume and maximum area 
##The network data file should have the form DepressionID|CatchmentID with a one-to-many relationship from DepressionID to CatchmentID
dep_data = pd.read_table(inpath + r'\ExampleNetworkDepressionData.txt',sep='\t')
network_data = pd.read_table(inpath + r'\ExampleNetworkJoinTable.txt',sep='\t')

##Sort the Pandas dataframes by DepressionID. This step is not necessary, but speeds up subsequent searches
dep_data.sort_values('DepressionID')
network_data.sort_values('DepressionID')

##Call the ConstructUpslopeSets method in the DepressionNetwork code. 
##This will construct the set of directly neighboring upstream contributing depressions for
##each depression in the network 
dn.ConstructUpstreamSets(dep_data,network_data)

##Call the ConstructDownstreamSets method in the DepressionNetwork code. 
##This will construct the set of directly neighboring downslope contributing depressions for
##each depression in the network 
dn.ConstructDownstreamSets(dep_data,network_data)

##Check for circular references in the network data (i.e. a depression does not have a flow connection to one or more outlet depressions)
dn.UpdateTypesAndFlags(dep_data)
dn.CheckNetworkStructure(dep_data)

##Create a new network named 'Example'
network = dn.DepressionalNetwork('Example')
network.construct_network(dep_data)

##Retreive the total storage volume and total contributing areas of the network (sum(Vmax) and sum(Amax))
Vt = network.total_storage_volume
At = network.total_contributing_area
##Retreive the size of the network
n = len(network.network_list)

##Calculate the network depressional specific storage index
Sd = Vt/At

##Create an empty list to append the IDs of cascading depressions
##This was used to help create figures 9, 10, and 13
casc = []

##Create an empty list to append record data for each depression
##as shown in the first for loop. This will be used to fill-in data tables for each 
##precipitation-dependent property listed in Green and Crumpton, 2023
records = []

##Create a dictionary of each depression and its primary characteristics
for dep in network.network_list:
    record = {'DepressionID':int(dep.identifier),'Am':dep.maximum_area,'Vm':dep.maximum_volume,'Ac':dep.local_contributing_area,'Nus':len(dep.upstream_depressions),'Nds':len(dep.downstream_depressions),'Nuscd':len(dep.upslope_network),'Ndsds':len(dep.downslope_network),'USAc':dep.upslope_contributing_area,'USVm':dep.upslope_maximum_storage_volume,'OutletID':[d.identifier for d in dep.outlet_depression],'Sdus':dep.upslope_specific_storage,'Strn':dep.strahler_number,'Shrn':dep.shreve_number}
    records.append(record)

##Creater unique dataframes for each precipitation-dependent property
#depressions = pd.DataFrame.from_records(records)
volume = pd.DataFrame.from_records(records)
export = pd.DataFrame.from_records(records)
contarea = pd.DataFrame.from_records(records)
usvol = pd.DataFrame.from_records(records)
numcusd = pd.DataFrame.from_records(records)
usexport = pd.DataFrame.from_records(records)
locrunoff = pd.DataFrame.from_records(records)

##Set the index of each dataframe to the depression ID. This allows for quick and easy querying.
volume = volume.set_index('DepressionID')
export = export.set_index('DepressionID')
contarea = contarea.set_index('DepressionID')
usvol = usvol.set_index('DepressionID')
numcusd = numcusd.set_index('DepressionID')
usexport = usexport.set_index('DepressionID')
locrunoff = locrunoff.set_index('DepressionID')

##Create a set of only outlet depressions. This will be used to "solve" each network for each precipitation amount
##The example network contains a single outlet, and so the algorithm will give one solution for each precipitation amount
##Since this code uses object-oriented structures to represent a network, each depression within a network and the entire network can be queried for network or depression-specific properties 
outlets = network.outlet_depressions

#Create an array of precipitation amounts to iterate over. Each value in the list will be fed into the routing routine to product a new result
P = np.arange(0.0,0.6,0.01)
P[0] = 0.001

#Create lists to capture outputs for plotting
Qr = []
Ar = []
Nr = []
Nf = []
Vs = []
Qt = []
Dr = []
Fr = []
Ff = []
FAr = []
FQr = []
Fv = []
Fs = []
Rs = []
Nhw = []
Nms = []
Nout = []
Niso = []

#Lists to store function call times
SolTime = []
CasTime = []
RecTime = []
TotTime = []

solution_sets = {}

k = 0
s1 = time.time()

for dep in network.network_list:
    dep.local_curvenumber = 100.0

for p in P:
    print('Rainfall: ' + str(p))
    start = time.time()
    
#Assign uniform precipitation to each depression
    for depression in outlets:
        depression.local_precipitation = p
    
#Time the solution
    mpt1 = time.time()
    solution = network.route_runoff_uniform_precipitation(p)
    mpt2 = time.time()
    SolTime.append(mpt2-mpt1)
    print('Solution Elapsed Time: ' + str(mpt2-mpt1))
    
    Qr.append(solution[0])
    Ar.append(solution[1])
    Nr.append(solution[2])
    Nf.append(solution[3])
    Vs.append(solution[4])
    Qt.append(solution[5])
    Dr.append({'P': p,'SolutionSet':solution[6]})
    
    Fr.append(solution[2]/float(n))
    Ff.append(solution[3]/float(n))
    FAr.append(solution[1]/At)
    FQr.append(solution[0]/solution[5] if solution[5] > 0 else 0)
    Fv.append(solution[4]/Vt)
    Fs.append(solution[4]/solution[5] if solution[5] > 0 else 0)
    
    Rs.append(Vt/solution[5] if solution[5] > 0 else 0)
    
    mpt3 = time.time()
    print('Solution Elapsed Time: ' + str(mpt3 - mpt2))

#This subsection illustrates how to use the script to tabulate the types of depressions contributing runoff to the outlet depression
    h = 0
    hw = 0
    ms = 0
    out = 0
    iso = 0
    for dep in solution[6]:
        record = {'P':p,'DepressionID':int(dep.identifier)}
        casc.append(record)
        
        if dep.depression_type == 'Headwater':
            hw += 1
        elif dep.depression_type == 'Interior':
            ms += 1
        elif dep.depression_type == 'Outlet':
            out += 1
        elif dep.depression_type == 'Disjunct':
            iso += 1
        
    Nhw.append(hw)
    Nms.append(ms)
    Nout.append(out)
    Niso.append(iso)
    
    
    mpt4 = time.time()
    print('Cascade Record Compilation Time: ' + str(mpt4 - mpt3))
    CasTime.append(mpt4-mpt3)
    
    nwlist = network.network_list
    
#Assign solution values to the solution values
    for dep in nwlist:
        depid = int(dep.identifier)
        
        volume.at[depid,'p'+str(p)] = dep.storage_volume
        export.at[depid,'p'+str(p)] = dep.export
        contarea.at[depid,'p'+str(p)] = dep.runoff_contributing_area
        numcusd.at[depid,'p'+str(p)] = len(dep.runoff_network)
        usvol.at[depid,'p'+str(p)] = dep.upslope_storage_volume
        usexport.at[depid,'p'+str(p)] = dep._upstream_export
        locrunoff.at[depid,'p'+str(p)]=dep.local_runoff
        
    k += 1
    stop = time.time()
    print('Depression Record Compilation: ' + str(stop - mpt4))
    print('Total Time elapsed: ' + str(stop - start))
    RecTime.append(stop - mpt4)
    TotTime.append(stop-start)

s2 = time.time()
print('Total Elapsed Time: '+str(s2-s1))
cascade = pd.DataFrame.from_records(casc)

#write results to tables
usvol.to_csv(outpath+r'\usvolume.txt',sep='\t')
volume.to_csv(outpath+r'\depvolume.txt',sep='\t')
export.to_csv(outpath+r'\depexport.txt',sep='\t')
contarea.to_csv(outpath+r'\contarea.txt',sep='\t')
numcusd.to_csv(outpath+r'\numcusd.txt',sep='\t')
usexport.to_csv(outpath+r'\usexport.txt',sep='\t')
locrunoff.to_csv(outpath+r'\localrunoff.txt',sep='\t')
cascade.to_csv(outpath+r'\cascade.txt',sep='\t')

#append solution to a dataframe and write to the specified output directory
solutions = pd.DataFrame({'P':P,'Qr':Qr,'Ar':Ar,'Nr':Nr,'Nf':Nf,'Vs':Vs,'Qt':Qt,'Fr':Fr,'Ff':Ff,'FAr':FAr,'FQr':FQr,'Fv':Fv,'Fs':Fs,'SolTime':SolTime,'CasTime':CasTime,'RecTime':RecTime,'TotTime':TotTime,'Nhw':Nhw,'Nms':Nms,'Nout':Nout,'Niso':Niso})
solutions.to_csv(outpath+r'\modeloutputs.txt',sep='\t')

#plot results
fig,axs = plt.subplots(2,1,sharex=True)
fig.set_size_inches(8,8,forward=True)
fig.subplots_adjust(hspace = 0.1, wspace = 0.0)
ax0tw = axs[0].twinx()
ln2 = axs[0].plot(P*1000.0,Qt,ls='--',color='black',lw=1,label=r'$Q_{tot}$')
ln1 = axs[0].plot(P*1000.0,Qr,ls='-',lw=1,color='black',label=r'$Q_{ro}$')
ln3 = axs[0].plot(P*1000.0,Vs,ls='-.',color='black',lw=1,label=r'$V_{sto}$')
ln4 = ax0tw.plot(P*1000.0,[A*0.0001 for A in Ar],ls=':',color='black',lw=1,label=r'$A_{ro}$')
axs[0].set_xlim([0,0.5*1000.0])
axs[0].set_ylim([0,1.0e7])
axs[0].set_ylabel(r'Volume $(\mathrm{m^3})$',fontsize=12)
ax0tw.set_ylabel(r'Area (ha)',fontsize=12)
axs[0].yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
ax0tw.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.1e'))
axs0 = ln1+ln2+ln3+ln4
labs0 = [ax.get_label() for ax in axs0]
leg = axs[0].legend(axs0,labs0,loc='upper left',fontsize=12)
leg.get_frame().set_alpha(0.9)

ff = axs[1].plot(P*1000.0,Ff,ls='-.',color='black',lw=1,label=r'$f_{fill}$')
fr = axs[1].plot(P*1000.0,Fr,ls='-',color='black',lw=1,label=r'$f_{ro}$')
far = axs[1].plot(P*1000.0,FAr,'--',color='black',lw=1,label=r'$f_{A_{ro}}$')
frs = axs[1].plot(P*1000.0,Fs,ls=(0, (3, 3, 10, 10)),color='black',lw=1,label=r'$f_{sto}$')
fQr = axs[1].plot(P*1000.0,FQr,':',color='black',lw=1,label=r'$f_{Q_{ro}}$')
axes1 = ff+fr+far+frs+fQr
labs = [ax.get_label() for ax in axes1]
axs[1].set_xlabel('Precipitation (mm)',fontsize=12)
axs[1].set_ylabel('Value (-)',fontsize=12)
axs[1].set_ylim([0,1.0])
leg2 = axs[1].legend(labs,loc='lower right',fontsize=12)
leg2.get_frame().set_alpha(0.9)
fig.tight_layout()

plt.savefig(outpath+'\Figure3.tif',dpi=300)
plt.savefig(outpath+'\Figure3.pdf',dpi=600)
