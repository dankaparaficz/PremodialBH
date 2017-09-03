#!/usr/bin/env python
########################################################################################################################
########################################################################################################################
#Creates Ligtcurves and performes test if the simulated light cuvres vary the same
# way as data

#CALLING:
#python search_plot_many_convolved.py input 'quasar_name' 'data_file'
#INPUT :

#input - from 1-11 -
#1: 100% dark matter in form of compact objects,
#2: 90% dark matter in form of compact objects,
#3: 80% dark matter in form of compact objects,
#4: 70% dark matter in form of compact objects,
#5: 60% dark matter in form of compact objects,
#6: 50% dark matter in form of compact objects,
#7: 40% dark matter in form of compact objects,
#8: 30% dark matter in form of compact objects,
#9: 20% dark matter in form of compact objects,
#10: 10% dark matter in form of compact objects,
#11: 1% dark matter in form of compact objects

#quasar_name - dicrectory where simulated quasar micorlesning data are e.g. 'HE0435_B_new'

#data_file - real datapoint of quasar lightcurve 'HE0435_all_astrofix_2016'

########################################################################################################################
########################################################################################################################

import numpy as np
import os.path
import math
import sys
import matplotlib.pyplot as plt
import scipy.ndimage as ndimage
import matplotlib.gridspec as gridspec
import argparse
import time
import sizes


parser = argparse.ArgumentParser()
parser.add_argument("input", help="Dark matter percentage.", type=int)
parser.add_argument("quasar_name", help="Name of quasar system.", type=str)
parser.add_argument("data_file", help="File where you strone light curve.", type=str)
parser.add_argument("--conv", help="Do you want convolution?", action='store_true')
parser.add_argument("--speed", help="What speed do you want?", type=int)
args = parser.parse_args()


########################################################################################################################
# List of statysitical methods for
########################################################################################################################
def chi2(a, b):
    return np.sum((a-b) * (a-b))
def norm_a_chi2(a, b):
    return np.sum((a-b) * (a-b) / a)
def norm_b_chi2(a, b):
    return np.sum((a-b) * (a-b) / b)
def l1(a, b):
    return  np.sum(np.abs((a-b)))
def std_dev(a):
    return np.std(a)
def mad(arr):
    """ Median Absolute Deviation: a "Robust" version of standard deviation.
        Indices variabililty of the sample.
        https://en.wikipedia.org/wiki/Median_absolute_deviation 
    """
    arr = np.ma.array(arr).compressed() # should be faster to not use masked arrays.
    med = np.median(arr)
    return np.median(np.abs(arr - med))


path=('/Users/danka/Documents/MY_PUBLICATIONS/PremodialBHMicrolensing/'+args.quasar_name+'/map_'+str(args.input))

########################################################################################################################
# Converts microlensing maps from gerlump into magnitude units and writes it down
########################################################################################################################
try:
    print path
    with open(path+'/mapmag.bin') as f:
        mapp = np.fromfile(f,'float')
except:
    print('MAP NOT found :(')
########################################################################################################################

########################################################################################################################
# Reshapes map into 2D and convolves it with source size, if necessary
########################################################################################################################

map2d = mapp.reshape([10000,10000])
if args.conv:
    try:
        with open(path+'/mapconv4.bin') as f:
            map2d_4 = np.fromfile(f,'float')
    except:
            map2d_4 = ndimage.gaussian_filter(map2d, sigma=float(4.)/1.18, order=0)
    try:
        with open(path+'/mapconv7.bin') as f:
            map2d_7 = np.fromfile(f,'float')
    except:
            map2d_7 = ndimage.gaussian_filter(map2d, sigma=float(7.)/1.18, order=0)
    try:
        with open(path+'/mapconv12.bin') as f:
            map2d_12 = np.fromfile(f,'float')
    except:
            map2d_12 = ndimage.gaussian_filter(map2d, sigma=float(12.)/1.18, order=0)
    map2d_4 = map2d_4.reshape([10000,10000])
    map2d_7 = map2d_7.reshape([10000,10000])
    map2d_12 = map2d_12.reshape([10000,10000])

# Performs statistics on real data
########################################################################################################################
source_list=np.loadtxt(open("/Users/danka/Documents/MY_PUBLICATIONS/PremodialBHMicrolensing/"+args.data_file,"rb"))
########################################################################################################################
# Velocity
########################################################################################################################

num_of_days = (max(source_list[:,0])-min(source_list[:,0]))

if args.data_file == "HE0435_all_astrofix_2016.rdb":
    re_size=(sizes.einstain(0.4546,1.689))/1000. #[km]
    velocity = 837.04 #[km/s] Mosquera&Kochanek 2011
if args.data_file == "RXJ1131_datapoints.txt":
    re_size=(sizes.einstain(0.295,0.654))/1000. #[km]
    velocity = 711.77 #[km/s] Mosquera&Kochanek 2011
    
    
total_dist=num_of_days*86400.*velocity  #[km]
pixel_size_km = (re_size*25.)/10000. #[km]
total_dist_in_pixel = total_dist/pixel_size_km

BC = source_list[:,3]-source_list[:,5]
AC = source_list[:,1]-source_list[:,5]
AB = source_list[:,1]-source_list[:,3]
DB = source_list[:,7]-source_list[:,3]

# Creats random lightcurves
########################################################################################################################
guesses_no = int(1E3)
guesses = np.random.randint(9997-total_dist_in_pixel, size=(guesses_no, 2))+1.
likelihood_methods = [std_dev, mad]
best_fits = [(-1, sys.float_info.max)] * len(likelihood_methods)
fit_values = []
for i in range(len(likelihood_methods)):
    fit_values.append([])

########################################################################################################################
distance=np.array(range(total_dist_in_pixel))
x0=guesses[0,0]
y0=guesses[0,1]
import matplotlib
matplotlib.rc('xtick', labelsize=11)
matplotlib.rc('ytick', labelsize=11)
import matplotlib.pyplot as plt

for iteration, (x0, y0) in enumerate(guesses):
    print iteration
    m_slope=(y0-9998.)/(distance.min()-distance.max())
    m_slope_negative=y0/(distance.min()-distance.max())
    finished = False
    while not finished:
        a=np.random.rand(1)*2.-1.#*m_slope
        ys = (distance+x0) * a + y0
        xs = distance+x0
        if max(ys) <9998. and min(ys)> 0. and max(xs) <9998.:
            finished = True
    model_table=np.zeros([guesses_no+1,len(xs)])
    model_table=np.zeros([guesses_no+1,len(xs)])
    model_table4=np.zeros([guesses_no+1,len(xs)])
    model_table7=np.zeros([guesses_no+1,len(xs)])
    model_table12=np.zeros([guesses_no+1,len(xs)])
    
    for index, (x, y) in enumerate(zip(xs, ys)):
        #print index,x, y,int(round(y)),int(round(x)),map2d.shape,map2d_4.shape,len(xs), len(ys)
        model_table[iteration+1,index]=map2d[int(round(y))][int(round(x))]
        if args.conv:
            model_table4[iteration+1,index]=map2d_4[int(round(y))][int(round(x))]
            model_table7[iteration+1,index]=map2d_7[int(round(y))][int(round(x))]
            model_table12[iteration+1,index]=map2d_12[int(round(y))][int(round(x))]
    for method_i, method in enumerate(likelihood_methods):
        fitness = method(model_table[iteration+1])
        fit_values[method_i].append(fitness)

gs = gridspec.GridSpec(5, 5)
f=plt.figure(figsize=(10,10))
ax1 = plt.subplot(gs[0:4,0:4])
ax2 = plt.subplot(gs[-1,0:4])
for tick in ax2.xaxis.get_major_ticks():
    tick.label.set_fontsize(11)
for tick in ax2.yaxis.get_major_ticks():
    tick.label.set_fontsize(11)
ax1.plot([min(xs), max(xs)], [ys[0], ys[len(ys)-1]], 'ro-',linewidth=3)
ax1.imshow(map2d)
ax1.axis('off')

ax2.plot(xs, model_table[1],'ro', markersize=2)

fileName=str('/Users/danka/Documents/MY_PUBLICATIONS/PremodialBHMicrolensing/'+args.quasar_name)+'/map_'+str(args.input)+'/paper'
ax2.set_xlabel("Days",  fontsize=11)
ax2.set_ylabel("Magnitude",  fontsize=11)
ax2.set_ylim([min(model_table[2]),max(model_table[2])])
plt.subplots_adjust(hspace=None)
plt.savefig(fileName, format="png")
        #plt.show()
    
########################################################################################################################
# Saves lists of statistics performed with simulated lightcurves
########################################################################################################################
if args.conv:
    
    #with open(str(quasar_name)+'/'+str(quasar_name)+str(input)+'_SHORTconvolved_'+str(convolution_value)+'model_table.bin', 'w') as outfile:
    #   outfile.write(model_table)
    np.save(path+'model_table4.bin',model_table4 )
    np.save(path+'model_table7.bin',model_table7 )
    np.save(path+'model_table12.bin',model_table12 )
else:
    #  with open(str(quasar_name)+'/'+str(quasar_name)+'_'+str(input)+'_SHORT_NOTconvolved_model_table.bin', 'w') as outfile:
    #outfile.write(model_table)
    np.save(path+'model_table.bin',model_table )

print 'INFO:',guesses_no,'INFO'
########################################################################################################################
########################################################################################################################
if args.conv:
    f = open(path+'convolution_value.txt', 'w')
else:
    f = open(path+'NOTconvolved.txt', 'w')

print '\nStatistics:'
for method, values in zip(likelihood_methods, fit_values):
    values = np.asarray(values)
    f.write(method.func_name)
    temp=('Mean =', np.mean(values), '\n')
    f.write(str(temp))
    temp=('Median =', np.median(values), '\n')
    f.write(str(temp))
    temp=('Mean best 100 =', np.mean(np.sort(values)[:100]), '\n')
    f.write(str(temp))
    if method.func_name == 'std_dev':
        les_0_05 = len(values[values<std_dev(BC)*2.])
        print np.mean(values)
        temp=('Less than 0.05 :', les_0_05, 'or', (100.*les_0_05/len(values)),'%')
        print np.mean(values)
        f.write(str(temp))
    if method.func_name == 'mad':
        les_0_032 = len(values[values<mad(BC)*2.])
        temp=('Less than 0.032 :', les_0_032, 'or', (100.* les_0_032/len(values)),'%')
        f.write(str(temp))
        print np.mean(values)

print '0.'+str(input)+'%'
print 'Less than :', std_dev(BC)*2., les_0_032, 'or', (100.* les_0_032/len(values)),'%',
print 'Less than :', mad(BC)*2., les_0_05, 'or', (100.*les_0_05/len(values)),'%'
print 'Less than :', mad(BC)*2., les_0_05, 'or', (100.*les_0_05/len(values)),'%'

f.close()


