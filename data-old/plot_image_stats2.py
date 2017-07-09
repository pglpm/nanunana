import numpy
import pylab
from mpl_toolkits.mplot3d import Axes3D
import os
import itertools
pylab.ion()
pylab.rc(('xtick','ytick','axes'), labelsize=20.0)
pylab.rcParams['font.size']=10.
pylab.rcParams['ytick.minor.pad'] = 10
pylab.rcParams['lines.markersize']=9
pylab.rcParams['lines.linewidth']=2
colors=['forestgreen','magenta','darkkhaki']
groups=[['C001_p0.001',
'C002_p0.001',
'C003_p0.001',
'C005_p0.001', 
'C006_p0.001', 
'C007_p0.001',
'C008_p0.001',
'C009_p0.001',
'C011_p0.001', 
'C012_p0.001',
'C016_p0.001', 
'C018_p0.001', 
'C021_p0.001', 
'C023_p0.001', 
'C024_p0.001', 
'C025_p0.001', 
'C026_p0.001', 
'C027_p0.001', 
'C028_p0.001', 
'C030_p0.001', 
'C031_p0.001', 
'C032_p0.001', 
'C033_p0.001', 
'C903_p0.001',
'C901_p0.001',
'C902_p0.001'],
['M002_p0.001',
'M008_p0.001',
'M010_p0.001',
'M015_p0.001',
'M016_p0.001',
'M019_p0.001', 
'M020_p0.001', 
'M021_p0.001', 
'M024_p0.001', 
'M029_p0.001', 
'M030_p0.001', 
'M032_p0.001',
'M038_p0.001',
'M040_p0.001', 
'M901_p0.001', 
'M902_p0.001'],
['A003_p0.001',
'A004_p0.001',
'A005_p0.001',
'A007_p0.001',
'A008_p0.001',
'A901_p0.001', 
'A902_p0.001', 
'A903_p0.001', 
'A905_p0.001', 
'A906_p0.001', 
'A907_p0.001', 
'A908_p0.001', 
'A909_p0.001', 
'A910_p0.001']]

#'M002_p0.001', 'M008_p0.001', 'M010_p0.001', 'M015_p0.001', 'M016_p0.001', 'M019_p0.001', 
#         'M020_p0.001', 'M021_p0.001', 'M024_p0.001', 'M029_p0.001', 'M030_p0.001', 'M032_p0.001',
#         'M040_p0.001'
group_names = ['Control','Alzheimer','MCI']
path = '/home/bachmann/data/data_Claudia/data/AD_fMRI/'
path2 = '/home/bachmann/AD_fMRI/results/'
dir2= '/stats/'



for groupid, group in enumerate(groups):
    mean1_all = []
    std1_all = []
    mean2_all = []
    std2_all = []
    mean3_all = []
    std3_all = []
    lenzn_all = []
    for mm,ii in enumerate(group):
        if os.path.isfile(path+ii+dir2+'mean1.npy'):
            len_nz = numpy.load(path+ii+'/stats/'+'len_nz.npy')
            sigma1 = numpy.load(path+ii+'/stats/'+'sig1.npy') # across time
            mean1 = numpy.load(path+ii+'/stats/'+'mean1.npy')
            sigma2 = numpy.load(path+ii+'/stats/'+'sig2.npy') # across voxels
            mean2 = numpy.load(path+ii+'/stats/'+'mean2.npy')
            sigma3 = numpy.load(path+ii+'/stats/'+'sig3.npy')
            mean3 = numpy.load(path+ii+'/stats/'+'mean3.npy')
            #print 'session: ', ii,  numpy.mean(mean1)
        if type(mean1)!=list:
            lenzn_all += [float(len_nz)]
            mean1_all += [float(numpy.mean(mean1))]
	    std1_all += [float(numpy.mean(sigma1))]
            mean2_all += [float(numpy.mean(mean2))]
	    std2_all += [float(numpy.mean(sigma2))]
            mean3_all += [float(mean3)]
	    std3_all += [float(sigma3)]
        else:
	    print 'file not found: ',path+ii+dir2+'mean1.npy'
    print 'groupid', groupid
    print 'mean1', mean1_all
    print 'std1', std1_all
    pylab.figure('mean1-sigma1')
    pylab.plot(mean1_all,std1_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('mean2-sigma2')
    pylab.plot(mean2_all,std2_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('mean3-sigma3')
    pylab.plot(mean3_all,std3_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('lenzn-sigma1')
    pylab.plot(lenzn_all,std1_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('lenzn-sigma2')
    pylab.plot(lenzn_all,std2_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('lenzn-sigma3')
    pylab.plot(lenzn_all,std3_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('lenzn-mean1')
    pylab.plot(lenzn_all,mean1_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('lenzn-mean2')
    pylab.plot(lenzn_all,mean2_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('lenzn-mean3')
    pylab.plot(lenzn_all,mean3_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('mean1-mean2')
    pylab.plot(mean1_all,mean2_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('mean1-mean3')
    pylab.plot(mean1_all,mean3_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('sigma1-sigma2')
    pylab.plot(std1_all,std2_all,'.',color=colors[groupid], label=group_names[groupid])
    pylab.figure('sigma1-sigma3')
    pylab.plot(std1_all,std3_all,'.',color=colors[groupid], label=group_names[groupid])
    if groupid==0:
        fig = pylab.figure('mean1-std1-lenzn')
        ax = fig.add_subplot(111, projection='3d')
    
    ax.scatter(mean1_all, std1_all, lenzn_all, color=colors[groupid], s=30)
    ax.set_xlabel('mean1')
    ax.set_ylabel('std1')
    ax.set_zlabel('non zero nodes')

pylab.figure('mean1-sigma1')
pylab.xlabel('mean across time')
pylab.ylabel('std across time')
pylab.legend()
pylab.savefig(path2+dir2+'mean1_std1.eps')

pylab.figure('mean2-sigma2')
pylab.xlabel('mean across voxels')
pylab.ylabel('std across voxels')
pylab.legend()
pylab.savefig(path2+dir2+'mean2_std2.eps')

pylab.figure('mean3-sigma3')
pylab.xlabel('mean across time and voxels')
pylab.ylabel('std across time and voxels')
pylab.legend()
pylab.savefig(path2+dir2+'mean3_std3.eps')

pylab.figure('lenzn-sigma1')
pylab.xlabel('number of nodes')
pylab.ylabel('std across time')
pylab.legend()
pylab.savefig(path2+dir2+'lenzn_std1.eps')

pylab.figure('lenzn-sigma2')
pylab.xlabel('number of nodes')
pylab.ylabel('std across voxels')
pylab.legend()
pylab.savefig(path2+dir2+'lenzn_std2.eps')

pylab.figure('lenzn-sigma3')
pylab.xlabel('number of nodes')
pylab.ylabel('std across time and voxels')
pylab.legend()
pylab.savefig(path2+dir2+'lenzn_std3.eps')

pylab.figure('lenzn-mean1')
pylab.xlabel('number of nodes')
pylab.ylabel('mean across time')
pylab.legend()
pylab.savefig(path2+dir2+'lenzn_mean1.eps')

pylab.figure('lenzn-mean2')
pylab.xlabel('number of nodes')
pylab.ylabel('mean across voxels')
pylab.legend()
pylab.savefig(path2+dir2+'lenzn_mean2.eps')

pylab.figure('lenzn-mean3')
pylab.xlabel('number of nodes')
pylab.ylabel('mean across time and voxels')
pylab.legend()
pylab.savefig(path2+dir2+'lenzn_mean3.eps')

pylab.figure('mean1-mean2')
pylab.xlabel('mean across time')
pylab.ylabel('mean across voxels')
pylab.legend()
pylab.savefig(path2+dir2+'mean1_mean2.eps')

pylab.figure('mean1-mean3')
pylab.xlabel('mean across time')
pylab.ylabel('mean across time and voxels')
pylab.legend()
pylab.savefig(path2+dir2+'mean1_mean3.eps')

pylab.figure('sigma1-sigma2')
pylab.xlabel('std across time')
pylab.ylabel('std across voxels')
pylab.legend()
pylab.savefig(path2+dir2+'sigma1_sigma2.eps')

pylab.figure('sigma1-sigma3')
pylab.xlabel('std across time')
pylab.ylabel('std across time and voxels')
pylab.legend()
pylab.savefig(path2+dir2+'sigma1_sigma3.eps')


## std of variance2
'''
pylab.figure(figsize=(20,10))
num_bins=20
ran = (0,40000.)
for groupid, group in enumerate(groups):
    graphweights_all=[]
    x = []
    for ii in group:
        #print 'hi', path+ii+directory+prop   
        if os.path.isfile(path+ii+'/stats/'+'sig2.npy'):
             graphweights = numpy.load(path+ii+'/stats/'+'sig2.npy')
             norm=1./len(graphweights)*numpy.ones(len(graphweights))     
	     
             graphweights_all += [graphweights]
	     x += [numpy.histogram(graphweights,range=ran,bins=num_bins, weights=norm)[0]]
             print ii, len(graphweights)
    ga=numpy.concatenate(graphweights_all)
    norm1=1./len(ga)*numpy.ones(len(ga))
    hist_all=numpy.histogram(ga,range=ran,bins=num_bins, weights=norm1)[0]
    yer=numpy.std(x,axis=0)
    xvalues = numpy.arange(len(hist_all)*(1+len(groups)))[::3]+groupid
    pylab.bar(xvalues, hist_all, color=colors[groupid], yerr=yer, label=group_names[groupid], ecolor='royalblue',error_kw={'elinewidth':3})
pylab.ylim((0,0.6))
pylab.xticks(xvalues[::2],numpy.arange(0,ran[1],ran[1]/num_bins)[::2],fontsize=20)

#pylab.xlim((3,xvalues[-1]+2))
pylab.legend(fontsize=20)
pylab.xlabel('std across voxels',fontsize=20)
pylab.ylabel('relative frequency',fontsize=20)
#pylab.savefig(path2+'edge_weight_pop.eps')
#pylab.savefig(path2+'edge_weight_pop.pdf')
#pylab.savefig(path2+'edge_weight_pop.jpg')
pylab.savefig(path2+dir2+'dist_std2.eps')
pylab.figure(figsize=(20,10))
num_bins=20
ran = (1.2e7,2.6e7)
for groupid, group in enumerate(groups):
    graphweights_all=[]
    x = []
    for ii in group:
        #print 'hi', path+ii+directory+prop   
        if os.path.isfile(path+ii+'/stats/'+'sig1.npy'):
             graphweights = numpy.load(path+ii+'/stats/'+'sig1.npy')
             norm=1./len(graphweights)*numpy.ones(len(graphweights))     
	     
             graphweights_all += [graphweights]
	     x += [numpy.histogram(graphweights,range=ran,bins=num_bins, weights=norm)[0]]
             print ii, len(graphweights)
    ga=numpy.concatenate(graphweights_all)
    norm1=1./len(ga)*numpy.ones(len(ga))
    hist_all=numpy.histogram(ga,range=ran,bins=num_bins, weights=norm1)[0]
    yer=numpy.std(x,axis=0)
    xvalues = numpy.arange(len(hist_all)*(1+len(groups)))[::3]+groupid
    pylab.bar(xvalues, hist_all, color=colors[groupid], yerr=yer, label=group_names[groupid], ecolor='royalblue',error_kw={'elinewidth':3})
pylab.ylim((0,0.4))
pylab.xticks(xvalues[::2],numpy.arange(ran[0],ran[1],(ran[1]-ran[0])/num_bins)[::2],fontsize=20)

#pylab.xlim((3,xvalues[-1]+2))
pylab.legend(fontsize=20)
pylab.xlabel('std across time',fontsize=20)
pylab.ylabel('relative frequency',fontsize=20)
pylab.savefig(path2+dir2+'dist_std1.eps')
pylab.figure(figsize=(20,10))
for groupid, group in enumerate(groups):
    graphweights_all=[]
    x = []
    for ii in group:
        #print 'hi', path+ii+directory+prop   
        if os.path.isfile(path+ii+'/stats/'+'sig1.npy'):
             graphweights = numpy.load(path+ii+'/stats/'+'sig1.npy')
             graphweights_all += [graphweights]
	     
    
    norm1=1./len(ga)*numpy.ones(len(ga))
    hist_all=numpy.mean(graphweights_all, axis=0)[::3]
    yer=numpy.std(graphweights_all, axis=0)[::3]
    xvalues = numpy.arange(len(hist_all)*(1+len(groups)))[::3]+groupid
    pylab.bar(xvalues, hist_all, color=colors[groupid], yerr=yer, label=group_names[groupid], ecolor='royalblue',error_kw={'elinewidth':3})
pylab.xlim((0,132))
#pylab.xticks(xvalues[::2],numpy.arange(ran[0],ran[1],(ran[1]-ran[0])/num_bins)[::2],fontsize=20)

#pylab.xlim((3,xvalues[-1]+2))
pylab.legend(fontsize=20)
pylab.ylabel('std across time',fontsize=20)
pylab.xlabel('time point',fontsize=20)
pylab.savefig(path2+dir2+'std1_time.eps')
### mean

pylab.figure(figsize=(20,10))
num_bins=20
ran = (0,40000.)
for groupid, group in enumerate(groups):
    graphweights_all=[]
    x = []
    for ii in group:
        #print 'hi', path+ii+directory+prop   
        if os.path.isfile(path+ii+'/stats/'+'mean2.npy'):
             graphweights = numpy.load(path+ii+'/stats/'+'mean2.npy')
             norm=1./len(graphweights)*numpy.ones(len(graphweights))     
	     
             graphweights_all += [graphweights]
	     x += [numpy.histogram(graphweights,range=ran,bins=num_bins, weights=norm)[0]]
             print ii, len(graphweights)
    ga=numpy.concatenate(graphweights_all)
    norm1=1./len(ga)*numpy.ones(len(ga))
    hist_all=numpy.histogram(ga,range=ran,bins=num_bins, weights=norm1)[0]
    yer=numpy.std(x,axis=0)
    xvalues = numpy.arange(len(hist_all)*(1+len(groups)))[::3]+groupid
    pylab.bar(xvalues, hist_all, color=colors[groupid], yerr=yer, label=group_names[groupid], ecolor='royalblue',error_kw={'elinewidth':3})
pylab.ylim((0,0.6))
pylab.xticks(xvalues[::2],numpy.arange(0,ran[1],ran[1]/num_bins)[::2],fontsize=20)

#pylab.xlim((3,xvalues[-1]+2))
pylab.legend(fontsize=20)
pylab.xlabel('mean across voxels',fontsize=20)
pylab.ylabel('relative frequency',fontsize=20)
#pylab.savefig(path2+'edge_weight_pop.eps')
#pylab.savefig(path2+'edge_weight_pop.pdf')
#pylab.savefig(path2+'edge_weight_pop.jpg')
pylab.savefig(path2+dir2+'dist_mean2.eps')

pylab.figure(figsize=(20,10))
num_bins=20
ran = (0.85e4,1.1e4)
for groupid, group in enumerate(groups):
    graphweights_all=[]
    x = []
    for ii in group:
        #print 'hi', path+ii+directory+prop   
        if os.path.isfile(path+ii+'/stats/'+'mean1.npy'):
             graphweights = numpy.load(path+ii+'/stats/'+'mean1.npy')
             norm=1./len(graphweights)*numpy.ones(len(graphweights))     
	     
             graphweights_all += [graphweights]
	     x += [numpy.histogram(graphweights,range=ran,bins=num_bins, weights=norm)[0]]
             print ii, len(graphweights)
    ga=numpy.concatenate(graphweights_all)
    norm1=1./len(ga)*numpy.ones(len(ga))
    hist_all=numpy.histogram(ga,range=ran,bins=num_bins, weights=norm1)[0]
    yer=numpy.std(x,axis=0)
    xvalues = numpy.arange(len(hist_all)*(1+len(groups)))[::3]+groupid
    pylab.bar(xvalues, hist_all, color=colors[groupid], yerr=yer, label=group_names[groupid], ecolor='royalblue',error_kw={'elinewidth':3})
pylab.ylim((0,0.4))
pylab.xticks(xvalues[::2],numpy.arange(ran[0],ran[1],(ran[1]-ran[0])/num_bins)[::2],fontsize=20)
pylab.savefig(path2+dir2+'dist_mean1.eps')
#pylab.xlim((3,xvalues[-1]+2))
pylab.legend(fontsize=20)
pylab.xlabel('mean across time',fontsize=20)
pylab.ylabel('relative frequency',fontsize=20)

pylab.figure(figsize=(20,10))
for groupid, group in enumerate(groups):
    graphweights_all=[]
    x = []
    for ii in group:
        #print 'hi', path+ii+directory+prop   
        if os.path.isfile(path+ii+'/stats/'+'mean1.npy'):
             graphweights = numpy.load(path+ii+'/stats/'+'mean1.npy')
             graphweights_all += [graphweights]
	     
    
    norm1=1./len(ga)*numpy.ones(len(ga))
    hist_all=numpy.mean(graphweights_all, axis=0)[::3]
    yer=numpy.std(graphweights_all, axis=0)[::3]
    xvalues = numpy.arange(len(hist_all)*(1+len(groups)))[::3]+groupid
    pylab.bar(xvalues, hist_all, color=colors[groupid], yerr=yer, label=group_names[groupid], ecolor='royalblue',error_kw={'elinewidth':3})
pylab.xlim((0,132))
#pylab.xticks(xvalues[::2],numpy.arange(ran[0],ran[1],(ran[1]-ran[0])/num_bins)[::2],fontsize=20)

#pylab.xlim((3,xvalues[-1]+2))
pylab.legend(fontsize=20)
pylab.ylabel('mean across time',fontsize=20)
pylab.xlabel('time point',fontsize=20)
pylab.savefig(path2+dir2+'mean1_time.eps')
'''

