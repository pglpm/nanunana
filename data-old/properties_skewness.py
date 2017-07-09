import os
import sys
import numpy
import pylab
from scipy import stats

pylab.rc(('xtick','ytick','axes'), labelsize=20.0)
pylab.rcParams['font.size']=10.
pylab.rcParams['ytick.minor.pad'] = 10
pylab.rcParams['lines.markersize']=9
pylab.rcParams['lines.linewidth']=2
capsize=5
try:
   ending=sys.argv[1]
   filex=sys.argv[2]
except:
   ending=''
   filex=sys.argv[1]
print 'ending', ending
print 'file', filex
path = '/home/bachmann/data/data_Claudia/data/AD_fMRI/'#data_Claudia/data/AD_fMRI/'

directory = '/RGS_graphprop%s/'%filex

if not os.path.isdir('/home/bachmann/AD_fMRI/results/%s/' %directory):
    os.mkdir('/home/bachmann/AD_fMRI/results/%s/' %directory)

prop = 'weights%s.npy'%ending
prop2 = 'shortest_path%s.npy'%ending#'weight_distribution.npy'
prop3 = 'clustering_coefficient%s.npy'%ending
prop4= 'degree%s.npy'%ending
prop5 = 'modul%s.npy'%ending


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
 


group_names = ['Control','MCI','Alzheimer']
colors=['darkslategray','forestgreen','darkkhaki']#['0.25','0.5','0.75']
fs=35

## weight distribution
pylab.figure('summarize',figsize=(20,8))
pylab.clf()
for groupid, group in enumerate(groups):
    graphweights_all=[]
    sp_all=[]
    cc_all=[]
    wd_all=[]
    nwd_all=[]
    ns_all=[]
    wcc_all=[]
    mod_all=[]
    mean_all = []
    std_all = []
    nz_all = []
    for ii in group:
        
        if os.path.isfile(path+ii+directory+prop2):
             ## graph weights
             graphweights = numpy.load(path+ii+directory+prop)
             graphweights=1./graphweights.flatten()
             graphweights_all += [stats.moment(graphweights, moment=3)]
             print ii, stats.moment(graphweights, moment=3)
             ## shortest path
             sp = numpy.load(path+ii+directory+prop2)
             sp=sp.flatten()
             sp_all += [stats.moment(sp, moment=3)]
             ## cluster coefficient
             cc = numpy.load(path+ii+directory+prop3)
             cc = cc.flatten()
             cc_all += [stats.moment(cc, moment=3)]
             ## weighted degree
             wd = numpy.load(path+ii+directory+prop4)
             wd=wd.flatten()
             wd_all += [stats.moment(wd, moment=3)]
             ## modularity
             mod = numpy.load(path+ii+directory+prop5)
             mod_all += [mod]  
             ## weighted closeness centrality
             sp = numpy.load(path+ii+directory+prop2)
             wcc=sp.flatten()
             gsize=int(numpy.sqrt(len(wcc)))
             wcc=numpy.reshape(wcc,(gsize,gsize))
             gw_avg=[]
             for k in range(gsize):
                
                gw_avg.append(gsize/numpy.nansum(wcc[k]))
             wcc_all += [stats.moment(gw_avg, moment=3)]   
             ## normalized weighted degree
             nwd = numpy.load(path+ii+directory+prop4)
             nwd=nwd.flatten()
             w=1./numpy.load(path+ii+directory+prop)
             s=len(nwd)
             w.sort()
             w=w[-s:-1]
             if sum(w)!=0.:
                 nwd= nwd*(1./(sum(w)))*numpy.ones(s)
             else:
                 nwd=numpy.zeros(len(nwd)) 
             nwd_all += [stats.moment(nwd, moment=3)]   
             ## number of nodes
             ns_all += [len(wd)]
        else:
	     print 'data %s does not exist' %path+ii+directory+prop4
    ro=10
    print group_names[groupid]
    print ' graphweights:', ', '.join(map(str, list(numpy.round(graphweights_all,ro).astype(str)))) 
    print ' shortest path:',', '.join(map(str, list(numpy.round(sp_all,ro).astype(str)))) 
    print ' cluster coeff.:',', '.join(map(str, list(numpy.round(cc_all,ro).astype(str)))) 
    print ' weighted degree:',', '.join(map(str, list(numpy.round(wd_all,ro).astype(str)))) 
    print ' normalized weighted degree',', '.join(map(str, list(numpy.round(nwd_all,ro).astype(str)))) 
    print ' node size:',', '.join(map(str, list(numpy.round(ns_all,ro).astype(str)))) 
    print ' closeness centrality:',', '.join(map(str, list(numpy.round(wcc_all,ro).astype(str)))) 
    print ' modularity:',', '.join(map(str, list(numpy.round(mod_all,ro).astype(str)))) 
            
   
    pylab.figure('summarize')
    pylab.subplot(1,8,1)
    pylab.plot(numpy.ones(len(ns_all))*groupid,ns_all,'.',color=colors[groupid])	
    
    pylab.errorbar(groupid,numpy.mean(ns_all), yerr=numpy.std(ns_all), fmt='k_',capsize=capsize)
    pylab.xlim((-1,3))  
    pylab.ylabel('number of nodes') 
    pylab.xticks([0.,1.,2.], ['C','MCI','AD'])

    pylab.subplot(1,8,2)
    graphweights_all=numpy.array(graphweights_all).astype(float)
    graphweights_all[numpy.where(numpy.isnan(graphweights_all)==True)]=0.
    pylab.plot(numpy.ones(len(graphweights_all))*groupid,graphweights_all,'.',color=colors[groupid])
#for en in ends:
#    for fi in filx:
#        os.system('ipython properties_mean.py ' + en+ ' '+str(fi))
	
    
    pylab.errorbar(groupid,numpy.mean(graphweights_all), yerr=numpy.std(graphweights_all), fmt='k_',capsize=capsize)
    pylab.xlim((-1,3))   
    pylab.ylabel('3rd moment averaged edge weights') 
    pylab.xticks([0.,1.,2.], ['C','MCI','AD'])
    
    pylab.subplot(1,8,3)
    
    pylab.plot(numpy.ones(len(sp_all))*groupid,sp_all,'.',color=colors[groupid])
    
    pylab.errorbar(groupid,numpy.mean(sp_all), yerr=numpy.std(sp_all), fmt='k_',capsize=capsize)
    pylab.xlim((-1,3))   
    pylab.ylabel('3rd moment averaged shortest path') 
    pylab.xticks([0.,1.,2.], ['C','MCI','AD'])

    pylab.subplot(1,8,4)
    cc_all=numpy.array(cc_all).astype(float)
    cc_all[numpy.where(numpy.isnan(cc_all)==True)]=0.

    pylab.plot(numpy.ones(len(cc_all))*groupid,cc_all,'.',color=colors[groupid])
    
    pylab.errorbar(groupid,numpy.mean(cc_all), yerr=numpy.std(cc_all), fmt='k_',capsize=capsize)	
    pylab.xlim((-1,3))
    pylab.ylabel('3rd moment averaged cluster coefficient') 
    pylab.xticks([0.,1.,2.], ['C','MCI','AD'])

    pylab.subplot(1,8,5)
    wd_all=numpy.array(wd_all).astype(float)
    wd_all[numpy.where(numpy.isnan(wd_all)==True)]=0.
    pylab.plot(numpy.ones(len(wd_all))*groupid,wd_all,'.',color=colors[groupid])
    
    pylab.errorbar(groupid,numpy.mean(wd_all), yerr=numpy.std(wd_all), fmt='k_',capsize=capsize)	
    pylab.xlim((-1,3))
    pylab.ylabel('3rd moment averaged weighted degree')      
    pylab.xticks([0.,1.,2.], ['C','MCI','AD'])

    pylab.subplot(1,8,6)
    pylab.plot(numpy.ones(len(mod_all))*groupid,mod_all,'.',color=colors[groupid])
    
    pylab.errorbar(groupid,numpy.mean(mod_all), yerr=numpy.std(mod_all), fmt='k_',capsize=capsize)	
    pylab.xlim((-1,3))
    pylab.ylabel('modularity')      
    pylab.xticks([0.,1.,2.], ['C','MCI','AD'])
 
    pylab.subplot(1,8,7)
    pylab.plot(numpy.ones(len(wcc_all))*groupid,wcc_all,'.',color=colors[groupid])
    
    pylab.errorbar(groupid,numpy.mean(wcc_all), yerr=numpy.std(wcc_all), fmt='k_',capsize=capsize)	
    pylab.xlim((-1,3))
    pylab.ylabel('3rd moment weighted closeness centrality')      
    pylab.xticks([0.,1.,2.], ['C','MCI','AD'])

    pylab.subplot(1,8,8)
    pylab.plot(numpy.ones(len(nwd_all))*groupid,nwd_all,'.',color=colors[groupid])
    
    pylab.errorbar(groupid,numpy.mean(nwd_all), yerr=numpy.std(nwd_all), fmt='k_',capsize=capsize)	
    pylab.xlim((-1,3))
    pylab.ylabel('3rd moment normalized weighted degree')      
    pylab.xticks([0.,1.,2.], ['C','MCI','AD'])


    
    
    pylab.figure('compare gw-sp')
    pylab.plot(graphweights_all,sp_all,'.', label=group_names[groupid],color=colors[groupid])
    pylab.xlabel('3rd moment edge weights') 
    pylab.ylabel('3rd moment shortest path') 
    pylab.savefig('/home/bachmann/AD_fMRI/results/%s/compare_gw_sp%s_skew.eps' %(directory,ending) )
    

    pylab.figure('compare gw-cc')
    pylab.plot(graphweights_all,cc_all,'.', label=group_names[groupid],color=colors[groupid])
    pylab.xlabel('3rd moment edge weights')
    pylab.ylabel('3rd moment cluster coefficient') 
    pylab.savefig('/home/bachmann/AD_fMRI/results/%s/compare_gw_cc%s_skew.eps' %(directory,ending) )

    pylab.figure('compare gw-wd')
    pylab.plot(graphweights_all,wd_all,'.', label=group_names[groupid],color=colors[groupid])
    pylab.xlabel('3rd moment edge weights')
    pylab.ylabel('3rd moment weighted degree')
    pylab.savefig('/home/bachmann/AD_fMRI/results/%s/compare_gw_wd%s_skew.eps' %(directory,ending) )
    
    pylab.figure('compare gw-mod')
    pylab.plot(graphweights_all,mod_all,'.', label=group_names[groupid],color=colors[groupid])
    pylab.xlabel('3rd moment edge weights')
    pylab.ylabel('modularity')
    pylab.savefig('/home/bachmann/AD_fMRI/results/%s/compare_gw_mod_skew%s.eps' %(directory,ending) )

    pylab.figure('compare gw-ns')
    pylab.plot(graphweights_all,ns_all,'.', label=group_names[groupid],color=colors[groupid])
    pylab.xlabel('3rd moment edge weights')
    pylab.ylabel('number of nodes') 
    pylab.savefig('/home/bachmann/AD_fMRI/results/%s/compare_gw_ns_skew%s.eps' %(directory,ending) )

    pylab.figure('compare gw-wcc')
    pylab.plot(graphweights_all,wcc_all,'.', label=group_names[groupid],color=colors[groupid])
    pylab.xlabel('3rd moment edge weights')
    pylab.ylabel('3rd moment weighted closeness centrality')    
    pylab.savefig('/home/bachmann/AD_fMRI/results/%s/compare_gw_wcc_skew%s.eps' %(directory,ending) )

    pylab.figure('compare gw-nwd')
    pylab.plot(graphweights_all,nwd_all,'.', label=group_names[groupid],color=colors[groupid])
    pylab.xlabel('3rd moment edge weights')
    pylab.ylabel('3th moment normalized weighted degree')
    pylab.savefig('/home/bachmann/AD_fMRI/results/%s/compare_gw_nwd_skew%s.eps' %(directory,ending) )
    '''
    pylab.figure('compare gw-mean')
    pylab.plot(graphweights_all,mean_all,'.', label=group_names[groupid])
    pylab.figure('compare gw-std')
    pylab.plot(graphweights_all,std_all,'.', label=group_names[groupid])
    print 'yes'
    '''

pylab.figure('summarize')    
pylab.tight_layout()
#pylab.show()
if not os.path.isdir('/home/bachmann/AD_fMRI/results/%s/' %directory):
    os.mkdir('/home/bachmann/AD_fMRI/results/%s/' %directory)
pylab.savefig('/home/bachmann/AD_fMRI/results/%s/result_overview%s_skew.eps' %(directory,ending) )
filx=[22,16]
ends=['','0_1','0_2','0_3','0_4','0_5','0_6','0_7','0_8','_richclub']
#for en in ends:
#    for fi in filx:
#        os.system('ipython properties_mean.py ' + en+ ' '+str(fi))

