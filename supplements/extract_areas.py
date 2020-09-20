import sys
import numpy
import pylab
import nibabel as nib
import h5py_wrapper as h5

data_path=sys.argv[1]
'''
list_areas=[('Precentral_L', 2001),
('Precentral_R', 2002),
('Frontal_Sup_L', 2101),
('Frontal_Sup_R', 2102),
('Frontal_Sup_Orb_L', 2111),
('Frontal_Sup_Orb_R', 2112),
('Frontal_Mid_L', 2201),
('Frontal_Mid_R', 2202),
('Frontal_Mid_Orb_L', 2211),
('Frontal_Mid_Orb_R', 2212),
('Frontal_Inf_Oper_L', 2301),
('Frontal_Inf_Oper_R' ,2302),
('Frontal_Inf_Tri_L' ,2311),
('Frontal_Inf_Tri_R' ,2312),
('Frontal_Inf_Orb_L', 2321),
('Frontal_Inf_Orb_R', 2322),
('Rolandic_Oper_L', 2331),
('Rolandic_Oper_R', 2332),
('Supp_Motor_Area_L', 2401),
('Supp_Motor_Area_R', 2402),
('Olfactory_L', 2501),
('Olfactory_R', 2502),
('Frontal_Sup_Medial_L', 2601),
('Frontal_Sup_Medial_R',2602),
('Frontal_Med_Orb_L', 2611),
('Frontal_Med_Orb_R', 2612),
('Rectus_L', 2701),
('Rectus_R', 2702),
('Insula_L', 3001),
('Insula_R', 3002),
('Cingulum_Ant_L', 4001),
('Cingulum_Ant_R', 4002),
('Cingulum_Mid_L', 4011),
('Cingulum_Mid_R', 4012),
('Cingulum_Post_L', 4021),
('Cingulum_Post_R', 4022),
('Hippocampus_L', 4101),
('Hippocampus_R', 4102),
('ParaHippocampal_L', 4111),
('ParaHippocampal_R', 4112),
('Amygdala_L', 4201),
('Amygdala_R', 4202),
('Calcarine_L', 5001),
('Calcarine_R', 5002),
('Cuneus_L', 5011),
('Cuneus_R', 5012),
('Lingual_L', 5021),
('Lingual_R', 5022),
('Occipital_Sup_L' ,5101),
('Occipital_Sup_R', 5102),
('Occipital_Mid_L', 5201),
('Occipital_Mid_R', 5202),
('Occipital_Inf_L', 5301),
('Occipital_Inf_R', 5302),
('Fusiform_L', 5401),
('Fusiform_R', 5402),
('Postcentral_L', 6001),
('Postcentral_R', 6002),
('Parietal_Sup_L', 6101),
('Parietal_Sup_R', 6102),
('Parietal_Inf_L', 6201),
('Parietal_Inf_R', 6202),
('SupraMarginal_L', 6211),
('SupraMarginal_R', 6212),
('Angular_L', 6221),
('Angular_R', 6222),
('Precuneus_L', 6301),
('Precuneus_R', 6302),
('Paracentral_Lobule_L', 6401),
('Paracentral_Lobule_R', 6402),
('Caudate_L', 7001),
('Caudate_R', 7002),
('Putamen_L', 7011),
('Putamen_R', 7012),
('Pallidum_L', 7021),
('Pallidum_R', 7022),
('Thalamus_L', 7101),
('Thalamus_R', 7102),
('Heschl_L', 8101),
('Heschl_R', 8102),
('Temporal_Sup_L', 8111),
('Temporal_Sup_R', 8112),
('Temporal_Pole_Sup_L', 8121),
('Temporal_Pole_Sup_R', 8122),
('Temporal_Mid_L', 8201),
('Temporal_Mid_R', 8202),
('Temporal_Pole_Mid_L', 8211),
('Temporal_Pole_Mid_R', 8212),
('Temporal_Inf_L', 8301),
('Temporal_Inf_R', 8302),
('Cerebelum_Crus1_L', 9001),
('Cerebelum_Crus1_R', 9002),
('Cerebelum_Crus2_L', 9011),
('Cerebelum_Crus2_R', 9012),
('Cerebelum_3_L', 9021),
('Cerebelum_3_R', 9022),
('Cerebelum_4_5_L', 9031),
('Cerebelum_4_5_R', 9032),
('Cerebelum_6_L', 9041),
('Cerebelum_6_R', 9042),
('Cerebelum_7b_L', 9051),
('Cerebelum_7b_R', 9052),
('Cerebelum_8_L', 9061),
('Cerebelum_8_R', 9062),
('Cerebelum_9_L', 9071),
('Cerebelum_9_R', 9072),
('Cerebelum_10_L', 9081),
('Cerebelum_10_R', 9082),
('Vermis_1_2', 9100),
('Vermis_3', 9110),
('Vermis_4_5', 9120),
('Vermis_6', 9130),
('Vermis_7', 9140),
('Vermis_8', 9150),
('Vermis_9', 9160),
('Vermis_10', 9170)]
'''
area_list=[ 1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17,
       18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34,
       35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51,
       52, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65, 66, 67, 68,
       69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85,
       86, 87, 88, 89, 90, 91, 92, 95, 96]
list_areas=[('area%d'%ii,ii) for ii in area_list]
# load MNI_low
#data_path='/home/bachmann/data/modafinil/sub01_pre1/'
MNI_atlas = nib.load('/usr/share/fsl/data/atlases/HarvardOxford/HarvardOxford-cortl-maxprob-thr50-2mm.nii.gz') 
MNI_data=MNI_atlas.get_data()
sp=MNI_data.shape
reg_dic={}
reg_dic2={}
for reg, value in list_areas:
    #print reg
    cord1=numpy.where(MNI_data==value)
    cord2=numpy.array(cord1).transpose()
    reg_dic[reg]={}
    for mm in cord2:
	reg_dic[reg]['%d,%d,%d'%(mm[0],mm[1],mm[2])]=numpy.zeros(141)
for num in range(10):
    data=nib.load(data_path+'/vol000%d_stand.nii.gz'%num).get_data() 
    for key1 in reg_dic.keys():
	#print key1
        for key2 in reg_dic[key1].keys():
            x,y,z=[int(co) for co in key2.split(',')]
            reg_dic[key1][key2][num]=data[x,y,z]
for num in range(10,100):
    data=nib.load(data_path+'/vol00%d_stand.nii.gz'%num).get_data() 
    for key1 in reg_dic.keys():
	#print key1
        for key2 in reg_dic[key1].keys():
            x,y,z=[int(co) for co in key2.split(',')]
            reg_dic[key1][key2][num]=data[x,y,z]
for num in range(100,140):
    data=nib.load(data_path+'/vol0%d_stand.nii.gz'%num).get_data() 
    for key1 in reg_dic.keys():
	print key1
        for key2 in reg_dic[key1].keys():
            x,y,z=[int(co) for co in key2.split(',')]
            reg_dic[key1][key2][num]=data[x,y,z]
for key1 in reg_dic.keys():
    reg_dic2[key1]=numpy.mean(reg_dic[key1].values(), axis=0)
print reg_dic2
h5.add_to_h5(data_path+'/dic_ROI_activity.h5',reg_dic2)


