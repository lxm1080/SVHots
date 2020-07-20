# -*- coding: utf-8 -*-
"""
Created on Tue Dec 24 15:48:18 2019

@author: Lv
"""


def indextable(cancer,element):
    
    abbreviate={
    'ACC':'Adrenocortical_carcinoma',
    'BLCA':'Bladder_urothelial_carcinoma',
    'BRCA':'Breast_invasive_carcinoma',
    'CESC':'Cervical_squamous_cell_carcinoma_and_endocervical_adenocarcinoma',
    'CHOL':'Cholangiocarcinoma',
    'COAD':	'Colon_adenocarcinoma',
    'ESCA'	:'Esophageal_carcinoma',
    'GBM'	:'Glioblastoma_multiforme',
    'HNSC'	:'Head_and_neck_squamous_cell_carcinoma',
    'KICH'	:'Kidney_chromophobe',
    'KIRC'	:'Kidney_renal_clear_cell_carcinoma',
    'KIRP'	:'Kidney_renal_papillary_cell_carcinoma',
    'LAML'	:'Acute_myeloid_leukemia',
    'LGG':    'Brain_lower_grade_glioma',
    'LIHC'	:'Liver_hepatocellular_carcinoma',
    'LUAD'	:'Lung_adenocarcinoma',
    'LUSC'	:'Lung_squamous_cell_carcinoma',
    'OV'	:'Ovarian_serous_cystadenocarcinoma',
    'PAAD'	:'Pancreatic_adenocarcinoma',
    'PCPG'	:'Pheochromocytoma_and_paraganglioma',
    'PRAD'	:'Prostate_adenocarcinoma',
    'READ'	:'Rectum_adenocarcinoma',
    'SARC'	:'Sarcoma',
    'SKCM'	:'Skin_cutaneous_melanoma',
    'STAD'	:'Stomach_adenocarcinoma',
    'TGCT'	:'Testicular_germ_cell_tumor',
    'THCA'	:'Thyroid_carcinoma',
    'THYM'	:'Thymoma',
    'UCEC'	:'Uterine_corpus_endometrial_carcinoma',
    'UCS' :'Uterine_carcinosarcoma',
    'UVM'	:'Uveal_melanoma'
    }
    

    
    #第一个列表是超级增强子，第二个列表是eQTL,第三个是gwas，第四个是enhancer，第五个是type
    index_table={
    'ACC':[  ['Adrenal Gland'],  ['Adrenal_Gland'],'adrenal',['kidney'],'aderal'],
    'BLCA':[['Bladder'],[],'bladder',['urinary_bladder'],'bladder'],
    'BRCA':[['HCC1954','MCF-7'],['Breast_Mammary_Tissue'],'breast',[],'breast'],
    'CESC':[['HeLa'],['Uterus.signifpairs'],'cervical',['uterus'],'cervical'],
    'CHOL':[[],[],'cholangio',['gallbladder'],'cholangio'],
    'COAD':[['Colon Crypt 1','Colon Crypt 2','Colon Crypt 3','Sigmoid Colon'],['Colon_Sigmoid','Colon_Transverse'],'colonic',['large_intestine'],'colon'],
    'ESCA':[['Esophagus'],['Esophagus_Gastroesophageal_Junction','Esophagus_Mucosa','Esophagus_Muscularis'],'esophageal',['esophagus'],'esophageal'],
    'GBM':[['Astrocytes'],['Brain_Cerebellar_Hemisphere'],'glioblastoma',['brain'],'brain'],
    'HNSC':[[],[],'head and neck squamous',['throat'],'head and neck squamous'],
    'KICH':[[],['Kidney_Cortex'],'kidney',['kidney'],'kidney'],
    'KIRC':[[],['Kidney_Cortex'],'kidney',['kidney'],'kidney'],
    'KIRP':[[],['Kidney_Cortex'],'kidney',['kidney'],'kidney'],
    'LAML':[['K562'],['Whole_Blood'],'acute myeloid leukemia',['blood'],'blood'],
    'LGG':[[],['Brain_Cerebellar_Hemisphere'],'glioma',['brain'],'brain'],
    'LIHC':[['HepG2'],['Liver'],'liver',['liver'],'liver'],
    'LUAD':[['H2171','GLC16','NCI-H82','NCI-H69','Lung'],['Lung'],'lung',['lung'],'lung'],
    'LUSC':[['H2171','GLC16','NCI-H82','NCI-H69','Lung'],['Lung'],'lung',['lung'],'lung'],
    'OV':[['Ovary'],['Ovary'],'ovarian',['female_gonad'],'ovarian'],
    'PAAD':[['Panc1','Pancreas','Pancreatic islets'],['Pancreas'],'pancreatic',['pancreas'],'pancreas'],
    'PCPG':[['Adrenal Gland'],[],'pheochromocytoma_and_paraganglioma',['kidney'],'kidney'],
    'PRAD':[['LNCaP'],['Prostate'],'prostate',['prostate_gland'],'prostate'],
    'READ':[['HCT-116','VACO 400','VACO 503','VACO 9m'],[],'rectum',['large_intestine'],'rectum'],
    'SARC':[[],[],'sarcoma',['adipose_tissue'],'sarcoma'],
    'SKCM':[['NHDF-Ad','NHEK'],['Skin_Sun_Exposed_Lower_leg'],'melanoma skin',['skin_of_body'],'skin'],
    'STAD':[['Gastric','Stomach Smooth Muscle'],['Stomach'],'gastric',['stomach'],'gastric'],
    'TGCT':[[],[],'testicular germ cell',['testis'],'testis'],
    'THCA':[[],['Thyroid'],'',['thyroid_gland'],'thyroid_gland'],
    'THYM':[['Thymus'],[],'thymus',['thymus'],'thymus'],
    'UCEC':[['HeLa'],['Uterus'],'uterine',['uterus'],'uterus'],
    'UCS':[['HeLa'],['Uterus'],'uterine',['uterus'],'uterus'],
    'UVM':[[],[],'melanoma',['eye'],'uveal_melanoma']
    }
    
    if(element is 'superenhancer'):
        return index_table[cancer][0],index_table[cancer][4]
    if(element is 'eQTL'):
        return index_table[cancer][1],index_table[cancer][4]
    if(element is 'gwas'):
        return index_table[cancer][2]
    if(element is 'enhancer'):
        return index_table[cancer][3],index_table[cancer][4]
    if(element is 'drivergene'):
        return abbreviate[cancer]


#eQTL,type=indextable('STAD','superenhancer')
#print(eQTL)
#print(type)
