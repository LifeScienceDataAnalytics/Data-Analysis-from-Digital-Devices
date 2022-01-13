#!/usr/bin/env python
# coding: utf-8

# # Functions and imports

# Imports and the functions that call the HI-VAE, modified from the paper version only to allow inputting s_codes and z_codes manually.

# In[6]:


import time
import numpy as np
import matplotlib
from matplotlib import pyplot as plt
from pandas.plotting import scatter_matrix
import os
import re
import pandas as pd
import numpy as np
from IPython.display import Audio
import seaborn as sns

import helpers # this is where the main training/decoding functions are, modified from teh original HIVAE main.py

#import warnings 
#warnings.filterwarnings('ignore') ########## NOTE: comment out for testing in case it's hiding problems

def set_settings(opts,nepochs=500,modload=False,save=True): # note: modload doesnt do anything right now, hardcoded in helpers.py
    'replace setting template placeholders with file info'
    inputf=re.sub('.csv','',opts['files'].iloc[0])
    missf=inputf+'_missing.csv'
    typef=inputf+'_types.csv'
    
    template = '--epochs NEPOCHS --model_name model_HIVAE_inputDropout --restore MODLOAD         --data_file data_python/INPUT_FILE.csv --types_file data_python/TYPES_FILE          --batch_size NBATCH --save NEPFILL --save_file SAVE_FILE        --dim_latent_s SDIM --dim_latent_z 1 --dim_latent_y YDIM         --miss_percentage_train 0 --miss_percentage_test 0         --true_miss_file data_python/MISS_FILE --learning_rate LRATE'
    
    # replace placeholders in template
    settings = re.sub('INPUT_FILE',inputf,template)
    settings = re.sub('NBATCH',str(opts['nbatch'].iloc[0]),settings)
    settings = re.sub('NEPOCHS',str(nepochs),settings)
    settings = re.sub('NEPFILL',str(nepochs-1),settings) if save else re.sub('NEPFILL',str(nepochs*2),settings)
    settings = re.sub('YDIM',str(opts['ydims'].iloc[0]),settings)
    settings = re.sub('SDIM',str(opts['sdims'].iloc[0]),settings)
    settings = re.sub('MISS_FILE',missf,settings) #if not 'medhist' in inputf else re.sub('--true_miss_file data_python/MISS_FILE','',settings)
    settings = re.sub('TYPES_FILE',typef,settings)
    settings = re.sub('SAVE_FILE',inputf,settings)
    settings = re.sub('LRATE',str(opts['lrates'].iloc[0]),settings)
    settings = re.sub('MODLOAD','1',settings) if modload else re.sub('MODLOAD','0',settings)
    
    return settings


# In[7]:


os.getcwd()


# In[8]:


sample_size=178
# get file list
files=[i for i in os.listdir('data_python/') if not '_type' in i and not '_missing' in i and i not in '.DS_Store']
files.sort()
print(files)


# In[10]:


best_hyper=pd.read_csv('results_Altoida.csv',  sep = ',')
best_hyper = best_hyper.drop(columns=['Unnamed: 0', 'Unnamed: 0.1'])
best_hyper.files = best_hyper.files.str.replace('\\_grid.*?\\_results', '')
best_hyper


# In[16]:


sample_size=178
# get file list
files=[i for i in os.listdir('data_python/') if not '_type' in i and not '_missing' in i and i not in '.DS_Store']
files.sort()
best_hyper = best_hyper.sort_values('files')
best_hyper = best_hyper.reset_index(drop=True)
#sds =[1,2,4,2,2,2,2,2,1,4,2,3,4,2,1] 
#sds =[1,2,4,2,2,2,2,2,1,4,2,3,4,1] 
#sds =[1,1,1,1,1,1,2,1,1,5,2,4,5,2,1] 

#sds =[1,1,1,1,1,1,2,1,1,5,2,4,3,2,1] 
##latest 10.02.2021
sds =[1,1,1,1,1,1,1,1,1,1,2,4,1,1,1,1,1,1,1,1,1,2,2, 2, 2, 2, 1, 4,2,3,4,1]

#sds = best_hyper['ydims']
sdims=dict(zip(files,sds))
if any(files!=best_hyper['files']):
    print('ERROR')
else:
    best_hyper['sdims']=sds
#sds =[1,1,1,1,1,1,2,1,1,4,2,3,3,2,1] 
best_hyper


# In[17]:


best_hyper['nbatch'] = best_hyper['nbatch'].astype(int) 
best_hyper['ydims'] = best_hyper['ydims'].astype(int) 
best_hyper['sdims'] = best_hyper['sdims'].astype(int) 


# In[18]:


best_hyper = best_hyper[['lrates', 'nbatch', 'wdecay', 'ydims', 'files', 'loss', 'sdims']]
best_hyper


# In[19]:


best_hyper.to_csv('best_hyper_ALTOIDA_processed.csv',index= False)


# # General settings

# sds is info about which files have what dimension of the "s_codes", that determine the number of mixture components in the "zcodes", our continuous embeddings used in the Bayes Net

# # Training

# In[20]:


import tensorflow.compat.v1 as tf
tf.disable_v2_behavior() 
for f in files:
    opts=dict(best_hyper[best_hyper['files'].copy()==f])
    settings=set_settings(opts,modload=False,save=True)
    helpers.train_network(settings)
wave = np.sin(2*np.pi*400*np.arange(10000*2)/10000)
Audio(wave, rate=10000, autoplay=True)


# # Get embeddings

# In[ ]:


dat=list()
dfs=list()
for f in files:
    # replace placeholders in template
    opts=dict(best_hyper[best_hyper['files'].copy()==f])
    opts['nbatch'].iloc[0]=sample_size
    settings=set_settings(opts,nepochs=1,modload=True,save=False)
    
    #run
    encs,encz,d=helpers.enc_network(settings)

    # make deterministic embeddings
    subj=pd.read_csv('python_names/'+re.sub('.csv','',f)+'_subj.csv')['x']
    sc=pd.DataFrame({'scode_'+re.sub('.csv','',f):pd.Series(np.array([i for i in encs])),'SUBJID':subj})  
    zc=pd.DataFrame({'zcode_'+re.sub('.csv','',f):pd.Series(np.array([i[0] for i in encz])),'SUBJID':subj})
    enc=pd.merge(sc, zc, on = 'SUBJID')
    
    # save out individual file's metadata
    enc.to_csv('Saved_Networks/'+re.sub('.csv','',f)+'_meta.csv',index = False)
    dfs.append(enc)
    dat.append(d)

# join metadata
enc_vars=[pd.read_csv('Saved_Networks/'+re.sub('.csv','',f)+'_meta.csv') for f in files]
meta=helpers.merge_dat(enc_vars)
meta[meta.columns[['Unnamed' not in i for i in meta.columns]]].to_csv('metaenc.csv',index= False)

dat_dic=dict(zip(files,dat))


# In[ ]:


meta = pd.read_csv("metaenc.csv", sep = ',')


# In[ ]:


print(len(meta.columns))


# In[ ]:


fig = scatter_matrix(
    meta[meta.columns.drop(list(meta.filter(regex='SUBJID|scode_')))],
    figsize  = [20, 20],
    marker   = ".",
    s        = 10,
    diagonal = "kde"
)
for ax in fig.ravel():
    ax.set_xlabel(re.sub('_VIS|zcode_','',ax.get_xlabel()), fontsize = 20, rotation = 90)
    ax.set_ylabel(re.sub('_VIS|zcode_','',ax.get_ylabel()), fontsize = 20, rotation = 90)
    
plt.suptitle('HI-VAE embeddings (deterministic)',fontsize=20)


# # RP decoding (Reconstruction)

# In[ ]:


meta = pd.read_csv('metaenc.csv')
sample_size = 105
recon=list()
recdfs=list()
for f in files:
    # replace placeholders in template
    opts=dict(best_hyper[best_hyper['files'].copy()==f])
    opts['nbatch'].iloc[0]=sample_size
    settings=set_settings(opts,nepochs=1,modload=True,save=False)
    
    #run
    zcodes=meta['zcode_'+re.sub('.csv','',f)]
    scodes=meta['scode_'+re.sub('.csv','',f)]
    rec=helpers.dec_network(settings,zcodes,scodes)
    recon.append(rec)
    
    subj=pd.read_csv('python_names/'+re.sub('.csv','',f)+'_subj.csv')['x']
    names=pd.read_csv('python_names/'+re.sub('.csv','',f)+'_cols.csv')['x']
    recd=pd.DataFrame(rec)
    recd.columns=names
    recd['SUBJID']=subj
    recdfs.append(recd)
    
recon_dic=dict(zip(files,recon))

data_recon=helpers.merge_dat(recdfs)
data_recon.to_csv('reconRP.csv',index=False)


# In[ ]:





# In[ ]:




