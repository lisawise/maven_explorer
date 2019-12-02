## version 20190502
    # fixed bug in ms2 fa tail dictionary mapping

## version 20190320
    # refined ms2 fragment clustering and aggregation into function

## version 20190312
    # changed naming of spectra files to include experiment conditions

import pandas as pd
import numpy as np
from pyteomics import mgf, mzxml, mass
import os
import matplotlib.pyplot as plt
## for aggregation
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy import cluster
from scipy.spatial.distance import pdist


## variables

spectra_threshold = 10e-6
#threshold to identify fragment
fragment_threshold = 6e-6
#threshold to mask intensities (set to 0 for no masking)
int_threshold = 0
# set number of top intensities (set to 10000 for no masking)
int_range = 10000
# set fraction of intensities (set to 1 for no masking)
topfrac = 1
# set rt threshold in minutes
rt_thresh = 1

def test_function():
    return os.curdir

def mzxml_import_untargeted(mdf, mzxml_dir, compound_col = 'compound'):

    # create file_list from mdf using regexp search for date followed by _
    file_list = mdf.filter(regex='\d{4,8}[_]').columns
    file_list = [s + '.mzxml' for s in file_list]
    print(file_list)
    
    ### below commented code will search all mzxml files in a folder
    #file_list = os.listdir(mzxml_dir)
    #file_list = [f for f in file_list if ('.mzXML' in f)]
    #print(file_list)

    dataframe_list = []
    compound = compound_col

    compound = mdf[compound_col]

    medMz = mdf['medMz']
    medRt = mdf['medRt']



    ## mdf is a matrix with the relevant info for getting list
    mdf = pd.concat([compound, medMz, medRt], axis=1)
    
    
    for file_to_open in file_list:

        data_path_1 = os.path.join(mzxml_dir, file_to_open)

        ### Search spectra in .mzxml file for precursor ions matching Maven list (in mdf matrix)

        # Create empty list to store results
        results = []
        try:
            # Load mass spec file using mzxml - returns a generator with results
            ms_file = mzxml.read(data_path_1, read_schema=True, iterative=False, use_index=False, dtype=None)
            
            print (file_to_open+' : Processing...')

            # Iterate through elements in mass spec data generator
            for i, spectra in enumerate(ms_file):

                # iterate through m/z values to look for in each mass spec peak entry
                for k, mz_target in enumerate(mdf[mdf.columns[1]]):
                    #rt_target = mdf.iloc[k,2]

                    # 'precursorMz' in keys indicate a 2nd-order mass spec result with a 2nd smaller dictionary as a value
                    if 'precursorMz' in spectra.keys():
                        mz_measured = spectra['precursorMz'][0]['precursorMz']

                        # if one of the recorded peaks is within 20ppm of the reference value,
                        # extract the information in that entry and append to results

                        # combine m/z array and intensity array (mzarray), also take relative intensity and sort descending (mzarray_norm)
                        if (abs((mz_measured-mz_target)/mz_target)<spectra_threshold) and (abs(mdf.iloc[k,2]-spectra['retentionTime'])<rt_thresh) :
                        #if abs((mz_measured-mz_target)/mz_target)<spectra_threshold and abs(spectra['retentionTime']-mdf.iloc[k,2])<rt_thresh :
                            q = spectra['intensity array']
                            maxvali = np.amax(q,axis=0)
                            q_norm = q/maxvali
                            r = spectra['m/z array']
                            mzarray = np.array(list(zip(r,q)))
                            mzarray_norm = np.array(list(zip(r,q_norm)))
                            mzarray_norm = mzarray_norm[(-mzarray_norm[:,1]).argsort()]

                            result_dict = {'precursorMz_m':spectra['precursorMz'][0]['precursorMz'],
                                          'retentionTime_m':spectra['retentionTime'],
                                          'mzarray':mzarray,  
                                          'mzarray_norm':mzarray_norm,
                                          'compound':mdf.iloc[k]['compound']}
                            results.append(pd.Series(result_dict))

            # concatenate list of results into a single dataframe
            results = pd.DataFrame(results)

            # preview results
            #print('Done')
            #results.head()

            if len(results)==0: continue 

            #index mdf and results to compound and then combine
            results_c = results.set_index('compound')
            mdf_c = mdf.copy()
            mdf_c = mdf_c.set_index('compound')

            #this combines
            data = mdf_c.join(results_c)

            # look for specific fragments in each m/z array    

            #for ref in ref_array:
            #    data[str(ref)] = data['mzarray_norm'].apply(check_in_array, ref_value = ref)
            data["sample"] = file_to_open

            dataframe_list.append(data)
            group_df = pd.concat(dataframe_list)
            ## NOTE: when same compound names are included with different RTs, some spectra match to wrong RT
            ## below filter removes spectra matched with incorrect RT
            ## remove spectra that are not within rt range
            group_df['diff'] = group_df['medRt'] - group_df['retentionTime_m']
            group_df = group_df.loc[~(abs(group_df['diff']) >=1)]
            group_df = group_df.drop(['diff'], axis=1)

        except:
            print(file_to_open+' not found')
            continue
            
    return group_df

## make a copy, reset index and remove na

def df_remove_na(group_df):
    data = group_df.reset_index().copy()
    filtered_data = data.dropna()

    return filtered_data

### make mdf matrix from maven file (have to change groupid to compound)
def make_mdf(data_dir, maven):

    data_path_2 = os.path.join(data_dir, maven)

    mdf = pd.read_csv(data_path_2)
    mdf = mdf.drop(['compound'], axis=1)
    mdf = mdf.rename(columns={'groupId': 'compound'})

    return mdf

def check_in_array(x, ref_value, fragment_threshold=6e-6):
    try:
        return np.min(abs(x[:,0] - ref_value))/ref_value < fragment_threshold
    except:
        return np.nan
    
def intensity_threshold(x, int_threshold = int_threshold):
    try:
        #x = data.iloc[10]['mzarray_norm']
        xi = x[:,1]
        mask = xi > int_threshold
        xim = xi[mask]
        yim = x[:,0][mask]

        return np.array(list(zip(yim, xim)))
            #x[:,0]
    except:
        return np.nan
    
def intensity_topx(x, topnum = int_range):
    try:
        if topnum > len(x):
            topnum == len(x)
        return np.array(x[0:topnum, :])
    except:
        return np.nan
    
def intensity_topfrac(x, topfrac=topfrac):
    try:
        y = np.round(len(x)*topfrac, 0).astype(int)
        return np.array(x[0:y, :])
    except:
        return

## This function searches mdf for all clusters with specified m/z in range set by ppm  
def mz_search(mz, mdf_mean_wide, ppm = 20):

    return mdf_mean_wide.loc[(mdf_mean_wide['mean_medMz']>(mz-mz*ppm/1000000))&(mdf_mean_wide['mean_medMz']<(mz+mz*ppm/1000000))]
    

 
def apply_mdf(mdf, filtered_data):
    include = filtered_data['compound'].drop_duplicates()

    # make a new column, spectra, where = 1 if spectra present in ANY sample for that compound, 0 otherwise
    mdf.insert(0, 'spectra', np.where(mdf['compound'].isin(include), 1, 0) , allow_duplicates=False)

    mdf2 = mdf.drop(['label', 'metaGroupId', 'maxQuality', 'note',
                    'compoundId', 'expectedRtDiff', 'ppmDiff', 'parent', 'goodPeakCount'], axis=1)

    return mdf2

def apply_mdf_melt(mdf, filtered_data):
## apply list to mdf (all data)
    
    exclude = make_exclude(filtered_data)
    mdff = mdf[~mdf['compound'].isin(exclude)]
    
    include = filtered_data['compound'].drop_duplicates()
    mdff2 = mdff[mdff['compound'].isin(include)].reset_index()


    ## drop unnecessary columns
    mdff2 = mdff2.drop(['label', 'metaGroupId', 'maxQuality', 'note',
                    'compoundId', 'expectedRtDiff', 'ppmDiff', 'parent', 'index'], axis=1)



    ## melt data to get signal for each sample x each peak (keep mz and rt and fragments as ID vars)

    mdf_melt = pd.melt(mdff2, id_vars=['compound', 'goodPeakCount', 'medMz', 'medRt'])
    
    mdf_melt['treatment'] = mdf_melt['variable'].str.extract(r'(b11_m|b11_v|nt_m|nt_v|blank|dmem)', expand = False)
    #mdf_melt['dilution'] = mdf_melt['variable'].str.extract(r'(dil)', expand = False)
    mdf_melt['phase'] = mdf_melt['variable'].str.extract(r'(supe|pellet)', expand = False)
    mdf_melt['charge'] = mdf_melt['variable'].str.extract(r'(pos|neg)', expand = False)
    mdf_melt['lc_col'] = mdf_melt['variable'].str.extract(r'(f5|hilic)', expand = False)
    #mdf_melt['rep'] = mdf_melt['variable'].str.extract(r'(_a_|_b_)', expand = False)
    
    return mdf_melt  

## This function returns a CSV and spectra plot of ms2 data for a specified compound.
## Options include to limit fragments of specified relative intensity (spectensity), save files to output_dir,
## and return specific index or all spectra
def spectra_plot_compound(df, compound, base_name, output_dir, index=0, spectensity=0.01, 
                            compound_col='compound', save = False,
                            run_all = False):

    df = df.copy()
    df = df.reset_index(drop = True)
    df = df.loc[df[compound_col] == compound].reset_index(drop = True)
    print('# of spectra: '+str(len(df)))
    
    count = 0
    index_loop = 0
    #for index_loop in range(len(df)):
    while count < len(df):
         
        if run_all == False:
            index_loop = index
            count = count+len(df)
        #print(index_loop)
        array = df.iloc[index_loop]
        mz_m = round(array['precursorMz_m'], 4)
        #print(mz_m)
        rt_m = round(array['retentionTime_m'], 4)
        sample = array['sample']
        #print(rt_m)
        print(str(compound)+': mz = '+str(mz_m)+', RT = '+str(rt_m)+' ('+str(sample)+')')
        array = array['mzarray_norm']
        arraypd = pd.DataFrame(data=array, columns=['m/z','relative intensity'])
        arraypd.insert(0, sample, pd.Series([]))
        array = intensity_threshold(array, spectensity)
        #return array

        ## m/z (x)
        x = array[:,0]

        ## relative intensity (y)
        y = array[:,1]


        fig, ax = plt.subplots()
        rects1 = ax.bar(x, y, width=0, linewidth=2, edgecolor='black', alpha=0.4)

        plt.title("["+str(index_loop)+'] Experimental spectra for compound: '+compound
                    +'\n m/z = '+str(mz_m)+' RT = '+str(rt_m)+'\n sample: '+sample)
        plt.xlabel('m/z, fragmentation spectra')
        plt.ylabel('Intensity, rel. units')
        #plt.ylim(0,1.5)
        #plt.xlim(0, 250)

        3## adds m/z labels
        autolabel(rects1, ax)
        plt.show()

        ## file name
        compound = compound.replace(':', '-')
        compound = compound.replace('(', '-')
        compound = compound.replace(')', '')
        spectra_name = compound+'_index_'+str(index_loop)+'_'+base_name
        ## below will change non alphanumeric characters to _ in filename
        spectra_name = "".join([x if x.isalnum() else "_" for x in spectra_name])
        print(spectra_name)
        if save == True:
            #fig.savefig(os.path.join(output_dir, spectra_name+'.pdf'), bbox_inches='tight')
            arraypd.to_csv(os.path.join(output_dir, spectra_name+'.csv'), index=False)
            #np.savetxt(spectra_name+'.csv', array, delimiter=",")
        print(array)
        count = count+1
        index_loop = index_loop+1

### show scatterplot of signal (input is group_id for compound)
def signal_plot(mdf_melt, group_id):
    compound = mdf_melt.loc[(mdf_melt['compound'] == group_id)].reset_index()
    sns.swarmplot(x = 'treatment', y = 'value', data = compound)
    plt.title('medMz: '+str(compound.loc[0, 'medMz'])+" ; groupid: "+str(group_id))
    plt.show()

def signal_plot_clust(mdf_melt, group_id):
    compound = mdf_melt.loc[(mdf_melt['new_clust'] == group_id)].reset_index()
    sns.swarmplot(x = 'treatment', y = 'value', data = compound, hue='date')
    plt.title('medMz: '+str(compound.loc[0, 'mean_medMz'])+" ; cluster: "+str(group_id))
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    plt.show()

### adds text labels to figures (external)    
def autolabel(rects, ax):
    # Get y-axis height to calculate label position from.
    (y_bottom, y_top) = ax.get_ylim()
    y_height = y_top - y_bottom

    for rect in rects:
        height = rect.get_height()
        if height > 0.01:
            mz = rect.get_x()
        # Fraction of axis height taken up by this rectangle
            p_height = (height / y_height)

        # If we can fit the label above the column, do that;
        # otherwise, put it inside the column.
            if p_height > 0.95: # arbitrary; 95% looked good to me.
                label_position = height - (y_height * 0.05)
            else:
                label_position = height + (y_height * 0.01)

            ax.text(rect.get_x() + rect.get_width()/2., label_position,
                "{0:.4f}".format(mz), ha='center', va='bottom')

### aggregation functions

def flatten_list(l):
  out = []
  for item in l:
    if isinstance(item, (list, tuple)):
      out.extend(flatten_list(item))
    else:
      out.append(item)
  return out

def partition(values, indices):
    idx = 0
    for index in indices:
        sublist = []
        #print(type(idx))
        while idx < len(values) and values[idx] < index:
            sublist.append(values[idx])
            idx +=1
        if sublist:
            yield sublist

def g(x):
    return pd.Series(dict(mean_medMz = x['medMz'].mean(), 
                           mean_Rt = x['medRt'].mean(),
                            spectra = x[x['spectra']>0].spectra.count()))





## This function combines untargeted data from different replicates by hierarchical aggregation along m/z and RT 
def aggregate_mz(X, indices):
    
    df = pd.concat(X, ignore_index=True, sort=False).sort_values(by=['medMz']).reset_index(drop=True)
    
    #df = df.loc[df['medMz']<1500]
    
    X2 = df[['medMz']].values.tolist()
    print(len(X2))
    
    
    X3 = flatten_list(X2)

    print("Minimum m/z val: "+str(np.round(min(X2), 4)))
    print("Maximum m/z val: "+str(np.round(max(X2), 4)))


    X4 = list(partition(X3, indices))

    # check ppm differences among sublist splits (should be less than 5)
    print('ppm difference check')
    for i in range(len(X4)-1):
        #while i < (len(X4)-1):
        max_x = max(X4[i])
        min_x = min(X4[i+1])
        print(str(i)+'| delta_ppm = '+str(abs((max_x-min_x)/((max_x+min_x)/2)*float(1000000)))+' ('+str(max_x)+', '+str(min_x)+')')
    
    y = np.array([np.expand_dims(np.array(xi),-1) for xi in X4])
    print('begin aggregation')
    max_clust=0
    dataframe_list = []
    for idx, subarray in enumerate(y):
        
        dist = ((max(subarray)-min(subarray))/2+min(subarray))*float(10/1000000)
        Y = pdist(subarray, metric='euclidean')

        Z = linkage(Y, 'ward')
        #print(str(len(Z)))
        
        clusters = fcluster(Z, dist, criterion='distance')
        ## the addition of max_clust here makes each cluster unique as loop progresses
        clusters = np.array(list(map(lambda x:x+max_clust, clusters)))
        #print(clusters)
        max_clust = max(clusters)
        
        ### concatenate together the m/z values and the clusters and turn to dataframe (convert clusters to "int")
        x2 = np.vstack((subarray.flatten(), clusters))
        x3 = pd.DataFrame(x2.transpose(), columns = ['medMz_2', 'mzclusters'])
        x3['mzclusters'] = x3['mzclusters'].astype(int)
        
        #print(x3['mzclusters'])
        
        ## through each iteration, append to "dataframe_list" and then make into full dataframe
        dataframe_list.append(x3)
        group_df = pd.concat(dataframe_list).reset_index(drop=True)
        print(str(idx)+'| dist criterion:'+str(dist)+", unique clust="+str(len(np.unique(clusters)))+', total compounds='+str(len(subarray)))
        
    df2 = pd.concat([group_df['mzclusters'], df], axis=1, sort=False).sort_values(by = ['mzclusters'])
    
    appended_data = []
    for mzclusters, new_df in df2.groupby('mzclusters'):
        
        X1 = new_df[['medRt']].values
        if len(X1) >1:
            Y = pdist(X1,metric='euclidean')

            Z = linkage(Y, 'ward')
            clusters = fcluster(Z, 1.5, criterion='distance')
            rtclusters = pd.Series(clusters, name='rtclusters', index=new_df.index)
            #x3 = pd.Series(clusters, name='rtclusters', index=new_df.index)
        else: rtclusters = pd.Series(np.array([1]), name='rtclusters', index=new_df.index)

        appended_data.append(rtclusters)
        #print(new_df)
    
    appended_data = pd.concat(appended_data, axis=0)
    
    df2['rtclusters'] = appended_data
    df2['new_clust'] = df2['mzclusters'].astype(str)+'mz_'+df2['rtclusters'].astype(str)+'rt'
    df2 = df2.sort_values(by=['new_clust'])
    
    ## df with identifying information for each mz/rt cluster, averaged, and count of each spectra
    df_new = df2.groupby('new_clust').apply(g)#.reset_index()
    ## remove extraneous columns (averaged in other df)
    df2 = df2.drop(columns=['compound', 'medMz', 'medRt', 'spectra', 'mzclusters', 'rtclusters'])
    
    ## melt df with mz/rt cluster as id, then pivot to re-form with all data columns
    df_melt = pd.melt(df2, id_vars=['new_clust'])
    df_group = df_melt.groupby(['new_clust', 'variable']).mean().reset_index()
    df_group = df_group.pivot(index='new_clust', columns = 'variable', values = 'value').reset_index().set_index('new_clust')
    
    
    result = pd.concat([df_new, df_group], axis=1)
    print(len(result))
    
    
    return result

### MS2 fragment aggregation functions

def aggregate_ms2(filtered_data, indices = [200.01, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1600], ppm_frag = 60):
    #ppm_frag = 80
    dataframe_list=[]
    for name, group in filtered_data.groupby(['compound']):
        new_df = group.reset_index().drop(columns=['index', 'mzarray']).reset_index()
        #print(len(new_df))

        dataframe_list.append(new_df)
        group_df = pd.concat(dataframe_list).reset_index(drop=True)

    group_df = group_df.rename(index=str, columns={"index": "spec_index"})
    df_frag2 = group_df[['spec_index','compound', 'medMz', 'medRt', 'mzarray_norm', 'sample']].copy()
    df_frag2

    rows = []
    _ = df_frag2.apply(lambda row: [rows.append([row['spec_index'], row['compound'], row['medMz'], row['medRt'], nn, row['sample']]) 
                             for nn in row.mzarray_norm[:,0:2]], axis=1)
    df_new = pd.DataFrame(rows, columns=df_frag2.columns)#.set_index(['compound'])

    df_new['mz frag'] = df_new['mzarray_norm'].apply(lambda x: x[0])
    df_new['intensity frag'] = df_new['mzarray_norm'].apply(lambda x: x[1])
    df_new = df_new.sort_values(by = ['mz frag']).reset_index(drop=True).drop(labels = 'mzarray_norm', axis = 1)

    print('fragment matrix length: '+str(len(df_new)))
    
    X2 = df_new[['mz frag']].values.tolist()

    X3 = flatten_list(X2)

    print("Minimum m/z val: "+str(np.round(min(X2), 4)))
    print("Maximum m/z val: "+str(np.round(max(X2), 4)))

    ## set indices to split along
    #indices = [0, 100.055, 200.2, 300, 400.18, 500, 600, 700, 800, 900]
    #indices = [200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200]
    # generate list of sublists splitting m/z along set indices
    # X4[0] = list of m/z with values from 0 to 100 
    X4 = list(partition(X3, indices))

    for i in range(len(X4)-1):
        #while i < (len(X4)-1):
        max_x = max(X4[i])
        min_x = min(X4[i+1])
        print(str(i)+'| delta_ppm = '+str(abs((max_x-min_x)/((max_x+min_x)/2)*float(1000000)))+' ('+str(max_x)+', '+str(min_x)+')')
        
    # check ppm differences among sublist splits (should be less than 5)
    y = np.array([np.expand_dims(np.array(xi),-1) for xi in X4])
    max_clust=0
    dataframe_list = []
    for idx, subarray in enumerate(y):
        
        dist = ((max(subarray)-min(subarray))/2+min(subarray))*float(ppm_frag/1000000)
        Y = pdist(subarray, metric='euclidean')

        Z = linkage(Y, 'ward')
        #print(str(len(Z)))
        
        clusters = fcluster(Z, dist, criterion='distance')
        ## the addition of max_clust here makes each cluster unique as loop progresses
        clusters = np.array(list(map(lambda x:x+max_clust, clusters)))
        #print(clusters)
        max_clust = max(clusters)
        
        ### concatenate together the m/z values and the clusters and turn to dataframe (convert clusters to "int")
        x2 = np.vstack((subarray.flatten(), clusters))
        x3 = pd.DataFrame(x2.transpose(), columns = ['medMz_2', 'mzclusters'])
        x3['mzclusters'] = x3['mzclusters'].astype(int)
        
        ## through each iteration, append to "dataframe_list" and then make into full dataframe
        dataframe_list.append(x3)
        group_df = pd.concat(dataframe_list).reset_index(drop=True)
        print(str(idx)+'| dist criterion:'+str(dist)+", unique clust="+str(len(np.unique(clusters)))+', total compounds='+str(len(subarray)))

    df2 = pd.concat([group_df['mzclusters'], df_new], axis=1, sort=False)
    return df2

def groupby_spectra(df2):
    dataframe_list = []
    for name, group in df2.groupby(['compound', 'spec_index']):
        new_df = group.groupby(["mzclusters"]).agg({

                           'compound': {'compound': 'first'},
                            'spec_index': {'spec_index':'first'},
                            'sample': {'sample':'first'},
                           'medRt': {'medRt':'first'},
                            'medMz': {'medMz':'first'}, 
                            'mzclusters': {'fragment_cluster' : 'first'},
                           'mz frag': 
                                      {'frag_mw': 'mean', 'count': 'count'},
                            'intensity frag': {'avg intensity': 'mean'}
                                                            })
        new_df.columns = new_df.columns.droplevel(0)
        #print(new_df)
        new_df = new_df.sort_values(by = 'count', ascending = False)
        dataframe_list.append(new_df)
        group_df_frag2 = pd.concat(dataframe_list).reset_index(drop=True)
    return group_df_frag2
        
def groupby_ms2frag(df):
    dataframe_list = []
    for name, group in df.groupby(['compound']):
        ## total spectra needs to be outside the fragment cluster aggregation
        group['total'] = group['spec_index'].nunique()
        #print(group)
        new_df = group.groupby(["fragment_cluster"]).agg({

                           'compound': {'compound': 'first'},
                           'medRt': {'medRt':'first'},
                            'medMz': {'medMz':'first'}, 
                            'total' : {'total': 'first'},
                            'fragment_cluster': {'fragment_cluster' : 'first'},
                           'frag_mw': 
                                      {'frag_mw': 'mean', 'count': 'count'},
                            'avg intensity': {'avg intensity': 'mean'}
                                                            })
        new_df.columns = new_df.columns.droplevel(0)
        #print(new_df)
        new_df = new_df.sort_values(by = 'count', ascending = False)
        dataframe_list.append(new_df)
        df = pd.concat(dataframe_list).reset_index(drop=True)
        

    print(len(df))
    df['rel freq'] = df['count']/df['total']
    return df
    
def create_frag_dict(df2, file, ppm = 10):
    cluster_mz = df2.groupby('mzclusters').mean()['mz frag']
    #cluster_mz.get_loc((cluster_mz > 255.23)&(cluster_mz < 255.24))
    #ppm = 10000
    for index, row in file.iterrows():
        mz = row['mz']
        #print(mz)
        name = row['name']
        #print(name)

        cluster_mz_x = cluster_mz[cluster_mz.between(mz-ppm/1000000*mz, mz+ppm/1000000*mz)]
        
        #print('cluster for '+name+' is: '+
        if cluster_mz_x.empty:
            continue
        elif cluster_mz_x.size == 1: 

            print('cluster for '+name+'(m/z = '+str(round(mz,4))+') is: '+str(cluster_mz_x.index[0])+\
                 ', m/z = '+str(cluster_mz_x.values[0].round(4)))
            file.iloc[index, 2] = cluster_mz_x.index[0].astype(int)
            file.iloc[index, 3] = cluster_mz_x.values[0]
        elif cluster_mz_x.size > 1:
            print('multiple cluster values for '+name)
            for i in range(cluster_mz_x.size):
                    print(str(i)+" cluster "+\
                          str(cluster_mz_x.index[i].astype(int))+\
                         ", m/z is: "+str(cluster_mz_x.values[i].round(4)))

    file = file.loc[file['mzclusters']>0]  
    file_dict = dict(zip(file['mzclusters'], file['name']))
    return file_dict