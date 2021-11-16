
import re
import pandas as pd 

with open('/home/a0073895/Athlates_2014_04_26/protein_annotation/gencode.v19.pc_translations.fa') as file:
    buf = file.read()
    ## generate a dict that uses the transcript id as key 
    seqDict = dict(zip([i.strip('>').split('|')[0].split('.')[0] for i in buf.split()[::2]], buf.split()[1::2]))
    
    ## key is used to convert ref seq ids to Ensembl protein ids
    key = pd.read_table('/home/a0073895/Athlates_2014_04_26/protein_annotation/gencode.v24.metadata.RefSeq',header=None)
    key.columns = ['Ensembl','RefM','RefP']
    key['Ensembl'] = key.Ensembl.apply(lambda x: x.split('.')[0])
    file.close()
    
def remove_amb(data):
    ### Uses the IUB ambiguity codes to convert bases
    #data.columns = ['chr','pos','pos1','ref','var']
    
    ambiguityDict = {'var':{'Y':'TC','W':'AT','K':'TG','R':'GA','M':'CA','S':'GC'}}

    ### replaces ambiguity bases to simpleBases , ATCG
    data = data.replace(to_replace=ambiguityDict)

    ### make a copy of the Cons column
    #data.loc[:,'var1'] = data.loc[:,'var']

    ### replaces ref bases in var column, giving only var bases
    #data['var1'] = data.var1.replace(refDict,'')
    for num, row in data.iterrows():
        data.ix[int(num),'var'] = data.ix[int(num),'var'].replace(data.ix[int(num),'ref'],'')

    ### Select columns and reorder positions
    #data = data.loc[:,['chr','pos','ref','var1','sample']]
    #data.columns = ['chromosome', 'position', 'genotype_of_normal', 'mutant_allele']
    return data

def get_original(df):
    try:
        """Generates a 17 amino acid reference sequence by using the position mutated predicted by oncotator"""
        ## get the amino acid for the transcript
        #aaSeq = seqDict[df.Transcript_id.split(r'.')[0]]
        aaSeq = seqDict[df.Ensembl_id]
        ## use the protein change position to find the mutated aa 
        pos = int(re.findall('\d+',df.Protein_Change)[0])

        ## generate a 17 mer using the reference amino acid 
        seventeenMer = aaSeq[pos-9:pos+8]
    
    except KeyError as e:
        seventeenMer = ''
        
    return seventeenMer

def get_mutated(df):
    """Generates a mutated 17 amino acid sequence from get_original """
    seventeenMer = df.Normal 
    mutated_aa = seventeenMer[:8]+df.Protein_Change[-1] + seventeenMer[9:]
    return mutated_aa

def get_AA(series):
    """Generates a set of 9 amino acid sequences from a 17mer sequence"""
    seventeenMer = series
    nineMers = [seventeenMer[i[0]:i[1]] for i in zip(range(0,9),range(9,len(seventeenMer)+1))]
    return nineMers

def predict_neoantigen(mutationCallFile):
    """Predicts a set of 9 amino acid neoantigens using 17 aa sequences from reference transcripts
    and predictions from oncotator"""
    
    ## get the sample name
    sampleName = mutationCallFile.strip('.txt').split('/')[-1]
    
    ## read input file and format columns for oncotator 
    inputFile = pd.read_table(mutationCallFile,index_col=False)
    
    ## select columns that are useful and rename them
    #output = inputFile[['Chromosome','Position','Position','ReferenceBase','ConsensusBase']]
    output = inputFile[['Chr','Pos','Pos','Ref','Cons']]
    output.columns = ['chr','start','end','ref_allele','alt_allele']
    
    ## output into a tsv file for running oncotator
    output.to_csv(sampleName+'_snv.txt',index=False,sep='\t')
    
    
    snvName = sampleName+'_snv.txt'
    ## annotated file already exists
    if os.path.isfile(snvName):
        pass
    
    else:
        ## run oncotator 
        !/home/a0073895/homeEnv/bin/oncotator -v --db-dir /gpfs/public/tools/oncotator_v1_ds_Jan262014 {sampleName+'_snv.txt'} {sampleName+'_annotated.txt'} hg19
    
    ### for testing 
    #     annotatedDf =  pd.read_table('./snv_data/130T_snv_annotate.txt',header=3)
    
    ## read the oncotator annotated output 
    annotatedDf = pd.read_table(sampleName+'_annotated.txt',header=3)
    
    ## select only missense variants, should be extended to frame shift and splice site mutations 
    cols = ['Chromosome','Start_position','transcript_id','Hugo_Symbol','Protein_Change']
    annotatedDf = annotatedDf.ix[annotatedDf.Variant_Classification=='Missense_Mutation',cols]
    annotatedDf.ix[:,'Chromosome'] = 'chr'+ annotatedDf.Chromosome.map(str) 
    annotatedDf.columns = ['Chromosome','Position','Transcript_id','Hugo_Symbol','Protein_Change']
    
    ## get the original 17mer sequence, the mutated 17mer aa sequence and a list of 9mer slices
    annotatedDf.ix[:,'Normal'] = annotatedDf.apply(get_original,axis=1 )
    annotatedDf.ix[:,'Mutated'] = annotatedDf.apply(get_mutated,axis=1)
    annotatedDf.ix[:,'refNonamer'] = annotatedDf.Normal.apply(get_AA)
    annotatedDf.ix[:,'mutNonamer'] = annotatedDf.Mutated.apply(get_AA)
    
    ## combine the table as output 
    originalDf = inputFile[['Chr','Pos','Ref','Cons','Ref Context','Var Context']]
    originalDf.columns = ['Chromosome','Position','Ref','Var','Ref Context','Var Context']
    combinedDf = annotatedDf.merge(originalDf,on=['Chromosome','Position'])
    combinedDf.to_csv('./'+sampleName+'_predicted.txt',sep='\t',index=False)
    
    refBuffer= ''
    for row in combinedDf['refNonamer']:
        refBuffer+='\n'.join(row)
        refBuffer+='\n'

    with open('./'+sampleName+'_ref_antigen.txt','w') as refFile:
        refFile.write(refBuffer)
        refFile.close()
        
        
    mutBuffer= ''
    for row in combinedDf['mutNonamer']:
        mutBuffer+='\n'.join(row)
        mutBuffer+='\n'

    with open('./'+sampleName+'_mut_antigen.txt','w') as mutFile:
        mutFile.write(mutBuffer)
        mutFile.close()

    return combinedDf


## read the oncotator annotated output 
#annotatedDf = pd.read_table(sampleName+'_annotated.txt',header=3)

## select only missense variants, should be extended to frame shift and splice site mutations 
#annotatedDf = annotatedDf.ix[annotatedDf.Variant_Classification=='Missense_Mutation',['Chromosome','Start_position','transcript_id','Hugo_Symbol','Protein_Change']]
#annotatedDf.ix[:,'Chromosome'] = 'chr'+ annotatedDf.Chromosome.map(str) 

def predict_neo_antigenV3(input,sampleName):
    
    ## input should be a formatted file with the few columns 
    ## ['Chr','Start','End','Ref','Alt','Gene','Sample','Type','Transcript','ProteinChange']
    
    annotatedDf = input[['Chr','Start','End','Ref','Alt','Gene','Sample','Type','Transcript','Protein_Change']]
        
    logBuffer = []

    annotatedDf['Ensembl_id'] = annotatedDf.Transcript.apply(lambda x : list(set(key.ix[key.RefM.str.contains(x),'Ensembl'].values.tolist()))) 
    annotatedDf['id_len'] = annotatedDf.Ensembl_id.apply(lambda x: len(x))

    ## check if any gene without NCBI-ensembl id conversion are left out 
    noHits = annotatedDf.ix[annotatedDf.id_len==0,'Gene'].unique()
    for gene in noHits:

        if annotatedDf.ix[(annotatedDf.Gene==gene)&(annotatedDf.id_len>=1),:].empty:
            print gene
        else:
            pass

    ## for multiple transcripts, milk what i can milk, use only those that have ensembl ids in the seqDict 

    seriesList = []

    for row, df in annotatedDf[annotatedDf.id_len>1].iterrows():

        try:
            ## check if ensembl ids are in seqDict

            withKeys = [i for i in df.Ensembl_id if i in seqDict.keys()]
            if withKeys ==[]:
                logBuffer.append("no ensembl ids found for %s, ids being: %s" % (df.Gene,df.Ensembl_id))

            else:
                for ids in withKeys:
                    #df['Protein_Change'] = df.Protein_Change.astype(str)
                    matchPos = re.match('(\w)(\d{0,6})\w',df.Protein_Change)
                    pos= int(matchPos.group(2))
                    ref = matchPos.group(1)

                    ## check if position in protein sequence is the same as the Protein_Change aa
                    if ref == seqDict[ids][pos-1]:

                        ## put it back into a list to be joined to id_len ==1 
                        df['Ensembl_id'] = ids
                        seriesList.append(df)

                    else:
                        logBuffer.append("AA residues dont match Protein_Change for %s %s , AA being: %s %s" % (df.Gene, df.Ensembl_id,ref,seqDict[ids][pos-1]))

        except Exception as e:
            print str(e) +'\n'
            
    ## change ensembl_id from list to string, and combine samples w 1 ensembl id with the fixed rows
    annotatedDf.ix[annotatedDf.id_len==1,'Ensembl_id'] = annotatedDf.ix[annotatedDf.id_len==1,'Ensembl_id'].apply(lambda x: x[0])
    
    if seriesList != []:
        ## fixedDf has each individual transcript as 1 row, with non existent keys removed and mismatch AA positions 
        ## removed as well
        fixedDf = pd.concat(seriesList,axis=1).T.drop_duplicates('Ensembl_id')  
        joinedDf = annotatedDf.ix[annotatedDf.id_len==1,:].append(fixedDf).drop('id_len',axis=1)   

    else:
        joinedDf = annotatedDf.ix[annotatedDf.id_len==1,:].drop('id_len',axis=1)   
    
    try:
        ## get the original 17mer sequence, the mutated 17mer aa sequence and a list of 9mer slices
        joinedDf.ix[:,'Normal'] = joinedDf.apply(get_original,axis=1 )
        joinedDf.ix[:,'Mutated'] = joinedDf.apply(get_mutated,axis=1)
        joinedDf.ix[:,'refNonamer'] = joinedDf.Normal.apply(get_AA)
        joinedDf.ix[:,'mutNonamer'] = joinedDf.Mutated.apply(get_AA)
    except KeyError as e:
        print e 
        pass

    # key errors are caught using lists in the Normal / mutated columns 

    #combinedDf = joinedDf.merge(originalDf,on=['Chr','Start'])
    combinedDf = joinedDf.drop_duplicates(['Normal','Mutated'])
    
    ## need to handle peptide sequences that did not give at least a 9 mer after finding the ref / mut position
    combinedDf = combinedDf.ix[combinedDf.Normal.str.len()>=9,:]
    
    combinedDf.to_csv('./'+sampleName+'_predicted.txt',sep='\t',index=False)
    #return [combinedDf,originalDf, joinedDf]

    refBuffer= ''
    for row in combinedDf['refNonamer']:
        refBuffer+='\n'.join(row)
        refBuffer+='\n'

    with open('./'+sampleName+'_ref_antigen.txt','w') as refFile:
        refFile.write(refBuffer)
        refFile.close()


    mutBuffer= ''
    for row in combinedDf['mutNonamer']:
        mutBuffer+='\n'.join(row)
        mutBuffer+='\n'

    with open('./'+sampleName+'_mut_antigen.txt','w') as mutFile:
        mutFile.write(mutBuffer)
        mutFile.close()

    return combinedDf