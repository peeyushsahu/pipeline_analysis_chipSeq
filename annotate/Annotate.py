__author__ = 'peeyush'
import pandas as pd
import alignment.commons as paths
Path = paths.path()
basepath = Path.basepath


def geneid_converter(listofids, input_identifier = None , output_identifier = None):
    '''
    Method will take gene identifier and convert into desired one. Identifiers availble:
    [u'accession', u'alias', u'biocarta', u'chr', u'end', u'ensemblgene', u'ensemblprotein', u'ensembltranscript',
    u'entrezgene', u'exons', u'flybase', u'generif', u'go', u'hgnc', u'homologene', u'hprd', u'humancyc', u'interpro',
    u'ipi', u'kegg', u'mgi', u'mim', u'mirbase', u'mousecyc', u'name', u'netpath', u'pdb', u'pfam', u'pharmgkb', u'pid',
    u'pir', u'prosite', u'ratmap', u'reactome', u'reagent', u'refseq', u'reporter', u'retired', u'rgd', u'smpdb',
    u'start', u'strand', u'summary', u'symbol', u'tair', u'taxid', u'type_of_gene', u'unigene', u'uniprot',
    u'wikipathways', u'wormbase', u'xenbase', u'yeastcyc', u'zfin']
    :param input_identifier: input identifier eg. entrezgene
    :param output_identifier: list of output identifier eg. ["symbol", "ensembl.gene"]
    :param listofids: list of ids to be mapped eg. ['1', '10', '10001']
    :return: DataFrame of mapped ids
    '''
    import mygene
    if input_identifier is None:
        input_identifier = "entrezgene"
    if output_identifier is None:
        output_identifier = ["symbol", "ensembl.gene"]
    mygene_object = mygene.MyGeneInfo()
    mapped_dataframe = mygene_object.querymany(listofids, scopes=input_identifier, fields=output_identifier, species="human", as_dataframe=True)
    return mapped_dataframe



def parse_gtf(path):
    from pandas import read_csv
    #gtf_file = read_csv('/ps/imt/genome/human/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf', sep='\t', header=0)
    colnames = ['chr', 'transcript_annotation', 'feature', 'start', 'stop', 'score', 'strand', 'frame', 'further_information']
    gtf_file = read_csv(path, sep='\t')
    gtf_file.columns = colnames
    gtf_file['chr'] = gtf_file['chr'].astype(str)
    return gtf_file


def vectors_4_chromosome():
    import numpy as np
    chr = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X','Y']
    position = [249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,141213431,135534747,135006516,
            133851895,115169878,107349540,102531392,90354753,81195210,78077248,59128983,63025520,48129895,51304566, 155270560, 59373566]
    chr_len = dict(zip(chr, position))
    chr_vector = {}
    for k, v in chr_len.iteritems():
        chr_vector[k] = np.zeros(v, dtype=np.uint32)
    return chr_vector

def exon_intron_vector(path):
    import re
    print "reading GTF file"
    gtf_file = parse_gtf(path)
    print "creating chromosome vector"
    chr_vector = vectors_4_chromosome()
    print "Annotating chr_vector"
    region = ['exon', 'transcript']
    for i in region:
        if i == 'exon':
            print 'exon'
            gene_name = ''
            index = 0
            stop = 0
            for k, v in gtf_file[gtf_file['region'] == i].iterrows():
                if len(str(v['chr'])) < 5 and v['chr'] != 'MT':
                    chr_vector[str(v['chr'])][v['start']: v['stop']] = 2
                    #print v['start'], v['stop']
                    name = re.split(';| |"', v['further_information']).index('gene_name')
                    new_gene_name = re.split(';| |"', v['further_information'])[name+2]
                    if index > 0 and new_gene_name == gene_name:
                    #print gene_name
                    #print index
                        chr_vector[str(v['chr'])][stop+1: v['start']-1] = 3
                        gene_name = new_gene_name
                        index += 1
                    stop = v['stop']
                    gene_name = new_gene_name
                    index += 1
        if i == 'transcript':
            print 'transcript'
            for k, v in gtf_file[gtf_file['region'] == i].iterrows():
                if len(str(v['chr'])) < 5 and v['chr'] != 'MT':
                    chr_vector[str(v['chr'])][v['start']] = 1
                    if str(v['chr'])=='X':
                        print k, 'Start', v['start']
    del(gtf_file)
    return chr_vector

def annotate_intronexon_junction(df, chr_vector):
    '''
    Calculate Intron-Exon juction will be called for the cookie-cut (+-250) peaks from summit
    :param df: Dataframe
    :param chr_vector: annotated chromosome vector
    :return:
    '''
    import numpy as np
    dataFrame = df
    dataFrame[['chr']] = dataFrame[['chr']].astype(str)
    dataFrame['exon_intron_junction'] = 0
    dataFrame['exon_intron_percent'] = 0
    for k, v in dataFrame.iterrows():
        if len(v['chr']) < 5:
            #print v
            #print type(v['chr']),type(v['start']),type(v['stop'])
            summit = v['start']+v['summit']
            unique_elements = set(chr_vector[v['chr']][summit-250:summit+250])
            unique_count = np.bincount(chr_vector[v['chr']][v['start']:v['stop']])
            zip_unique = dict(zip(unique_elements, unique_count))
            if unique_elements.issuperset([2,3]) and zip_unique[2]>0 and zip_unique[3]>0:
                if unique_elements.issuperset([0]) and zip_unique[0]<10 or unique_elements.issuperset([0]) == False:
                    if float(zip_unique[3])/float(zip_unique[2]) < 2 and float(zip_unique[3]) / float(zip_unique[2]) > 0.5:
                        print zip_unique
                        dataFrame.loc[k, 'exon_intron_junction'] = 1
                        dataFrame.loc[k, 'exon_intron_percent'] = float(zip_unique[3])/float(zip_unique[2])
    dataFrame.to_csv(basepath + '/exon_intron_junction_PRMT6_RA.csv', sep=",", encoding='utf-8', ignore_index=True, index=False)
    return dataFrame


def next5genes_annotator(dataframe, path):
    '''
    This method will take a datatframe and compute 5 nearest genes from the summit of the peak.
    :param dataframe:
    :param path: To GTF file
    :return:
    '''
    import sys
    import timeit
    print 'Process: Reading GTF file'
    start = timeit.default_timer()
    gtfFile = parse_gtf(path)
    #print(gtfFile.head())
    gtfFile = gtfFile[gtfFile['feature'] == 'transcript']
    gtfFile.index = range(0, len(gtfFile))
    #print(gtfFile.head())
    gtfFile = gtfFile.sort(['chr', 'start'], ascending=True)
    chr_group = gtfFile.groupby('chr')['start']
    dataframe = dataframe.sort(['chr'], ascending=True)
    nearGene = pd.DataFrame(columns=['chr', 'start', 'stop', 'Next transcript strand', 'Next transcript gene name', 'GenomicPosition TSS=1250 bp, upstream=5000 bp', 'ngene', 'distance2tss', 'nstrand'])
    count = 0
    print 'Process: Annotating peaks with genes'
    for k, v in dataframe.iterrows():
        sys.stdout.write("\rNumber of peaks annotated:%d" % count)
        sys.stdout.flush()
        count += 1
        List = chr_group.get_group(v['chr'])
        key = v['start'] + v['summit']
        #print(List.shape, key, min(List.index), max(List.index))
        loc = gtf_binary_search(List, key, min(List.index), max(List.index))
        #print(v)
        #print(loc)
        if loc != 'KEY OUT OF BOUND':
            outdf = next5genes(gtfFile, loc, v)
            nearGene = nearGene.append(outdf, ignore_index=True)

    del(gtfFile)
    stop = timeit.default_timer()
    print '\nTime elapsed:', stop-start,' sec'
    nearGene.to_csv(
            basepath + '/further_analysis/PRMT6_old+new_only_Enhancer_bound_nearest_6_genes.txt',
            sep="\t", ignore_index=True, header=True)
    return nearGene


def gtf_binary_search(List, key, imin, imax):
    '''
    Binary search function uses recursive method to search the position/index of peak-summit in the GTF database.
    :param List: DF of sorted GTF file only with 'start'.
    :param key: Position to be searched in the data
    :param imin: First index
    :param imax: Last index
    :return:
    '''
    if key < List[imin] or key > List[imax]:
        return 'KEY OUT OF BOUND'
    if imax < imin:
        return 'KEY NOT FOUND'
    elif (imax - imin) == 1 and key > List[imin] and key < List[imax]:
        return imin, imax
    else:
        imid = imin + ((imax - imin)/2)
        if List[imid] > key:
            #print('upper', List[imid], key, imin, imid)
            # key is in upper subset
            return gtf_binary_search(List, key, imin, imid)
        elif List[imid] < key:
            #print('lower', List[imid], key)
            # key is in lower subset
            return gtf_binary_search(List, key, imid, imax)
        elif List[imid] == key:
            #print('found', List[imid], key)
            # key is in lower subset
            return imid, imid+1

def next5genes(gtffile, position, row):
    '''
    This method actually search for five nearest genes in GTF file.
    :param gtffile:
    :param position:
    :return:
    '''
    import re
    geneList = []
    upstreamindex = position[0]
    downstreamindex = position[1]
    nearGene = pd.DataFrame(columns=['chr', 'start', 'stop', 'Next transcript strand', 'Next transcript gene name', 'GenomicPosition TSS=1250 bp, upstream=5000 bp', 'ngene', 'distance2tss', 'nstrand'])

    while len(geneList) < 6:
        #print geneList, upstreamindex, downstreamindex
        ## Traversing to upstream
        minus = re.split(';| |"', gtffile.iloc[upstreamindex]['further_information'])
        mgene_pos = minus.index('gene_name')
        mGene = minus[mgene_pos+2]
        if mGene not in geneList:
            #print mGene
            newrow = {'chr': row['chr'],
                    'start': int(row['start']),
                    'stop': int(row['stop']),
                    'Next transcript strand': row['Next transcript strand'],
                    'Next transcript gene name': row['Next transcript gene name'],
                    'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                    'ngene': mGene,
                    'nstrand': gtffile.iloc[upstreamindex]['strand'],
                    'distance2tss': ((row['start'] + row['summit']) - gtffile.iloc[upstreamindex]['stop'])*-1}
            nearGene = nearGene.append(pd.Series(newrow), ignore_index=True)
            #pos = str(gtffile.iloc[upstreamindex]['chr'])+':'+str(gtffile.iloc[upstreamindex]['start'])+':'+str(gtffile.iloc[upstreamindex]['stop'])
            geneList.append(mGene)
        upstreamindex -= 1

        ## Traversing to downstream
        plus = re.split(';| |"', gtffile.iloc[downstreamindex]['further_information'])
        pgene_pos = plus.index('gene_name')
        pGene = plus[pgene_pos+2]
        if pGene not in geneList:
            #print pGene
            newrow = {'chr': row['chr'],
                    'start': int(row['start']),
                    'stop': int(row['stop']),
                    'Next transcript strand': row['Next transcript strand'],
                    'Next transcript gene name': row['Next transcript gene name'],
                    'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                    'ngene': pGene,
                    'nstrand': gtffile.iloc[downstreamindex]['strand'],
                    'distance2tss': gtffile.iloc[downstreamindex]['start'] - (row['start'] + row['summit'])}

            nearGene = nearGene.append(pd.Series(newrow), ignore_index=True)
            geneList.append(pGene)
            #pos = str(gtffile.iloc[downstreamindex]['chr'])+':'+str(gtffile.iloc[downstreamindex]['start'])+':'+str(gtffile.iloc[downstreamindex]['stop'])
        downstreamindex += 1
    return nearGene















