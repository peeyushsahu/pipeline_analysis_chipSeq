__author__ = 'peeyush'
import pandas as pd
import alignment.commons as paths
import re, os
import numpy as np
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
    print('Parsing gtf file...')
    #gtf_file = read_csv('/ps/imt/genome/human/Homo_sapiens_Ensembl_GRCh37/Homo_sapiens/Ensembl/GRCh37/Annotation/Genes/genes.gtf', sep='\t', header=0)
    colnames = ['chr', 'transcript_annotation', 'feature', 'start', 'stop', 'score', 'strand', 'frame', 'further_information']
    gtf_file = read_csv(path, sep='\t')
    gtf_file.columns = colnames
    gtf_file['chr'] = gtf_file['chr'].astype(str)
    return gtf_file


def vectors_4_chromosome():
    import numpy as np
    chr_lengths = pd.read_csv('', header=None, sep='\t')
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


class AnnotateNearestGenes:
    '''
    This is will create an object for next gene annotation for peaks.
    '''
    def __init__(self, dataframe, path, maxdist=1000, maxgenes=2):
        self.dataframe = dataframe
        self.gtfpath = path
        self.maxdist = maxdist
        self.maxgenes = maxgenes

    def next_genes_annotator(self):
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
        gtfFile = parse_gtf(self.gtfpath)
        #print(gtfFile.head())
        gtfFile = gtfFile[gtfFile['feature'] == 'transcript']
        gtfFile.index = range(0, len(gtfFile))
        #print(gtfFile.head())
        gtfFile = gtfFile.sort(['chr', 'start'], ascending=True)
        chr_group = gtfFile.groupby('chr')['start']
        dataframe = self.dataframe.sort(['chr'], ascending=True)
        nearGene_df = pd.DataFrame(columns=['chr', 'start', 'stop', 'Next transcript strand', 'Next transcript gene name',
                                         'GenomicPosition TSS=1250 bp, upstream=5000 bp', 'gene_name', 'transcript_name', 'distance2tss', 'strand'])
        count = 0
        print 'Process: Annotating peaks with genes'
        for k, v in dataframe.iterrows():
            if len(v['chr']) < 4:
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
                    outdf = self.next_genes(gtfFile, loc, v)
                    nearGene_df = nearGene_df.append(outdf, ignore_index=True)

        del(gtfFile)
        stop = timeit.default_timer()
        print '\nTime elapsed:', stop-start,' sec'
        return nearGene_df

    def next_genes(self, gtffile, position, row):
        '''
        This method actually search for five nearest genes in GTF file.
        :param gtffile:
        :param position:
        :param maxdist: maximum distance from peak in kbs
        :return:
        '''
        import re
        maxdist = self.maxdist * 1000  # converting kb into basepairs
        #print(maxdist)
        geneList = []
        plusindex = position[0]
        minusindex = position[1]
        neargene = pd.DataFrame(columns=['chr', 'start', 'stop', 'Next transcript strand', 'Next transcript gene name', 'GenomicPosition TSS=1250 bp, upstream=5000 bp', 'gene_name', 'transcript_name', 'distance2tss', 'strand'])

        while len(geneList) < self.maxgenes:
            # print geneList, upstreamindex, downstreamindex
            dist_list = []
            for index in [plusindex, minusindex]:
                gtf_anno = re.split(';| |"', gtffile.iloc[index]['further_information'])
                strand = gtffile.iloc[index]['strand']
                gene_pos = gtf_anno.index('gene_name')
                gene_name = gtf_anno[gene_pos+2]
                tns_pos = gtf_anno.index('transcript_name')
                tns_name = gtf_anno[tns_pos+2]

                if strand == '+':
                    dist2tss = int((row['start'] + row['summit']) - gtffile.iloc[index]['start'])
                else:
                    dist2tss = int((row['start'] + row['summit']) - gtffile.iloc[index]['stop'])
                dist_list.append(abs(dist2tss))
                if (gene_name not in geneList) and (abs(dist2tss) < maxdist):
                    #print mGene
                    newrow = {'chr': row['chr'],
                            'start': int(row['start']),
                            'stop': int(row['stop']),
                            'Next transcript strand': row['Next transcript strand'],
                            'Next transcript gene name': row['Next transcript gene name'],
                            'GenomicPosition TSS=1250 bp, upstream=5000 bp':row['GenomicPosition TSS=1250 bp, upstream=5000 bp'],
                            'gene_name': gene_name,
                            'transcript_name': tns_name,
                            'strand': strand,
                            'distance2tss': int(dist2tss)*-1}
                    neargene = neargene.append(pd.Series(newrow), ignore_index=True)
                    geneList.append(gene_name)
                    if len(geneList) >= self.maxgenes:
                        break
            if (dist_list[0] > maxdist) and (dist_list[1] > maxdist):
                break
            plusindex += 1
            minusindex -= 1
        return neargene


class Generate_annotate_vectors:

    def __init__(self, gtf_path, tss_distance_up=1500, tss_distance_down=500, upstream=5000):
        self.gtf_path = gtf_path
        self.tssdistance_up = tss_distance_up
        self.tssdistance_down = tss_distance_down
        self.upstream = upstream
        self.gtf_file = parse_gtf(self.gtf_path)
        self.chr_lengths = self.get_chr_size()
        self.longest_transcript = self.get_longest_transcript()

    def get_chr_size(self):
        print('Reading chr sizes...')
        chr_size_path = os.path.join('/ps/imt/f/Genomes/geneAnotations', 'chrom_sizes')
        chr_lengths = pd.read_csv(chr_size_path, header=None, sep='\t')
        return chr_lengths

    def get_longest_transcript(self):
        print('Estimating longest transcript per gene...')
        longest_tr_path = os.path.join('/ps/imt/f/Genomes/geneAnotations', 'longest_transcript_annotation.db')
        longest_tr = pd.read_csv(longest_tr_path, header=0, sep='\t')
        return longest_tr

    def chromosome_vactors(self):
        '''
        create chromosome zero vector for region annotaion.
        '''
        chr_vector = {}
        #print(self.chr_lengths[0])
        for k, v in self.chr_lengths.iterrows():
            chr_vector[str(v[0])] = np.zeros(v[1]+1, dtype=np.uint32)
        return chr_vector

    def save_chr_vector(self, chr_vector):
        for chro, vec in chr_vector.items():
            f = open("/ps/imt/f/Genomes/geneAnotations/chr_vec/"+str(chro),  "wb")
            np.save(f, vec)
            f.close()

    def annotate_chrvectors(self):
        '''
        fill chromosome vector with respective region annotation
        '''
        chr_vector = self.chromosome_vactors()
        # we will annotate tss, exon, upstream
        gene_gtf = self.gtf_file
        gene_gtf = gene_gtf[(gene_gtf['chr'].str.len() < 4) | (gene_gtf['chr'].str.contains('GL'))]

        for feature, sym in zip(['exon', 'transcript'],[1, 3]):
            gtffile = gene_gtf[gene_gtf['feature'] == feature]
            gtffile = gtffile.sort_values(['chr', 'start'], ascending=True)
            #print(gtffile)
            gtf_chr_group = gtffile.groupby('chr')
            if feature == 'exon':
                for chr, df in gtf_chr_group:
                    vector = chr_vector[str(chr)]
                    for ind, row in df.iterrows():
                        start = row['start'] - 1  # converting into o based system
                        stop = row['stop']
                        vector[start: stop] = sym

            if feature == 'transcript':
                for chr, df in gtf_chr_group:
                    vector = chr_vector[str(chr)]
                    for ind, row in df.iterrows():
                        start = row['start'] - 1
                        stop = row['stop']
                        if row['strand'] == '-':
                            vector[stop-self.tssdistance_down: stop+self.tssdistance_up] = sym
                            tr_vec = np.unique(vector[stop+self.tssdistance_up: stop+self.upstream])
                            if len(tr_vec) == 1 and tr_vec[0] == 0:
                                vector[stop+self.tssdistance_up: stop+self.upstream] = 4
                        else:
                            vector[start-self.tssdistance_up: start+self.tssdistance_down] = sym
                            tr_vec = np.unique(vector[start-self.upstream: start-self.tssdistance_up])
                            if len(tr_vec) == 1 and tr_vec[0] == 0:
                                vector[start-self.upstream: start-self.tssdistance_up] = 4
        # now we will annotate intron
        longest_transcript_group = self.longest_transcript.groupby('chr')
        for key, df in longest_transcript_group:
            print(key, list(df['chr'])[0])
            vector = chr_vector[str(key)]
            df = df.sort_values(['chr', 'start'], ascending=[True, True])
            for ind, row in df.iterrows():
                tr_vec = vector[row['start']-1: row['stop']]
                tr_vec[tr_vec == 0] = 2
                vector[row['start']-1: row['stop']] = tr_vec
                #for pos in range(row['start']-1, row['stop']-1, 1):
                #    if vector[pos] == 0:
                #        vector[pos] = 2
        self.save_chr_vector(chr_vector) # save chr vector in bin files
        #return chr_vector


class AnnotatePeaks:

    def __init__(self, gtf_path, path_chr_vector='/ps/imt/f/Genomes/geneAnotations/chr_vec', dataframe=None):
        self.gtf_path = gtf_path
        self.dataframe = dataframe
        self.gtf_file = parse_gtf(self.gtf_path)
        self.path_chr_vector = path_chr_vector
        self.chr_vector = None

    def get_chr_vector(self):
        chr_vector = {}
        files = os.listdir(self.path_chr_vector)
        for chro in files:
            if os.path.isfile(os.path.join(self.path_chr_vector, chro)):
                chr_vector[chro] = np.load(os.path.join(self.path_chr_vector, chro))
        self.chr_vector = chr_vector

    def annotate_region(self):
        '''
        This will annotate region for the peaks (intron, exon)
        '''
        annotated_df = pd.DataFrame(columns=['chr', 'start', 'stop', 'summit', 'region', 'next_transcript',
                                             'next_transcript_tss_distance', 'strand'])
        gtffile = self.gtf_file
        gtffile = gtffile.sort(['chr', 'start'], ascending=True)
        gtf_chr_group = gtffile.groupby('chr')['start']
        peak_df = self.dataframe.sort(['chr'], ascending=True)

        # iterating peak data frame
        for ind, row in peak_df.iterrows():
            List = gtf_chr_group.get_group(row['chr'])
            key = row['start'] + row['summit']
            loc = gtf_binary_search(List, key, min(List.index), max(List.index))
            if loc != 'KEY OUT OF BOUND':
                outdf = self.genomic_region(gtffile, loc, key)

        return

    def genomic_region(self, loc, key):
        gtffile = self.gtf_file
        region = ''
        pos1 = loc[0]
        pos2 = loc[1]
        for pos in [pos1, pos2]:
            if (key < gtffile.iloc[pos]['stop']) & (key > gtffile.iloc[pos]['start']):
                region = re.split(';| |"', gtffile.iloc[pos]['further_information'])

        return region


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


def create_longestTranscriptDB(path):
    '''
    This will filter longest transcript per gene.
    :param path: path of transcript_annotation.db from def create_gtf2transcriptDB()
    :return:
    '''
    ## Parsing gtf file and storing gene names with position.
    gene_gtf = pd.read_csv(path+'/transcript_annotation.db', sep='\t', header=0, index_col=None)
    gene_gtf['chr'] = gene_gtf['chr'].astype(str)
    gene_gtf = gene_gtf[(gene_gtf['chr'].str.len() < 4) | (gene_gtf['chr'].str.contains('GL'))]
    #print(gene_gtf['chr'].value_counts())
    #gene_gtf = gene_gtf
    transcript_group = gene_gtf.groupby('gene_name')
    ## for intron select min(list(start) and max(list(stop)))
    with open(os.path.join(path, 'longest_transcript_annotation.db'), 'a') as annotationDB:
        annotationDB.write('chr'+'\t'+'start'+'\t'+'stop'+'\t'+'strand'+'\t'+'gene_name'+'\t'+'gene_id\n')
        for gene, gtf in transcript_group:
            gene_name = gene
            chr = list(gtf['chr'])[0]
            gene_id = list(gtf['gene_id'])[0]
            strand = list(gtf['strand'])[0]
            start = gtf['start'].min()
            stop = gtf['stop'].max()
            annotationDB.write(str(chr)+'\t'+str(start)+'\t'+str(stop)+'\t'+str(strand)+'\t'+str(gene_name)+'\t'+str(gene_id)+'\n')
    annotationDB.close()














