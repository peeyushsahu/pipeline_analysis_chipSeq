import multiprocessing
import os
import plotsAndseq.plots as HeatMap
import pandas as pd
from pandas import read_csv
import timeit
import alignment.commons as paths

__author__ = 'peeyush'


class Consumer(multiprocessing.Process):

    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        process_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # We will look for a POISON PILL :D
                print('%s:Exiting' %process_name)
                self.task_queue.task_done()
                break
            print('%s' %(process_name))
            result = next_task.computeheatmap()
            # print('place:',self.task_queue.task_done())
            self.task_queue.task_done()
            self.result_queue.put(result)
        return


class BamHeatmapFromPeak(object):

    def __init__(self, peak_df, bam_path, name, outdir):
        self.peak_df = peak_df
        self.bam_path = bam_path
        self.name = name
        self.outdir = outdir

    def computeheatmap(self):
        overlap_dict = {}
        print(self.name)
        try:
            normdf, df = HeatMap.overlapping_peaks_distribution(self.name, self.bam_path, self.peak_df, self.outdir)
        except:
            raise ValueError(self.name, 'FAILED!!!')
        overlap_dict[self.name] = normdf
        return overlap_dict


def join_result_dict_into_df(peak_df, df_dict):
    new_df = pd.DataFrame()
    sample_names = []
    for key, df in df_dict.items():
        new_df = pd.concat([new_df, df], axis=1)
        sample_names.append(key)

    new_df = pd.concat([peak_df[['chr', 'start', 'stop', 'Next transcript strand', 'summit']], new_df], axis=1)
    print(', '.join(sample_names))
    return new_df, sample_names


if __name__ == '__main__':

    start = timeit.default_timer()
    db_path = '/ps/imt/e/Encode_data_all/ENCODE_bam'
    out_dir = '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/ENCODE_heatmaps'
    paths.ensure_path(out_dir)

    peak_df = read_csv('/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/further_analysis/H3R2me2a_analysis/H3R2me2a_E9_RA vs IgG_E9_RA filtered/H3R2me2a_E9_RA vs IgG_E9_RA filtered.txt',header=0, sep='\t')
    peak_df = peak_df[peak_df['chr'].str.len() < 4]

    meta_df_bam = read_csv(os.path.join(db_path, 'metadata.tsv'), sep='\t', header=0)
    meta_df_bam['Experiment target'] = meta_df_bam['Experiment target'].map(lambda x: x.split('-')[0].strip())
    #print(meta_df_bam['Experiment target'])

    # Include samples from other source
    inhouse_sample = {
        'H3R2me2a_E9_RA': '/ps/imt/e/20141009_AG_Bauer_peeyush_re_analysis/results/AlignedLane/H3R2me2a_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup/aligned_unique_H3R2me2a_E9_RA__aligned_with_bowtie2_against_EnsemblGenome_Homo_sapiens_74_37_dedup.bam',
    }

    # enter desired cell lines and targets
    desired_cellline = ['H1-hESC', 'NT2/D1']
    desired_targets = ['EP300', 'H3K27ac']
    print('Size of ENCODE database:', len(meta_df_bam))
    if len(desired_cellline) > 0 and len(desired_targets) > 0:
        meta_df_bam = meta_df_bam[(meta_df_bam['Biosample term name'].isin(desired_cellline)) & (meta_df_bam['Experiment target'].isin(desired_targets))]
    elif len(desired_cellline) > 0 and len(desired_targets) == 0:
        meta_df_bam = meta_df_bam[(meta_df_bam['Biosample term name'].isin(desired_cellline))]
    if not len(meta_df_bam) > 0:
        raise ValueError('After filtering no encode sample left for analysis!!!')

    print('ENCODE files after filtering:', len(meta_df_bam))
    meta_df_bam.index = range(0, len(meta_df_bam))
    print(meta_df_bam.head())

    # writing information for analysed encode samples
    with open(os.path.join(out_dir, 'analysed_encode_sample_info.txt'), 'w') as File:
        File.write('File accession'+'\t'+'target'+'\t'+'Experiment accession'+'\t'+'target_type'+'\t'+'cell_line'+'\n')
        for ind, row in meta_df_bam.iterrows():
            sample_name = row['File accession'].strip()
            experiment_accession = row['Experiment accession'].strip()
            target = row['Experiment target'].split('-')[0].strip()
            cell_line = row['Biosample term name'].strip()
            target_type = row['Target type'].strip()
            file_format = row['File format'].strip()
            File.write(sample_name+'\t'+target+'\t'+experiment_accession+'\t'+target_type+'\t'+cell_line+'\n')
        File.close()

    # Establish communication queues
    task = multiprocessing.JoinableQueue()
    result = multiprocessing.Queue()

    # Start consumer
    new_consumer = int(multiprocessing.cpu_count()/2)
    print('Creating %d consumers' %new_consumer)

    consumers = [Consumer(task, result) for i in range(new_consumer)]

    for con in consumers:
        con.start()

    # Enqueue ENCODE jobs
    num_of_job = 4
    for ind, row in meta_df_bam.iterrows():
        sample_path = os.path.join(db_path, row['File accession']+'.bam')
        name = '_'.join([row['Experiment target'], row['Biosample term name'], row['File accession']])
        if os.path.exists(sample_path):
            task.put(BamHeatmapFromPeak(peak_df, sample_path, name, out_dir))
        else:
            print('Sample path does not exist:', sample_path)

    # Enqueue in house jobs
    for name, path in inhouse_sample.items():
        if os.path.exists(path):
            task.put(BamHeatmapFromPeak(peak_df, path, name, out_dir))
        else:
            print('Sample path does not exist:', path)

    # Add a poison pill
    for i in range(new_consumer):
        task.put(None)
    # join the task in the end
    task.join()

    # join the result
    result_dict = {}
    size = result.qsize()
    while size:
        results = result.get()
        #print(results)
        result_dict.update(results)
        size -= 1
    out_df, smaple_names = join_result_dict_into_df(peak_df, result_dict)
    stop = timeit.default_timer()
    print(out_df.head())
    out_df.to_csv(os.path.join(out_dir, 'Heatmap_results.txt'), sep='\t', header=True)
    print('Time consumed with '+str(new_consumer)+' processors:', ((stop-start)/60)/60, 'hrs')
