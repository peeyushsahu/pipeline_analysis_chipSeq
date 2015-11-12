
import alignment.lanes
import alignment.aligner
import rna_seq.bam_processing as bamPrcessing

genome = alignment.aligner.human_GRCh37_74()


raw_lanes = [
    alignment.lanes.Lane('HL60_GFP3_1', '/ps/imt/f/christine/20151015/Sample_HL60_GFP3_1'),
    alignment.lanes.Lane('HL60_GFP3_2', '/ps/imt/f/christine/20151015/Sample_HL60_GFP3_2'),
    alignment.lanes.Lane('HL60_GFP3_3', '/ps/imt/f/christine/20151015/Sample_HL60_GFP3_3'),
    alignment.lanes.Lane('HL60_2_10_1', '/ps/imt/f/christine/20151015/Sample_HL60_Klon2_10_1'),
    alignment.lanes.Lane('HL60_2_10_2', '/ps/imt/f/christine/20151015/Sample_HL60_Klon2_10_2'),
    alignment.lanes.Lane('HL60_2_10_3', '/ps/imt/f/christine/20151015/Sample_HL60_Klon2_10_3')
]

## Aligner = Bowtie2, Tophat2, STAR
raw_lanes = dict((x.name, x) for x in raw_lanes)
print raw_lanes
alignedLane = {}
for lanes in raw_lanes:
    raw_lanes[lanes].join_multiple_fq()
    alignedLane[lanes] = raw_lanes[lanes].do_alignment(genome, 'Tophat2')
'''
## RNASeq diffcalling CuffDiff
controlSamples = ['HL60_GFP3_1', 'HL60_GFP3_2', 'HL60_GFP3_3']
conditionSamples = ['HL60_2_10_1', 'HL60_2_10_2','HL60_2_10_3']
bamPrcessing.cuffDiff(alignedLane, controlSamples, conditionSamples, genome)

## Count data for features using HTSeq-count
samples = ['HL60_GFP3_1', 'HL60_GFP3_2', 'HL60_GFP3_3', 'HL60_2_10_1', 'HL60_2_10_2', 'HL60_2_10_3']
bampaths = []
outpath = None
for name in samples:
    bampaths.append(alignedLane[name].bampath)
    if outpath is None: outpath = alignedLane[name].resultdir
## calling HTSeq-count
bamPrcessing.count_Data_featureCounts(bampaths, genome.gtfFile, stranded=1)
bamPrcessing.count_data(bampaths, genome.gtfFile, outpath)
'''

