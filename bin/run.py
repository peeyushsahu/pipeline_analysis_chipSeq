
import alignment.lanes
import alignment.aligner

genome = alignment.aligner.human_GRCh37()


raw_lanes = [
    alignment.lanes.Lane( 'Sample_8C9', '/ps/imt/f/20150320/150303_C00113_0099_AHFHHGADXX/Sample_8C9_test/'),
    #alignment.lanes.Lane( 'Sample_18F3', '/ps/imt/f/20150320/150303_C00113_0099_AHFHHGADXX/Sample_18F3_0/'),
    #alignment.lanes.Lane( 'Sample_18F3_RA', '/ps/imt/f/20150320/150303_C00113_0099_AHFHHGADXX/Sample_18F3_1/')
]

for lanes in raw_lanes:
    lanes.join_multiple_fq()
    lanes.do_alignment(genome)