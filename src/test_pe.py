import SICER_MS
import HTSeq 

genome_data = SICER_MS.get_genome_data('../genomes/hg19')
window_size = 200 
fragment_size = 150
scaling_factor = 1000000.
iterator = HTSeq.BAM_Reader('../input/sample_small_sorted.sam')

#bam_file = 'rmdup.sam'

rc_rf, wcd_rf, nwa_rf, tr_rf = SICER_MS.get_window_counts_pe(iterator, genome_data, window_size, scaling_factor)

#SICER_MS.paired_to_single_read('../input/' + bam_file)
#paired_to_single_bed_file = bam_file[:-4] + '_single' + '.bed'
#bed_iterator = HTSeq.BED_Reader('../input/' + paired_to_single_bed_file)
#fragment_size = 1


#read_dict, window_dict, tr_og = SICER_MS.make_dict_of_reads_and_windows(bed_iterator, genome_data, fragment_size, window_size)

#wcd_og, nwa_og = SICER_MS.get_window_counts_and_normalize(window_dict, read_dict, genome_data,scaling_factor, tr_og, window_size)

wcd_rf['chr1'].sort()

#mismatch = []
#for i in range(len(wcd_og['chr1'])):
 #   start_og = wcd_og['chr1'][i][0]
  #  count_og = wcd_og['chr1'][i][1]
   # window_og = HTSeq.GenomicInterval('chr1', start_og, start_og + window_size )
    #
   # start_rf = wcd_rf['chr1'][i][0]
#    count_rf = wcd_rf['chr1'][i][1]
 #   window_rf = HTSeq.GenomicInterval('chr1', start_rf, start_rf + window_size )
  #  if (start_og != start_rf) or (count_og != count_rf):
    
   #     print list(nwa_og[window_og].steps())
    #    print list(nwa_rf[window_rf].steps())
     #   print window_og
      #  print window_rf
       # mismatch.append(i)
      #  print 'i=' + str(i)
       # print ''


