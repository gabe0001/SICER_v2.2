import SICER_MS
import HTSeq 

genome_data = SICER_MS.get_genome_data('../genomes/hg19')
window_size = 200 
fragment_size = 150
scaling_factor = 1000000.
iterator = HTSeq.BAM_Reader('../input/sample_small_sorted.bam')
score_island_file_name = 'test.scoreisland'
graph_file_name = 'test.graph'


genome_fraction = .74
gap_size = 600
e_value = 1000

rc_rf, wcd_rf, nwa_rf, tr_rf = SICER_MS.get_window_counts(iterator, genome_data, window_size, fragment_size, scaling_factor)

read_dict, window_dict, tr_og = SICER_MS.make_dict_of_reads_and_windows(iterator, genome_data, fragment_size, window_size)

wcd_og, nwa_og = SICER_MS.get_window_counts_and_normalize(window_dict, read_dict, genome_data,scaling_factor, tr_og, window_size)


island_array, islands_list = SICER_MS.find_islands(wcd_rf, tr_rf, gap_size, window_size, genome_data,
                                                genome_fraction, e_value, score_island_file_name,
                                                graph_file_name, 2)
                                                
island_list, tot_chip, tot_cont = SICER_MS.count_reads_in_islands(islands_list, read_dict, read_dict)

islands_list_ref, tot_chip_ref, tot_cont_ref = SICER_MS.count_reads_in_islands_ref(islands_list, window_size, rc_rf, rc_rf)

