#! /usr/bin/env python

from run import *
import argparse
from version import __version__

def create_parser():

    parser = argparse.ArgumentParser(prog = "trans_savnet")

    parser.add_argument("--version", action = "version", version = "%(prog)s " + __version__)

    parser.add_argument("sample_list_file", metavar = "sample_list.txt", default = None, type = str,
                        help = "tab-delimited file of the cohort with the followin columns (the order can be arbitrary): \
                                Sample_Name: Sample labels, \
                                Mutation_Info: IDs for somatic mutations (e.g., SF3B1:R625H), \
                                SJ_File: Path to splicing junction data (generated by STAR), \
                                Weight: the path to intron retention count file (generated by intron_retention_utils simple_count)")

    parser.add_argument("output_prefix", metavar = "output_prefix", default = None, type = str, 
                        help = "the prefix of the output")

    parser.add_argument("--genome_id", choices = ["hg19", "hg38", "mm10"], default = "hg19",
                        help = "the genome id used for selecting UCSC-GRC chromosome name corresponding files (default: %(default)s)")

    parser.add_argument("--SJ_pooled_control_file", default = None, type = str,
                        help = "the path to control data created by junc_utils merge_control (default: %(default)s)")

    # parser.add_argument("--IR_pooled_control_file", default = None, type = str,
    #                     help = "the path to control data created by intron_retention_utils merge_control (default: %(default)s)")

    parser.add_argument("--SJ_num_thres", type = int, default = 5,
                        help = "extract splicing junctions whose supporting numbers are equal or more than this value \
                        at least one sample in the cohort (default: %(default)s)")

    # parser.add_argument("--keep_annotated", default = False, action = 'store_true',
    #                     help = "do not remove annotated splicing junctions")

    # parser.add_argument("--IR_num_thres", type = int, default = 3,
    #                     help = "extract intron retentions whose supporting numbers are equal or more than this value \
    #                     and the ratio is equal or more than IR_ratio_thres at least one sample in the cohort (default: %(default)s)")

    # parser.add_argument("--IR_ratio_thres", type = int, default = 0.05,
    #                    help = "extract intron retentions whose ratios (Intron_Retention_Read_Count / Edge_Read_Count) \
    #                     is equal or more than this value and supporting numbers are equal or more than IR_num_thres \
    #                     at least one sample in the cohort (default: %(default)s)")

    parser.add_argument("--permutation_num", type = int, default = 10,
                        help = "the number of permutation for calculating false discovery rate")

    parser.add_argument("--alpha0", type = float, default = 1.0,
                        help = "the shape parameter of prior Gamma Distribution for inactive states")

    parser.add_argument("--beta0", type = float, default = 1.0,
                        help = "the rate parameter of prior Gamma Distribution for inactive states") 

    parser.add_argument("--alpha1", type = float, default = 1.0,
                        help = "the shape parameter of prior Gamma Distribution for active states") 

    parser.add_argument("--beta1", type = float, default = 0.01,
                        help = "the rate parameter of prior Gamma Distribution for active states") 

    parser.add_argument("--log_BF_thres", type = float, default = 3.0,
                        help = "the threshould of logaraithm of Bayes Factor (default: %(default)s)")

    parser.add_argument("--effect_size_thres", type = float, default = 3.0,
                        help = "the thresould of effect size estimator used for simple edge pruning (default: %(default)s")

    # parser.add_argument("--zero_filter_prob", type = float, default = 0.5,
    #                     help = "if the fraction of non-zero support reads sample is above this value \
    #                     among those without any somatic mutation, this node will be pruned out (default: %(defaults)s")

    parser.add_argument("--debug", default = False, action = 'store_true', help = "keep intermediate files")
    
    
    return parser
 

