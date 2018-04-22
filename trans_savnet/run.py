#! /usr/bin/env python

import sys, os, re, subprocess
from collections import namedtuple
from logger import get_logger
logger = get_logger(__name__)

try:
    import cPickle as pickle
except:
    import pickle

from savnet import preprocess, analysis_network, utils
from savnet.network import *
from savnet.utils import is_tool

Link_info_mut = namedtuple("Link_info_mut", ("Mutation_Key", "Motif_Pos", "Mutation_Type", "Is_Canonical", "Splicing_Key", "Splicing_Class", "Is_Inframe"))


def trans_savnet_main(args):

    ##########
    # check if the executables exist
    is_tool("bedtools")
    is_tool("tabix")
    is_tool("bgzip")
    is_tool("junc_utils")
    
    ##########
    output_prefix_dir = os.path.dirname(args.output_prefix)
    if output_prefix_dir != "" and not os.path.exists(output_prefix_dir):
       os.makedirs(output_prefix_dir)

    sample_list = []
    weight_vector = []
    mutation_status = {}
    link2info = {}
    sj_file_list = []

    ##########
    # organize mutation information
    logger.info("Organizing mutation data.")
    mut2ind = {}
    tind = 0
    sample_ind = 0
    with open(args.sample_list_file, 'r') as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            sample_list.append(F[header2ind["Sample_Name"]])
            weight_vector.append(float(F[header2ind["Weight"]]))
            sj_file_list.append(F[header2ind["SJ_File"]])

            mutation_info = F[header2ind["Mutation_Info"]]
            if mutation_info != "None":
                if mutation_info not in mut2ind: 
                    mut2ind[mutation_info] = tind
                    mutation_status[tind] = []
                    tind = tind + 1

                mutation_status[mut2ind[mutation_info]].append(sample_ind)

            sample_ind = sample_ind + 1                

    ##########
    logger.info("Merging splicing junction data.")
    preprocess.merge_SJ2(sj_file_list, args.output_prefix + ".SJ_merged.txt", args.SJ_pooled_control_file, args.SJ_num_thres, True)

    logger.info("Adding annotation to splicing junction data.")
    annotate_commands = ["junc_utils", "annotate", args.output_prefix + ".SJ_merged.txt", args.output_prefix + ".SJ_merged.annot.txt"]
    subprocess.call(annotate_commands)
    ##########

    ##########
    logger.info("Creating pickles of splicing association network instances.")
    out_s = open(args.output_prefix + ".splicing_mutation.network.pickles", 'wb')
    with open(args.output_prefix + ".SJ_merged.annot.txt", 'r') as hin:
        header2ind = {}
        header = hin.readline().rstrip('\n').split('\t')
        for (i, cname) in enumerate(header):
            header2ind[cname] = i

        for line in hin:
            F = line.rstrip('\n').split('\t')
            if F[header2ind["Splicing_Class"]] not in ["Exon skipping", "Alternative 3'SS", "Alternative 5'SS",
                                                       "Intronic alternative 3'SS", "Intronic alternative 5'SS"]: continue

            genes = F[header2ind["Gene_1"]].split(';') + F[header2ind["Gene_2"]].split(';')
            genes = sorted(list(set(genes)))

            if "---" in genes: genes.remove("---")
            if len(genes) > 0:
                genes_nm = filter(lambda x: x.find("(NM_") > 0, genes)
                if len(genes_nm) > 0: genes = genes_nm

            if len(genes) > 0:
                genes_single = filter(lambda x: x.find("-") == -1, genes)
                if len(genes_single) > 0: genes = genes_single
 
            gene = genes[0]
            gene = re.sub(r"\(N[MR]_\d+\)", "", gene)

            splicing_key = F[header2ind["SJ_1"]] + ':' + F[header2ind["SJ_2"]] + '-' + F[header2ind["SJ_3"]]
            splicing_counts = [[int(x) for x in F[header2ind["SJ_4"]].split(',')]]

            link2info = {}
            for k, v in mut2ind.items():

                # extract samples with the mutation of the link in consideration
                mut_vector = [0] * len(sample_list)
                for sample_id in mutation_status[v]:
                    mut_vector[sample_id] = 1

                link2info[(v, 0)] = Link_info_mut(k, "---", "---", "---", splicing_key, F[header2ind["Splicing_Class"]], F[header2ind["Is_Inframe"]])

            if len(link2info) > 0:
                network = Network(gene, mutation_status, splicing_counts, link2info, sample_list, weight_vector)
                pickle.dump(network, out_s)


    out_s.close()


    ##########
    logger.info("Extracting splicing associated variants.")
    sav_list_target = analysis_network.extract_sav_list(args.output_prefix + ".splicing_mutation.network.pickles", 
                                                        args.effect_size_thres, 0.5, 0.9, args.log_BF_thres, 2, 
                                                        args.alpha0, args.beta0, args.alpha1, args.beta1, permutation = False)

    ##########
    logger.info("Extracting of splicing associated variants on permutation pairs to estimate false positive ratios.")
    sav_lists_permutation = []
    for i in range(args.permutation_num):
        print >> sys.stderr, "Permutation: " + str(i)
        temp_sav_list = analysis_network.extract_sav_list(args.output_prefix + ".splicing_mutation.network.pickles", 
                                                          args.effect_size_thres, 0.5, 0.9, args.log_BF_thres, 2,
                                                          args.alpha0, args.beta0, args.alpha1, args.beta1, permutation = True)

        sav_lists_permutation.append(temp_sav_list)

    ##########
    logger.info("Adding Q-values to splicing associated variants.")
    analysis_network.add_qvalue_to_sav_list(sav_list_target, sav_lists_permutation)


    ##########
    logger.info("Generating final outputs.")
    with open(args.output_prefix + ".savnet.result.txt", 'w') as hout:
        print >> hout, Sav.print_header_mut
        for sav in sav_list_target:
            print >> hout, '\n'.join(sav.print_records(sv_mode = False, with_fdr = True))


    with open(args.output_prefix + ".splicing_mutation.count_summary.anno.perm_all.txt", 'w') as hout:
        for i in range(len(sav_lists_permutation)):
            for sav in sav_lists_permutation[i]:
                print >> hout, '\n'.join([str(i) + '\t' + x for x in sav.print_records(sv_mode = False, with_fdr = False)])


    if args.debug == False:

        subprocess.call(["rm", "-rf", args.output_prefix + ".SJ_merged.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".SJ_merged.annot.txt"])
        subprocess.call(["rm", "-rf", args.output_prefix + ".splicing_mutation.network.pickles"]) 

