import os
import sys
import argparse
import logging.config
import pysam
from spliceai.utils import Annotator, get_delta_scores
from pkg_resources import resource_filename

try:
    from sys.stdin import buffer as std_in
    from sys.stdout import buffer as std_out
except ImportError:
    from sys import stdin as std_in
    from sys import stdout as std_out


os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # FATAL

import yaml
logger_path resource_filename(__name__,'config_files/spliceai_logging_config.yml')
with open(logger_path, 'r') as config:
    logging.config.dictConfig(yaml.safe_load(config))

logger = logging.getLogger('cloaked_chatter')


def get_options():

    parser = argparse.ArgumentParser(description='Version: 1.3.1')
    parser.add_argument('-I', metavar='input', nargs='?', default=std_in,
                        help='path to the input VCF file, defaults to standard in')
    parser.add_argument('-O', metavar='output', nargs='?', default=std_out,
                        help='path to the output VCF file, defaults to standard out')
    parser.add_argument('-R', metavar='reference', required=True,
                        help='path to the reference genome fasta file')
    parser.add_argument('-A', metavar='annotation', required=True,
                        help='"grch37" (GENCODE V24lift37 canonical annotation file in '
                             'package), "grch38" (GENCODE V24 canonical annotation file in '
                             'package), or path to a similar custom gene annotation file')
    parser.add_argument('-D', metavar='distance', nargs='?', default=50,
                        type=int, choices=range(0, 5000),
                        help='maximum distance between the variant and gained/lost splice '
                             'site, defaults to 50')
    parser.add_argument('-M', metavar='mask', nargs='?', default=0,
                        type=int, choices=[0, 1],
                        help='mask scores representing annotated acceptor/donor gain and '
                             'unannotated acceptor/donor loss, defaults to 0')
    args = parser.parse_args()

    return args


def main():

    args = get_options()

    if None in [args.I, args.O, args.D, args.M]:
        logger.error('Usage: spliceai [-h] [-I [input]] [-O [output]] -R reference -A annotation '
                      '[-D [distance]] [-M [mask]]')
        exit()

    try:
        vcf = pysam.VariantFile(args.I)
    except (IOError, ValueError) as e:
        logger.error('{}'.format(e))
        exit()

    header = vcf.header
    header.add_line('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 variant '
                    'annotation. These include delta scores (DS) and delta positions (DP) for '
                    'acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). '
                    'Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">')

    try:
        output = pysam.VariantFile(args.O, mode='w', header=header)
    except (IOError, ValueError) as e:
        logger.error('{}'.format(e))
        exit()

    ann = Annotator(args.R, args.A)

    for record in vcf:
        scores = get_delta_scores(record, ann, args.D, args.M)
        if len(scores) > 0:
            record.info['SpliceAI'] = scores
        output.write(record)

    vcf.close()
    output.close()


if __name__ == '__main__':
    main()
