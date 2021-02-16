import argparse
import hail as hl
import os

"""
Simple script to read VCF(s) file(s) and write as Hail MatrixTable.
VCFs requirements: version 4.2, multi sample-called and BGZ compressed.

Usage example:

nohup python vcf2mt.py \
      --vcf_path "file:///home/ubuntu/data/test_vcfs/*.vcf.gz" \
      --output_path "file:///home/ubuntu/data/test.mt" \
      --force_bgz \
      --overwrite \
      > nohup.log 2>&1 &
"""


def main(args):

    # Start Hail on local mode
    hl.init(default_reference='GRCh38')

    # getting list of VCF files from given path
    # vcf_files_list = get_files_names(args.vcf_path, ext='vcf.gz')

    # import VCF(s) as Hail MatrixTable
    mt = hl.import_vcf(path=args.vcf_path, 
                       force_bgz=args.force_bgz)
    
    if args.split_multi:
        mt = hl.split_multi_hts(mt)
    
    # write MatrixTable
    mt.write(output=args.output_path,
             overwrite=args.overwrite)

    # Stop Hail
    hl.stop()
    
    print("Finished!")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--vcf_path', help='VCF files path')
    parser.add_argument('--output_path', help='Output MatrixTable path')
    parser.add_argument('--force_bgz', help='Forces BGZ decoding for VCF files with .gz extension.',
                        action='store_true')
    parser.add_argument('--split_multi', help='Split multi-allelic variants. See Hail function split_multi_hts for details.',
                        action='store_true')
    parser.add_argument('--overwrite', help='Overwrite if MatrixTable exists already.',
                        action='store_true')

    args = parser.parse_args()

    main(args)
