import os
import subprocess
import sys

from utils.runner import load_config, get_project_root, run_external_command

project_root = get_project_root()
CONFIG_PATH = os.path.join(project_root, "config.yaml")
RUNTIME_CONFIG_PATH = os.path.join(project_root, "src", "runtime_config.sh")

def methylation_pileup(config):
    alignment_output_dir = config['paths']['alignment_output_dir']
    aligned_sorted_file = os.path.join(
        project_root,
        alignment_output_dir,
        config['paths']['aligned_bam_name']
    )
    methylation_dir = config['paths']['methylation_dir']
    methylation_bed_name = config['paths']['methylation_bed_name']
    output_bed = os.path.join(
        project_root,
        methylation_dir,
        methylation_bed_name
    )

    pileup_cmd = [
        'modkit', 'pileup',
        aligned_sorted_file,
        output_bed
    ]

def main():
    config = load_config()
    methylation_pileup(config)

if __name__ == '__main__':
    main()
