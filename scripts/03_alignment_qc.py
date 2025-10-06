import os
import subprocess
import sys

from utils.runner import load_config, get_project_root, run_external_command

project_root = get_project_root()
CONFIG_PATH = os.path.join(project_root, "config.yaml")
RUNTIME_CONFIG_PATH = os.path.join(project_root, "scripts", "runtime_config.sh")


def alignment_qc(config):
    alignment_output_dir = config['paths']['alignment_output_dir']
    aligned_sorted_file = os.path.join(
        project_root,
        alignment_output_dir,
        config['paths']['aligned_bam_name']
    )
    qc_dir = config['paths']['qc_dir']
    flagstat_report = os.path.join(
        project_root,
        qc_dir,
        config['paths']['alignment_flagstat_name']
    )

    flagstat_cmd = [
        "samtools", "flagstat",
        aligned_sorted_file
    ]

    print(f"Running command: {' '.join(flagstat_cmd)}")

    try:
        result = subprocess.run(
            flagstat_cmd,
            check=True,
            capture_output=True,  # Capture stdout/stderr
            text=True  # Decodes output as strings
        )

        with open(flagstat_report, 'w') as f:
            f.write(result.stdout)

        print(f"Flagstat report saved to {flagstat_report}")

    except FileNotFoundError:
        print(f"CRITICAL: 'samtools' not found.", file=sys.stderr)
        sys.exit(1)
    except subprocess.CalledProcessError as e:
        print(f"CRITICAL: samtools flagstat failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(1)


def main():
    config = load_config()
    alignment_qc(config)


if __name__ == '__main__':
    main()
