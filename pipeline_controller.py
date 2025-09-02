import subprocess
import argparse
import shlex
import shutil


def run_command(command, step_name=""):
    """
    This function runs an external command, prints it, streams it output in real time and checks for errors
    :param command: The bash command to use
    :param step_name: (optional) the name of the current step
    :return:
    """

    if step_name:
        print(f"\n--- Starting Step: {step_name} ---")

    use_shell = isinstance(command, str) and any(c in command for c in ">|&;")

    try:
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            shell=use_shell,
            text=True,
            bufsize=1,
            executable="/bin/bash" if use_shell else None
        )

        for line in iter(process.stdout.readline, ""):
            print(line, end="")

        process.wait()

        if process.returncode != 0:
            raise subprocess.CalledProcessError(process.returncode, process.args)

        if step_name:
            print(f"--- Finished Step: {step_name} ---")

    except subprocess.CalledProcessError as e:
        print(f"\n\n--- ERROR ---")
        print(f"Error during '{step_name}': Command failed with exit code {e.returncode}")
        print(f"Command was: {e.cmd}")
        exit(1)

    except FileNotFoundError:
        print(f"\n\n--- ERROR ---")
        print(f"Error during '{step_name}': Command '{command.split()[0]}' not found.")


def main(args):
    print(f"Starting pipeline with {args.threads} threads")
    # Check here if the input files exist and create the output directories if necessary.

    dorado_cmd = [
        "dorado", "basecaller",
        "--device", "cuda:0",
        "batchsize", str(args.batchsize),
        args.model,
        args.input_pod5,
        ">", args.output_bam
    ]

    run_command(" ".join(dorado_cmd), "Basecalling")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Nanopore Pipeline")
    args = parser.parse_args()

    bash_path = shutil.which("bash.exe")
    subprocess.run([bash_path, 'script.sh'], check=True)

    #main(args)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
