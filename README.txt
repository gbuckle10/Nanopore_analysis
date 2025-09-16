
Before running the programme, activate the conda environment:
    conda activate nanopore_analysis
Sometimes, before running conda commands in a new WSL terminal you'll have to run:
    eval "$(path/to/conda shell.bash hook)"

Sample data can be gotten from s3://ont-open-data/. There might be a better way, but you can browse the data using the command:
    - aws s3 ls {folder_path} --no-sign-request
    - This will show the contents of the folder you provided the path for.


When you're downloading a reference genome for alignment, you will need to download a reference genome with USCS-style
headers. e.g. chr1.

Download a full aligned "plate" using:
    aws s3 cp "s3://ont-open-data/gm24385_mod_2021.09/extra_analysis/alignment/20210510_1600_X1_FAQ32172_f02f2d1c.bam" "data/alignment_output/taken_alignment.bam" --no-sign-request

You need to install awscliv2:
    - Download the zip file using command: curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"
    - Unzip using command: unzip awscliv2.zip
    - Install using command: sudo ./aws/install
This should be run on Linux/Ubuntu. On colab it's fine for testing purposes, but if you're running this on a Windows machine
Setup for WSL:
    Installation:
        run wsl --install in powershell
        Restart and open WSL. Create a username and password.
        Download and install miniconda. In the command line run:
         - wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
         - bash Miniconda3-latest-Linux-x86_64.sh
         - Accept the licence agreement. Say no when it asks you whether you want to always open in the conda environment
         - If you accidentally said yes to automatically initialise conda, you can undo that (after restarting WSL) with the command: conda config --set auto_activate_base false
        Close WSL to finish installation
    Setup the project to run in WSL:
        Open WSL again and go to the project folder
        Check that the environment.yml file is there, then run conda env create -f environment.yml to set up the Conda environment.
    If you have changed the environment.yml file or want to update the environment for some other reason, use the following command:
        - conda env update --file environment.yml --prune
    It should be now possible to run the pipeline by running main.py. There are some possible issues though:
        Errors like $'\r': command not found:
            I edited the bash code in Windows, but it's running from Ubuntu which means that the line endings are different.
            If you're lucky you can change the line endings in your code editor. For example, in Pycharm:
                On the bottom right status bar there is a button which has CRLF/LF/CR. Change this to LF and it'll have the right file endings.
            If not, you can do it manually:
                In the WSL terminal, install dos2unix so that you can convert from Windows to Unix-style LF line endings. Run the following in the terminal:
                    sudo apt-get update (connects to configured package repositories and gets the latest package lists)
                    sudo apt-get install dos2unix (installs dos2unix)
                When it's been installed, convert every bash file by running the following:
                    find . -type f -name "*.sh" -exec dos2unix {} \+;
                        find . - find something
                        -type f - find files
                        -name "*.sh" - find files which end with ".sh". These are bash files.
                        -exec dos2unix - execute the dos2unix command
                        {} - this is replaced with the list of filenames found
                        \+ - tells find to add all of the filenames
                        ; - end of the line.

Useful base calling models:
    -