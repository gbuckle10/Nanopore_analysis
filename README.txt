Sample data can be gotten from s3://ont-open-data/. There might be a better way, but you can browse the data using the command:
    - aws s3 ls {folder_path}
    - This will show the contents of the folder you provided the path for.

First thing to do is run script.sh to download the fast_5 data from the Nanopore sample data.
    I think script.sh needs to be run using git bash


Basecalling - Use Dorado to do the basecalling
    Install Dorado (install it into the folder for this project, but if I put this on git it should be in .gitignore):
    NOTE ----- I HAD TO DOWNLOAD AN OLD VERSION OF DORADO TO DO THE BASECALLING FOR THE SAMPLE DATASET

You need to install awscli.
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