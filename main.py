import subprocess
import shutil


def start():
    bash_path = shutil.which("bash.exe")
    subprocess.run([bash_path, 'script.sh'], check=True)


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    start()

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
