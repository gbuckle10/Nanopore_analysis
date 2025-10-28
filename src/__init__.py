from pathlib import Path

SRC_INIT_PATH = Path(__file__).resolve()

SRC_DIR = SRC_INIT_PATH.parent

PROJECT_ROOT = SRC_DIR.parent

DATA_DIR = PROJECT_ROOT / "data"