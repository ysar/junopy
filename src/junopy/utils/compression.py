import shutil
from pathlib import Path

from isal import igzip


def create_gz_zip_file(file_path: Path, remove_original: bool = True) -> None:
    """
    Compresses a file using igzip into a .gz file.

    Args:
        file_path(Path):  Path to file to be compressed.  This method will
            create a new file with the same name with the '.gz' extension.
        remove_original(bool):  If True, the original file will be deleted.
    """
    with open(file_path, "rb") as f_in:
        suffix = file_path.suffix
        with igzip.open(file_path.with_suffix(suffix + ".gz"), "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

    if remove_original:
        file_path.unlink()
