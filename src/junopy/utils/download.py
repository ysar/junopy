# junopy - A Python API for the Juno spacecraft's particle and fields datasets.
# Copyright (C) 2018  Yash Sarkango

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""
Contains methods to download files.  This is a more lower-level method that
what most users need.  If you need to download data files, using the
appropriate `IndexFile` class is much better.
"""

import urllib.error
import urllib.request
import uuid
from pathlib import Path
from typing import Optional

from junopy import config

DEFAULT_CHUNK_SIZE = 32 * 1024 * 1024


def download_file(
    url: str,
    output_file: Path,
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    tmp_dir: Optional[Path] = None,
    show_progress: bool = True,
) -> bool:
    """
    Downloads a file from a URL and saves it to disk.

    Args:
        url (str): The URL of the file to download.
        output_file (Path): The destination file path.
        chunk_size (int, optional): Download chunk size in bytes (default: 16MB).
        tmp_dir (Path, optional): Temporary directory for storing files.
            If None, uses the default defined in `JunoPy.config`
        show_progress (bool, optional): If True, prints download progress.

    Returns:
        bool: True if download succeeded, False otherwise.
    """
    if tmp_dir is None:
        tmp_dir = config.ConfigVars.get("TMP_DIR")

    assert tmp_dir is not None

    tmp_dir.mkdir(parents=True, exist_ok=True)

    tmp_file = tmp_dir / str(uuid.uuid4())
    downloaded_bytes = 0

    try:
        # Open URL connection
        with urllib.request.urlopen(url) as response:
            total_bytes = int(response.headers.get("Content-Length", 0))

            # Download in chunks
            with open(tmp_file, "wb") as f:
                while True:
                    chunk = response.read(chunk_size)
                    if not chunk:
                        break
                    f.write(chunk)
                    downloaded_bytes += len(chunk)

                    if show_progress and total_bytes > 0:
                        progress = (downloaded_bytes / total_bytes) * 100
                        print(f"Downloading... {progress:.1f}%", end="\r")

            if total_bytes != downloaded_bytes:
                raise Exception(
                    "Download incomplete. "
                    f"{downloaded_bytes} of {total_bytes} bytes downloaded. "
                    "This error is usually a result of a PDS server limit. "
                    " Try downloading using the 'zip' option."
                )

        # Move temp file to final destination
        tmp_file.replace(output_file)
        if show_progress:
            print(f"\nDownloaded: {output_file}")
        return True

    except urllib.error.URLError as e:
        print(f"\nFailed to download {url}: {e}")
        if tmp_file.exists():
            tmp_file.unlink()  # Clean up temp file
        raise

    except Exception as e:
        print(f"\nUnexpected error downloading {url}: {e}")
        if tmp_file.exists():
            tmp_file.unlink()  # Clean up temp file
        raise


def download_multiple_files(
    urls: list[str],
    output_files: list[Path],
    chunk_size: int = DEFAULT_CHUNK_SIZE,
    num_threads: int = 4,
) -> bool:
    """
    Use threading to download multiple files in parallel.

    Args:
        urls (list[str]): List of URLs of the files to be downloaded
        output_files (list[Path]): List of output file paths for the
            downloaded files
        chunk_size (int, optional): Chunk size for downloads (default=16 MB)
        num_threads (int, optional): Number of threads to be used (default=4)

    Returns:
        bool: True if all downloads were successful.
    """
    from threading import Thread

    assert len(urls) > 0, "URL list cannot be empty"
    assert len(urls) == len(output_files), (
        "URLs and output file name lists must have the same length"
    )

    def download_wrapper(url, fOut):
        try:
            download_file(url, fOut, chunk_size=chunk_size, show_progress=False)
        except Exception as e:
            print(f"Failed to download {url}: {e}")

    total_files = len(urls)
    for i in range(0, total_files, num_threads):
        threads = []
        chunk_urls = urls[i : i + num_threads]
        chunk_files = output_files[i : i + num_threads]

        for url, fOut in zip(chunk_urls, chunk_files):
            thread = Thread(target=download_wrapper, args=(url, fOut))
            thread.start()
            threads.append(thread)

        for thread in threads:
            thread.join(timeout=300)

        progress = min((i + num_threads) / total_files * 100, 100)

        print(
            f"Parallel download (num_threads={num_threads}) = {progress:.1f}% "
            f"({i + len(chunk_urls)} / {total_files} files)",
            end="\r",
            flush=True,
        )

    print("\n ... Download complete.")

    # Check if all files are present.
    success = True
    for p in output_files:
        if not p.exists():
            print(f"File missing - {p}")
            success = False

    return success
