# Copyright (c) 2018 Yash Sarkango
# Licensed under the MIT License
"""Module to work with INDEX.TAB files present in the PDS3 archive."""

import logging
import pandas as pd
from pathlib import Path
from datetime import datetime

from junopy import config
from junopy import utils
from junopy.constants import DATA_IDS, PDS_URL, EXTENSIONS, TIME_FORMATS


class IndexFile:
    """
    Operates on PDS repository index files to fetch and parse data file URLs.

    Args:
        dataset (str): Dataset identifier (e.g., 'mag', 'waves-survey').
        do_download (bool): If True, forces download of index file even if it exists.
    """

    def __init__(
        self,
        dataset: str = "",
        do_download: bool = False,
    ):
        if not dataset:
            raise ValueError("Dataset must be specified.")
        if dataset not in DATA_IDS:
            raise ValueError(
                f"Invalid dataset: {dataset}. Valid options: {list(DATA_IDS.keys())}"
            )

        self.dataset = dataset

        # Set up paths
        self.outfolder = config.ConfigVars.get("INDEX_DIR")
        self.outfolder.mkdir(parents=True, exist_ok=True)
        self.index_file = self.outfolder / f"index_{self.dataset}.csv"

        # Download index file if needed
        if do_download or not self.index_file.exists():
            self.download_index_file()

        # Parse the index file
        self.df: pd.DataFrame = self.parse_index_file()

    def download_index_file(self) -> None:
        """
        Downloads the index file from PDS.
        """
        logging.info(f"Downloading index file for {self.dataset.upper()}...")
        temp_file = self.outfolder / "INDEX.TAB"

        try:
            utils.download.download_file(
                url=f"{PDS_URL}/{DATA_IDS[self.dataset]}/INDEX/INDEX.TAB",
                output_file=temp_file,
                show_progress=False,
            )
            temp_file.rename(self.index_file)  # Atomic rename
            logging.info(f"Index file saved to {self.index_file}")
        except Exception as e:
            if temp_file.exists():
                temp_file.unlink()
            raise RuntimeError(f"Failed to download index file: {e}")

    def parse_index_file(self) -> pd.DataFrame:
        """
        Parses the index file into a DataFrame.

        Returns:
            pd.DataFrame: The parsed index file as a pandas DataFrame.
        """
        logging.info(f"Parsing index file for {self.dataset.upper()}...")
        columns = [
            "VOLUME_ID",
            "SID",
            "DATA_SET_ID",
            "PRODUCT_ID",
            "START_TIME",
            "STOP_TIME",
            "FILE_SPECIFICATION_NAME",
            "CR_DATE",
            "PRODUCT_LABEL_MD5CHECKSUM",
        ]

        try:
            df = pd.read_csv(
                self.index_file,
                names=columns,
                parse_dates=["START_TIME", "STOP_TIME"],
                date_format=TIME_FORMATS[self.dataset],
                skiprows=1,
            )

            df["START_TIME"] = pd.to_datetime(df["START_TIME"])
            df["STOP_TIME"] = pd.to_datetime(df["STOP_TIME"])

            return df
        except Exception as e:
            raise RuntimeError(f"Failed to parse index file: {e}")

    def get_files_between_dates(
        self,
        start: str | datetime,
        stop: str | datetime,
        sid: str = "",
        do_download: bool = False,
        get_label: bool = False,
        do_parallel: bool = False,
        chunk_size: int = utils.download.DEFAULT_CHUNK_SIZE,
    ) -> list[Path]:
        """
        Gets files within a date range, optionally downloading them.

        Args:
            start (datetime): Start time (datetime or ISO-format string).
            stop (datetime): End time (datetime or ISO-format string).
            sid (str): Optional SID filter (e.g., '60s' for magnetometer data).
            do_download (bool): If True, downloads missing files.
            get_label (bool): If True, also downloads .LBL files.
            do_parallel (bool): If True, uses parallel downloads.
            chunk_size (int): Download chunk size in bytes.

        Returns:
            list[Path]: List of Path objects for relevant files.  These files
                will be downloaded only if `do_download` is `True`.
        """
        # Convert string dates to datetime
        if isinstance(start, str):
            start = pd.to_datetime(start)
        if isinstance(stop, str):
            stop = pd.to_datetime(stop)

        # Filter by date range
        mask = (self.df["START_TIME"] <= stop) & (self.df["STOP_TIME"] >= start)
        if sid:
            mask &= self.df["SID"].str.contains(sid)

        if self.dataset == "waves-survey":
            for s in ["BINT", "BUNC", "EINT", "EUNC"]:
                mask &= ~(self.df["PRODUCT_ID"].str.contains(s))

        df_filtered = self.df.loc[
            mask,
            [
                "FILE_SPECIFICATION_NAME",
                "START_TIME",
                "STOP_TIME",
                "PRODUCT_LABEL_MD5CHECKSUM",
            ],
        ]
        assert isinstance(df_filtered, pd.DataFrame)

        # Remove duplicate entries by using the MD5 checksum
        df_filtered.drop_duplicates(subset="PRODUCT_LABEL_MD5CHECKSUM", inplace=True)

        # Prepare file lists
        data_dir = config.ConfigVars.get("DATA_DIR") / self.dataset
        file_urls = []
        file_paths = []

        extension = EXTENSIONS[self.dataset]
        data_id = DATA_IDS[self.dataset]

        for row in df_filtered.itertuples():
            # Build URL and output path
            file_spec = row.FILE_SPECIFICATION_NAME.strip()
            file_url = f"{PDS_URL}/{data_id}/{file_spec[:-4]}{extension}"
            file_path = (
                data_dir
                / row.START_TIME.strftime("%Y")
                / row.START_TIME.strftime("%j")
                / file_url.split("/")[-1]
            )

            file_urls.append(file_url)
            file_paths.append(file_path)

        if not do_download:
            return file_paths

        # Download missing files
        download_urls = []
        download_paths = []
        lbl_urls = []
        lbl_paths = []

        for url, path in zip(file_urls, file_paths):
            path.parent.mkdir(parents=True, exist_ok=True)

            # Skip if file exists in any allowed format
            if any(
                p.exists()
                for p in [
                    path,
                    path.with_suffix(".hdf5"),
                    path.with_suffix(path.suffix + ".gz"),
                ]
            ):
                continue

            download_urls.append(url)
            download_paths.append(path)

            # Handle LBL files if requested
            if get_label and not path.with_suffix(".LBL").exists():
                lbl_urls.append(url[:-4] + ".LBL")
                lbl_paths.append(path.with_suffix(".LBL"))

        # Combine all downloads
        all_urls = download_urls + lbl_urls
        all_paths = download_paths + lbl_paths

        if not all_urls:
            logging.info("All files already exist. Nothing to download.")
            return file_paths

        logging.info(f"Downloading {len(all_urls)} files...")
        if do_parallel:
            utils.download.download_multiple_files(
                urls=all_urls,
                output_files=all_paths,
                chunk_size=chunk_size,
            )
        else:
            for url, path in zip(all_urls, all_paths):
                logging.info(f"Downloading {url} to {path}")
                utils.download.download_file(url, path, chunk_size=chunk_size)

        return file_paths
