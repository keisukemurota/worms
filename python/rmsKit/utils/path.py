"""Utilities for path manipulation."""
import re
from pathlib import Path
import logging
from typing import Dict, List, Tuple, Union, cast
import pandas as pd
import numpy as np
from .functions import get_loss

logger = logging.getLogger(__name__)


def extract_to_rmskit(path: Path) -> Path:
    """
    Extracts the path up to and including the 'rmsKit' directory from the given path.

    Traverses the parents of the given path to find a directory named 'rmsKit' and returns the path up to that directory.

    Args:
        path (Path): The path from which to extract the 'rmsKit' directory.

    Returns:
        Path: The path up to and including the 'rmsKit' directory.

    Raises:
        ValueError: If the 'rmsKit' directory is not found in the path's ancestors.
    """
    for parent in path.parents:
        if parent.name == 'rmsKit':
            return parent
    raise ValueError("The given path is not under 'rmsKit'.")

def find_info_txt_files(directory_path: Union[str, Path]) -> List[Path]:
    """
    Finds all 'info.txt' files within the given directory path.

    Recursively searches the given directory for files named 'info.txt' and returns a list of Paths to these files.

    Args:
        directory_path (Union[str, Path]): The directory path in which to search for 'info.txt' files.

    Returns:
        List[Path]: A list of Paths to the found 'info.txt' files.

    Raises:
        ValueError: If the given directory_path is not a directory.
    """
    dir_path = Path(directory_path)
    if not dir_path.is_dir():
        raise ValueError(f"The given path {directory_path} is not a directory.")
    info_txt_files = list(dir_path.rglob('info.txt'))
    return info_txt_files

def extract_info_from_txt(file_path: Path) -> Dict[str, Union[float, Path]]:
    """Extract information from the given info.txt file."""
    patterns = {
        "best_loss": r"best loss: ([\d.e-]+)",
        "initial_loss": r"initial loss: ([\d.e-]+)",
        "hamiltonian_path": r"hamiltonian was saved to ([/\w._-]+)",
        "unitary_path": r"best loss was saved to ([/\w._-]+)"
    }
    res = {}

    # Read the file and extract information
    content = file_path.read_text()

    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if match:
            if "path" in key:
                res[key] = Path(match.group(1))
            else:
                res[key] = float(match.group(1))

    if None in res.values():
        missing_keys = [key for key, value in res.items() if value is None]
        raise ValueError(f"""
        Could not extract all required information from the given file. Missing: {', '.join(missing_keys)}
        The file content is:
        {content}
        """)

    return res

def top_k_upath(search_path: Path, k: int = 10) -> List[tuple[float, Path]]:

    """
    Get the top k unitary paths and their corresponding loss values.
    """


    if not search_path.exists():
        raise ValueError(f"The given search path {search_path} does not exist.")
    if not search_path.is_dir():
        raise ValueError(f"The given path {search_path} is not a directory.")
    
    # get all the paths to files that contain "u"
    unitary_paths = list(search_path.rglob("u"))
    loss_vals = [float(re.search(r"[\d.]+", path.parent.name).group()) for path in unitary_paths if re.search(r"[\d.]+", path.parent.name)]

    #argsort loss_vals
    sorted_indices = np.argsort(loss_vals)
    return [(loss_vals[i], unitary_paths[i]) for i in sorted_indices[:k]]


def get_worm_path(
        search_path: Path,
        ) -> Tuple[float, float, Path, Path, Path]:
    """
    Extracts the paths to the model hamiltonian and optimized unitary based on information found in 'info.txt' files within a given directory.

    This function searches for 'info.txt' files within the specified directory, extracts the paths to the model hamiltonian and optimized unitary,
    and returns these paths along with the best and initial loss values. Optionally, it can also return the path to the 'info.txt' file used.

    Args:
        search_path (Path): The directory path in which to search for 'info.txt' files.
        return_info_path (bool, optional): If True, the function also returns the path to the 'info.txt' file used. Defaults to False.

    Returns:
        Tuple[float, float, Path, Path, Path]: A tuple containing the best loss, initial loss, path to the optimized unitary,
        path to the model hamiltonian, and path to the 'info.txt' file.

    Raises:
        ValueError: If the search_path does not exist, is not a directory, or if no 'info.txt' file is found.
        ValueError: If the paths extracted from 'info.txt' are not valid Path objects or if the referenced files do not exist.
    """
    if not search_path.exists():
        raise ValueError(f"The given search path {search_path} does not exist.")
    if not search_path.is_dir():
        raise ValueError(f"The given path {search_path} is not a directory.")

    # Get info.txt files under the hamiltonian path
    info_txt_files = find_info_txt_files(search_path)
    if not info_txt_files:
        raise ValueError(f"No info.txt file found under the search path: {search_path}")
    if len(info_txt_files) > 1:
        logger.warning(f"Multiple info.txt files found under the search path. {info_txt_files[0]} will be used.")
    info_txt_file = info_txt_files[0]
    extracted_info = extract_info_from_txt(info_txt_file)

    if not isinstance(extracted_info["unitary_path"], Path):
        raise ValueError(f"The unitary path in the info.txt file is not a Path object: {extracted_info['unitary_path']}")
    if not isinstance(extracted_info["hamiltonian_path"], Path):
        raise ValueError(f"The hamiltonian path in the info.txt file is not a Path object: {extracted_info['hamiltonian_path']}")

    # Ensure paths are Path objects
    unitary_path = cast(Path, extracted_info["unitary_path"])
    hamiltonian_path = cast(Path, extracted_info["hamiltonian_path"])

    if not unitary_path.exists():
        raise ValueError(f"The path to the optimized unitary {unitary_path} does not exist.")
    if not hamiltonian_path.exists():
        raise ValueError(f"The path to the model hamiltonian {hamiltonian_path} does not exist.")

    loss = cast(float, extracted_info["best_loss"])
    initial_loss = cast(float, extracted_info["initial_loss"])

    return loss, initial_loss, unitary_path, hamiltonian_path, info_txt_file

def find_summary_files(directory_path: Union[str, Path]) -> List[Dict[str, Path]]:
    dir_path = Path(directory_path).resolve()
    if not dir_path.is_dir():
        raise ValueError(f"The given path {directory_path} is not a directory.")
    sim_res_folder = [path.parent for path in dir_path.rglob('summary/')]
    res = []
    for sim_res in sim_res_folder:
        summary_folder = sim_res / "summary"
        info_file = sim_res / "info.txt"
        if not summary_folder.is_dir():
            raise ValueError(f"Expected {summary_folder} to be a directory.")
        if not info_file.exists():
            raise ValueError(f"Info file {info_file} does not exist.")
        for path in summary_folder.rglob('*.csv'):
            dic = {"summary": path, "info": info_file}
            res.append(dic)
    return res


def get_df_from_summary_files(summary_files: List[Dict[str, Path]], N: int) -> pd.DataFrame:
    """
    Compile a DataFrame from summary CSV files that match a specific sweep number.

    Args:
    - summary_files (List[Dict[str, Path]]): A list of dictionaries, each containing paths to a summary file and its corresponding info file.
    - N (int): The sweep number to filter summary files by.

    Returns:
    - pd.DataFrame: A DataFrame compiled from the filtered summary files, enriched with initial loss information.
    """
    dfs = []
    for res_dict in summary_files:
        sum_file = res_dict["summary"]
        info_file = res_dict["info"]
        info = extract_info_from_txt(info_file)
        if f"sweeps_{N}" in sum_file.resolve().as_posix():
            try:
                df = pd.read_csv(sum_file)
                if not df.u_path.isnull().values.any():
                    df = get_loss(df)
                    # if np.abs(df.loss.values - info["best_loss"]).mean() > 1E-7:
                    #     raise ValueError("best_loss is not consistent with summary data")
                df["init_loss"] = info["initial_loss"]
                dfs.append(df)
            except Exception as e:
                logging.warning(f"Could not concatenate or read {sum_file}")
                logging.warning(e)
    if dfs:
        df = pd.concat(dfs).reset_index(drop=True)
    else:
        df = pd.DataFrame()
    return df


def get_sim_result(directory_path: Path, N: int) -> pd.DataFrame:
    """Get the simulation result from the given directory."""
    summary_files = find_summary_files(directory_path)
    df = get_df_from_summary_files(summary_files, N)
    return df
