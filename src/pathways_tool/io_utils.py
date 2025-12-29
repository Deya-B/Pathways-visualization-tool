import logging
import pandas as pd
import os

def read_csv(path, sep="\t", encodings=("utf-8", "utf-16", "cp1252")):
    """Read a delimited text file trying multiple encodings.

    The function checks that the file exists, then attempts to load it
    with each encoding in the given sequence until one succeeds. On
    success, it logs the chosen encoding; if all attempts fail, it logs
    a summary error listing the encodings tried and raises the last
    caught exception.

    Parameters
    ----------
    path : str
        Filesystem path to the TSV/CSV file.
    sep : str, optional
        Field delimiter passed to ``pandas.read_csv``, by default "\\t".
    encodings : tuple of str, optional
        Candidate text encodings to try in order, by default
        (\"utf-8\", \"utf-16\", \"cp1252\").

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the parsed contents of the file.

    Raises
    ------
    FileNotFoundError
        If the file does not exist at ``path``.
    UnicodeError
        If decoding fails for all encodings tried.
    Exception
        Re-raises the last error from ``pandas.read_csv`` if all
        attempts fail for reasons other than Unicode decoding.
    """
    if not os.path.isfile(path):
        logging.error(f"File not found: {path}")
        raise FileNotFoundError(path)

    last_err = None
    failed = []  # (encoding, error_message)

    for enc in encodings:
        try:
            df = pd.read_csv(path, sep=sep, encoding=enc)
            logging.info(f"Successfully read...\n  {path}\n  Encoding {enc}")
            return df
        except (UnicodeError, UnicodeDecodeError) as e:
            failed.append((enc, str(e)))
            last_err = e
        except Exception as e:
            failed.append((enc, str(e)))
            last_err = e

    # If none of them work, log error and raise last exception
    tried = ", ".join(enc for enc, _ in failed) or ", ".join(encodings)
    logging.error(
        f"Could not read {path} with any encoding (tried: {tried})"
    )
    # Debug lines for each failed encoding
    for enc, msg in failed:
        logging.debug(f"Encoding {enc} failed for {path}: {msg}")

    raise last_err
