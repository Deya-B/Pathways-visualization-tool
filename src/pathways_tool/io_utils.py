import logging
import pandas as pd
import os

def read_csv(path, sep="\t", encodings=("utf-8", "utf-16", "cp1252")):
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
