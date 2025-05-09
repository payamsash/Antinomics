import argparse
from eeg.eeg_processing import process

def main():
    parser = argparse.ArgumentParser(description=("""
    ************************************
    Processing of the raw eeg recordings
    The process could be fully or semi automatic based on user choice.



    Notes
    -----
    This script is mainly designed for Antinomics / TIDE projects, however could be used for other purposes.
    """
    ),
    formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument("subject_id", help="file path to the eeg recording.")
    parser.add_argument("subjects_dir", help="The subject name, if subject has MRI data as well, should be FreeSurfer subject name, then data from both modality can be analyzed at once.")
    parser.add_argument("paradigm", default="rest", help="Name of the EEG paradigm, must be one of the: ['rest', 'rest_v1', 'rest_v2', 'gpias', 'xxxxx', 'xxxxy', 'omi', 'regularity']")
    parser.add_argument("--config_file", default=None, help="The yaml file containing site specific options.")
    parser.add_argument("--overwrite", default="warn", help="must be one of the ['ignore', 'warn', 'raise']")
    parser.add_argument("--verbose", default="ERROR", help="Control verbosity of the logging output. If None, use the default verbosity level.")

    args = parser.parse_args()
    process(**vars(args))

if __name__ == "__main__":
    main()