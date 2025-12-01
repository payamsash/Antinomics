import os
from pathlib import Path
from tqdm import tqdm
import numpy as np
import pandas as pd
from scipy.io import loadmat

data_dir = Path("/Volumes/G_USZ_ORL$/Research/ANTINOMICS/data")
eeg_dir = data_dir / "eeg"
audio_dir = data_dir / "audiometry"
files = [f for f in eeg_dir.iterdir() if f.is_file() and f.name.endswith("_rest.vhdr")]
files = sorted(files, key=os.path.getctime)
subject_tide_ids = [str(i) for i in range(70001, 70001 + len(files))]
subject_antinomics_ids = [file.stem[:4] for file in files]
subject_to_tide_map = dict(zip(subject_antinomics_ids, subject_tide_ids))

# The frequencies you want as columns
target_freqs = [125, 250, 500, 1000, 2000, 4000, 6000, 8000, 12000]

freq_cols = []
for hemi in ["L", "R"]:
    for f in target_freqs:
        freq_cols.append(f"{hemi}_{f} Hz")

rows = []
for subject in tqdm(subject_antinomics_ids):
    
    subject_audio_dir = audio_dir / subject
    tide_id = subject_to_tide_map.get(subject, np.nan)
    row = {"antinomics_id": subject, "tide_id": tide_id}

    # Initialize all possible columns with NaN so missing ones remain NaN
    for col in freq_cols:
        row[col] = np.nan

    if subject_audio_dir.exists():
        for fname in subject_audio_dir.iterdir():
            for hemi in ["L", "R"]:
                if fname.name.lower().startswith(f"{subject.lower()} {hemi.lower()}"):

                    data = loadmat(fname)
                    freqs = data["betweenRuns"]["var1Sequence"][0][0][0]
                    thrs  = data["betweenRuns"]["thresholds"][0][0][0]

                    # Sort
                    order = np.argsort(freqs)
                    freqs_sorted = freqs[order]
                    thrs_sorted  = thrs[order]

                    for f, thr in zip(freqs_sorted, thrs_sorted):
                        f_int = int(f)
                        col_name = f"A{hemi}E_{f_int}"
                        if col_name in row:
                            row[col_name] = thr

    rows.append(row)

df = pd.DataFrame(rows)
df = df[["antinomics_id", "tide_id"] + freq_cols]
df.to_csv("here.csv")