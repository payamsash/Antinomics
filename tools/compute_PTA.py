from pathlib import Path
import pandas as pd
import numpy as np
from scipy.io import loadmat

# --- Settings ---
df = pd.read_excel("../data/master.xlsx")
base_dir = Path("/Users/payamsadeghishabestari/temp_folder/audiometry")
df["PTA_L"] = np.nan
df["PTA_R"] = np.nan
pta_freqs = [500, 1000, 2000]

for subject in df["subject_ID"].values:
    subj_dir = base_dir / subject   

    for hemi in ["L", "R"]:
        files = [
                    f for f in subj_dir.glob("*.mat")
                    if f.name.lower().startswith(f"{subject.lower()} {hemi.lower()}")
                ]
        if not files:
            raise ValueError(f"{subject} has missing audiometry files.")

        fname = files[0]
        data = loadmat(fname)
        freqs = data["betweenRuns"]["var1Sequence"][0][0][0]
        thrs  = data["betweenRuns"]["thresholds"][0][0][0]
        order = np.argsort(freqs)
        freqs_sorted = freqs[order]
        thrs_sorted  = thrs[order]
        mask = np.isin(freqs_sorted, pta_freqs)
        if subject == "typy":
            mask = np.isin(freqs_sorted, [1000, 2000])

        if not mask.any():
            raise ValueError(f"{subject} has missing frequencies.")

        pta = thrs_sorted[mask].mean()
        df.loc[df["subject_ID"] == subject, f"PTA_{hemi}"] = pta

df["PTA"] = 0.5 * (df["PTA_L"] + df["PTA_R"])
df = df[['subject_ID', 'group', 'age', 'sex', 'PTA']]
df.to_csv("../data/master.csv")


################# create desing, and contrast matrixes
df2 = df.copy()

# --- 1) Dummy-code groups (one-hot) ---
df2['group_T'] = (df2['group'] == 'T').astype(int)
df2['group_C'] = (df2['group'] == 'C').astype(int)
df2['age_m'] = df2['age'] - df2['age'].mean()
df2['PTA_m'] = df2['PTA'] - df2['PTA'].mean()

design_cols = ['group_C', 'group_T', 'age_m', 'sex', 'PTA_m']
design = df2[design_cols]
contrast_C_T = np.array([1, -1, 0, 0, 0])
contrast_T_C = np.array([-1, 1, 0, 0, 0])

df2["subject_ID"] = df2["subject_ID"].astype(str) + ".mif"
subjects = df2['subject_ID']

subjects.to_csv('subjects.txt', index=False, header=False)
design.to_csv('design_matrix.txt', sep=' ', index=False, header=False, float_format='%.6f')
np.savetxt('contrast_C_T.txt', contrast_C_T.reshape(1, -1), fmt='%d', delimiter=' ')
np.savetxt('contrast_T_C.txt', contrast_T_C.reshape(1, -1), fmt='%d', delimiter=' ')