from pathlib import Path
import shutil

def create_subjects_dir(
                        subject_id,
                        raws_dir,
                        target_dir
                        ):
    """
    Creates a folder structure with the following hierarchy:
    
    Antinomics/
        ├── raws/
            ├── sMRI_T1/
            ├── fMRI/
            │   ├── session_1/
            │   ├── session_2/
            ├── dMRI/
            ├── sMRI_T2/
        ├── subjects_fs_dir/
    
    Parameters
    ----------
    subject_id : str
        Subject ID used in Antinomics or TIDE projects.
    raws_dir: str │ Path
        The root directory, where all raw MR recordings of the subject is stored.
    """

    subject_dir = Path(raws_dir) / f"{subject_id}_antinomics"
    target_dir = Path(target_dir)

    for subfolder in ["sMRI_T1", "sMRI_T2", "fMRI/s1", "fMRI/s2", "dMRI"]:
        (target_dir / subfolder).mkdir(parents=True, exist_ok=True)

    for fname in subject_dir.iterdir():
        name = fname.name
        if not name.startswith(subject_id[:2]):
            raise ValueError(f"{subject_id} files are not correctly mapped.")

        if name.endswith("_t1w_3d_tfe_t1.nii"):
            shutil.move(fname, target_dir / "sMRI_T1" / f"{subject_id}.nii")

        if name.endswith("_3dt2_tra.nii"):
            shutil.move(fname, target_dir / "sMRI_T2" / f"{subject_id}.nii")

        if name.endswith("_3_1_fmri.nii") or name.endswith("_5_1_fmri.nii"):        
            shutil.move(fname, target_dir / "fMRI" /  "s1" / f"{subject_id}.nii")      
                
        if name.endswith("_4_1_fmri.nii") or name.endswith("_6_1_fmri.nii"):        
            shutil.move(fname, target_dir / "fMRI" /  "s2" / f"{subject_id}.nii")  

        if name.endswith("_dti_32.par"):
            shutil.move(fname, target_dir / "dMRI" / f"{subject_id}.par")
        
        if name.endswith("_dti_32.rec"):
            shutil.move(fname, target_dir / "dMRI" / f"{subject_id}.rec")

if __name__ == "__main__":
    raws_dir = Path("/home/ubuntu/volume/Antinomics/mri")
    target_dir = Path("/home/ubuntu/volume/Antinomics/raws")

    for folder in raws_dir.iterdir():
        if folder.is_dir() and not folder.name.startswith("."):
            subject_id = folder.stem.replace("_antinomics", "")
            create_subjects_dir(subject_id, raws_dir, target_dir)
