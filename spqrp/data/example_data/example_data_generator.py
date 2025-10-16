import pandas as pd
import numpy as np

np.random.seed(42)

# Parameters
n_patients = 30
samples_per_patient = 2
proteins = ["ProtA", "ProtB", "ProtC", "ProtD", "ProtE"]
noise_level = 0.05  # small variation between samples of the same patient

data = []

for pid in range(1, n_patients + 1):
    patient_id = f"P{pid:03d}"

    # generate a baseline intensity vector for this patient
    patient_baseline = np.random.rand(len(proteins))
    for s in range(1, samples_per_patient + 1):
        sample_id = f"{pid}_S{s}"

        # add small noise for each sample
        sample_intensity = patient_baseline + np.random.normal(
            0, noise_level, len(proteins)
        )
        sample_intensity = np.clip(
            sample_intensity, 0, 1
        )  # keep between 0 and 1, can be adjusted to represent e.g. log_2 transformed data

        for prot, intensity in zip(proteins, sample_intensity):
            data.append([sample_id, patient_id, prot, np.round(intensity, 2)])

df_generated = pd.DataFrame(
    data, columns=["Sample_ID", "Patient_ID", "Protein", "Intensity"]
)

df_generated.to_csv("example_input_cohort_df.csv", index=False)

df_generated.head(10)
