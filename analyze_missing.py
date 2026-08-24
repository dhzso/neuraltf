import pandas as pd
import numpy as np

df = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\runs\pipeline_run\rank_neural.csv')
STREAMS = ['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity']
W_DEFAULT = np.array([0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105])

# Show all candidates with their stream values and missingness
for _, row in df.iterrows():
    streams = [row[s] for s in STREAMS]
    n_missing = sum(pd.isna(v) for v in streams)
    n_avail = 7 - n_missing
    print(f'{row["gene_id"]} | {row["gene_name"]} | avail={n_avail}/7 | integrated={row["integrated_score"]:.6f} | streams={streams}')

print("\n\n=== Detailed analysis for dd14115 ===")
row = df[df['gene_id'] == 'dd_Smed_v6_14115_0_1'].iloc[0]
streams_vals = [row[s] for s in ['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity']]
print(f"Streams: {streams_vals}")
print(f"Gene: {row['gene_name']}")
print(f"Integrated score (from CSV): {row['integrated_score']:.6f}")

# Fixed-weight calculation
S = np.array([v if not pd.isna(v) else np.nan for v in row[['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity']].values])
mask = ~np.isnan(S)
fixed_score = np.sum(S[mask] * W_DEFAULT[mask]) / np.sum(W_DEFAULT[mask])
print(f"\nFixed-weight (renormalized): {fixed_score:.6f}")

# Verify it matches
print(f"CSV integrated_score: {row['integrated_score']:.6f}")
print(f"Match: {abs(fixed_score - row['integrated_score']) < 1e-6}")

# Now Dirichlet - one draw
rng = np.random.default_rng(2024)
k = 40
alpha = W_DEFAULT * 40 + 1e-9
w = np.random.default_rng(2024).gamma(alpha, 1.0)
w = w / w.sum()
print(f"\nSampled weights: {w}")
print(f"Sum: {w.sum():.6f}")

S_arr = np.array([v if not pd.isna(v) else np.nan for v in row[['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity']].values])
mask = ~np.isnan(S)
S_zero = np.where(np.isnan(S), 0, S)

# Current Dirichlet (zero imputation)
dirichlet_score = S_zero @ w
print(f"\nCurrent Dirichlet (zero imputation): {dirichlet_score:.6f}")

# Corrected Dirichlet (renormalized)
w_avail = w[mask] / w[mask].sum()
corrected_score = S[mask] @ w_avail
print(f"Corrected Dirichlet (renormalized): {corrected_score:.6f}")

# Show available streams
print("\nAvailable streams:")
for i, (s, v) in enumerate(zip(['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity'], row[['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity']].values)):
    status = "AVAIL" if not pd.isna(v) else "MISSING"
    print(f"  {s}: {v} ({status})")

# Weight sum for available streams
print(f"\nSum of all weights: {w.sum():.6f}")
print(f"Sum of available weights: {w[mask].sum():.6f}")
print(f"Missing weight mass: {w[~mask].sum():.6f}")