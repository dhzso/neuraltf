import pandas as pd
import numpy as np

df = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\runs\pipeline_run\rank_neural.csv')
STREAMS = ['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity']
W_DEFAULT = np.array([0.211, 0.105, 0.158, 0.158, 0.105, 0.158, 0.105])

N_DRAWS = 1000
K_DIR = 40
SEED = 2024

# Compute for all candidates
results = []
for _, row in df.iterrows():
    gene_id = row['gene_id']
    gene_name = row['gene_name']
    
    # Get stream values
    streams = np.array([row[s] if not pd.isna(row[s]) else np.nan for s in ['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity']])
    mask = ~np.isnan([row[s] for s in ['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity']])
    n_avail = sum(mask)
    n_missing = 7 - sum(mask)
    
    # Fixed score (renormalized)
    S = np.array([v if not pd.isna(v) else np.nan for v in [row['expression'], row['specificity'], row['reproducibility'], row['rnai'], row['correlation'], row['neural_enriched'], row['neural_specificity']]])
    mask_arr = ~np.isnan(np.array([row['expression'], row['specificity'], row['reproducibility'], row['rnai'], row['correlation'], row['neural_enriched'], row['neural_specificity']]))
    fixed_score = np.sum(np.where(mask_arr, S, 0) * W_DEFAULT) / np.sum(W_DEFAULT[mask_arr]) if mask_arr.any() else 0
    
    # Dirichlet sensitivity - 1000 draws
    rng = np.random.default_rng(2024)
    S_arr = np.array([row[s] if not pd.isna(row[s]) else np.nan for s in ['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity']])
    mask_arr = ~np.isnan([row['expression'], row['specificity'], row['reproducibility'], row['rnai'], row['correlation'], row['neural_enriched'], row['neural_specificity']])
    
    # Fixed score
    fixed_score = np.sum(np.where(mask_arr, S, 0) * W_DEFAULT) / np.sum(W_DEFAULT[mask_arr]) if mask_arr.any() else 0
    
    # Dirichlet centered (1000 draws)
    rng = np.random.default_rng(2024)
    n_draws = 1000
    k = 40
    n = 1
    all_scores_dir = np.zeros(1000)
    all_scores_unif = np.zeros(1000)
    
    S_arr = np.array([row['expression'], row['specificity'], row['reproducibility'], row['rnai'], row['correlation'], row['neural_enriched'], row['neural_specificity']])
    S_arr = np.where(pd.isna(S_arr), np.nan, S_arr)
    mask = ~np.isnan(S_arr)
    
    # Fixed score (renormalized)
    fixed_score = np.sum(S[mask] * W_DEFAULT[mask]) / W_DEFAULT[mask].sum() if mask.any() else 0
    
    # Dirichlet centered
    all_scores_dir = np.zeros(1000)
    for d in range(1000):
        alpha = W_DEFAULT * 40 + 1e-9
        w = np.random.default_rng(2024).gamma(W_DEFAULT * 40 + 1e-9, 1.0)
        w = w / w.sum()
        score = np.where(mask, np.array([row['expression'], row['specificity'], row['reproducibility'], row['rnai'], row['correlation'], row['neural_enriched'], row['neural_specificity']]), 0) @ w
        all_scores_dir[d] = score
    
    dirichlet_median = np.median(all_scores_dir)
    
    # Uniform Dirichlet
    all_scores_unif = np.zeros(1000)
    for d in range(1000):
        w = np.random.default_rng(2024).dirichlet(np.ones(7))
        score = np.where(mask, np.array([row['expression'], row['specificity'], row['reproducibility'], row['rnai'], row['correlation'], row['neural_enriched'], row['neural_specificity']]), 0) @ w
        all_scores_unif[d] = score
    uniform_median = np.median(all_scores_unif)
    
    # Corrected Dirichlet (renormalized per draw)
    corrected_scores = np.zeros(1000)
    for d in range(1000):
        w = np.random.default_rng(2024).gamma(W_DEFAULT * 40 + 1e-9, 1.0)
        w = w / w.sum()
        # Renormalize over available streams
        mask = ~np.isnan([row['expression'], row['specificity'], row['reproducibility'], row['rnai'], row['correlation'], row['neural_enriched'], row['neural_specificity']])
        w_avail = w[np.where(mask)[0]]
        w_renorm = w[mask] / w[mask].sum() if w[mask].sum() > 0 else w[mask]
        vals = np.array([row['expression'], row['specificity'], row['reproducibility'], row['rnai'], row['correlation'], row['neural_enriched'], row['neural_specificity']])
        vals = np.where(np.isnan(vals), 0, vals)
        score = vals[mask] @ w[mask] / w[mask].sum()
        corrected_scores[d] = score
    
    corrected_median = np.median(corrected_scores)
    
    n_avail = sum(mask)
    
    results.append({
        'gene_id': row['gene_id'],
        'gene_name': row['gene_name'],
        'n_avail': int(n_avail),
        'fixed_score': fixed_score,
        'dirichlet_median': dirichlet_median,
        'uniform_median': uniform_median,
        'corrected_median': corrected_median,
        'diff_current_corrected': dirichlet_median - corrected_median
    })

results_df = pd.DataFrame(results)
results_df.to_csv('D:/Bioinformatics/analysis_results.csv', index=False)
print("Done!")
print(results_df.groupby('n_avail').agg({'diff_current_corrected': ['mean', 'median', 'std'], 'fixed_score': 'mean', 'dirichlet_median': 'mean', 'corrected_median': 'mean'}))