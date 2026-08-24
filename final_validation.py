import pandas as pd
import numpy as np

df_centered = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\dirichlet_overall_top10_byscore.csv')
df_uniform = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\dirichlet_uniform_overall_top10.csv')
df_fixed = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\top10_neural_tfs_prioritized.csv')
df_old = pd.read_csv(r'D:\Bioinformatics\analysis_results.csv')

df_rank = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\runs\pipeline_run\rank_neural.csv')

results = []
for _, row in pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\runs\pipeline_run\rank_neural.csv').iterrows():
    gid = row['gene_id']
    gene_name = row['gene_name']
    
    n_avail = 7 - sum(pd.isna(row[s]) for s in ['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity'])
    
    fixed = df_fixed[df_fixed['gene_id_v6']==gid]['composite_score'].values
    fixed = fixed[0] if len(fixed) > 0 else np.nan
    
    centered_new = df_centered[df_centered['gene_id_v6']==gid]['dirichlet_median_score'].values
    centered_new = centered_new[0] if len(centered_new) > 0 else np.nan
    
    uniform_new = df_uniform[df_uniform['gene_id_v6']==gid]['uniform_median_score'].values
    uniform_new = uniform_new[0] if len(uniform_new) > 0 else np.nan
    
    old_centered = df_old[df_old['gene_id']==gid]['dirichlet_median'].values
    old_centered = old_centered[0] if len(old_centered) > 0 else np.nan
    
    old_uniform = df_old[df_old['gene_id']==gid]['uniform_median'].values
    old_uniform = old_uniform[0] if len(old_uniform) > 0 else np.nan
    
    old_corrected = df_old[df_old['gene_id']==gid]['corrected_median'].values
    old_corrected = old_corrected[0] if len(old_corrected) > 0 else np.nan
    
    centered_new = df_centered[df_centered['gene_id_v6']==gid]['dirichlet_median_score'].values
    centered_new = centered_new[0] if len(df_centered[df_centered['gene_id_v6']==gid]) > 0 else np.nan
    
    uniform_new = df_uniform[df_uniform['gene_id_v6']==gid]['uniform_median_score'].values
    uniform_new = uniform_new[0] if len(uniform_new) > 0 else np.nan
    
    old_centered = df_old[df_old['gene_id']==gid]['dirichlet_median'].values
    old_centered = old_centered[0] if len(old_centered) > 0 else np.nan
    
    old_uniform = df_old[df_old['gene_id']==gid]['uniform_median'].values
    old_uniform = old_uniform[0] if len(old_uniform) > 0 else np.nan
    
    old_corrected = df_old[df_old['gene_id']==gid]['corrected_median'].values
    old_corrected = old_corrected[0] if len(old_corrected) > 0 else np.nan
    
    n_avail = 7 - sum(pd.isna(row[s]) for s in ['expression', 'specificity', 'reproducibility', 'rnai', 'correlation', 'neural_enriched', 'neural_specificity'])
    
    diff_centered_new_vs_old = centered_new - old_centered if not np.isnan(centered_new) and not np.isnan(old_centered) else np.nan
    diff_new_vs_corrected = centered_new - old_corrected if not np.isnan(centered_new) and not np.isnan(old_corrected) else np.nan
    
    match_centered = np.isclose(centered_new, old_corrected, rtol=1e-2) if not np.isnan(centered_new) and not np.isnan(old_corrected) else False
    
    results.append({
        'gene_id': gid,
        'gene_name': gene_name,
        'n_avail': n_avail,
        'centered_new': centered_new,
        'centered_old': old_centered,
        'centered_old_corrected': old_corrected,
        'match_centered': match_centered,
        'diff_centered_new_vs_old': centered_new - old_centered if not np.isnan(centered_new) and not np.isnan(old_centered) else np.nan,
        'diff_new_vs_corrected': centered_new - old_corrected if not np.isnan(centered_new) and not np.isnan(old_corrected) else np.nan,
    })

df_results = pd.DataFrame(results)
df_results.to_csv('D:/Bioinformatics/corrected_comparison.csv', index=False)

print('=== CORRECTED DIRICHLET VALIDATION ===')
print('Total candidates:', len(df_results))
print('Candidates with 7/7 streams:', sum(df_results['n_avail']==7))
print('Candidates with 6/7 streams:', sum(df_results['n_avail']==6))
print('Candidates with 5/7 streams:', sum(df_results['n_avail']==5))
print(f'Candidates with 4/7 streams: {sum(df_results["n_avail"]==4)}')
print()

matches = sum(1 for r in results if r['match_centered'])
print(f'Candidates where new centered matches old corrected: {sum(1 for r in results if r["match_centered"])}/{len(results)}')
print()

diffs = [r['centered_new'] - r['centered_old_corrected'] for r in results if not np.isnan(r['centered_new']) and not np.isnan(r['centered_old_corrected'])]
print(f'Diff (new vs corrected): mean={np.mean(diffs):.6f}, max={np.max(diffs):.6f}')

for n in [7, 6, 5, 4]:
    subset = [r for r in results if r['n_avail'] == n]
    if subset:
        diffs = [r['centered_new'] - r['centered_old_corrected'] for r in subset if not np.isnan(r['centered_new']) and not np.isnan(r['centered_old_corrected'])]
        if diffs:
            print(f'{n}/7 avail: n={len(subset)}, mean_diff={np.mean(diffs):.6f}, max={max(diffs):.6f}')

fixed_top10 = set(pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\top10_neural_tfs_prioritized.csv')['gene_id_v6'])
centered_top10 = set(pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\dirichlet_overall_top10_byscore.csv')['gene_id_v6'])
uniform_top10 = set(df_uniform[df_uniform['gene_id_v6']==gid]['uniform_median_score'].values)

fixed_top10 = set(pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\top10_neural_tfs_prioritized.csv')['gene_id_v6'])
centered_top10 = set(df_centered['gene_id_v6'])
uniform_top10 = set(df_uniform['gene_id_v6'])

print('\nFixed top-10:', len(fixed_top10))
print('Centered top-10:', len(centered_top10))
print('Uniform top-10:', len(uniform_top10))
print(f'Fixed intersect Centered: {len(fixed_top10 & centered_top10)}/10')
print(f'Fixed intersect Uniform: {len(fixed_top10 & uniform_top10)}/10')
print(f'Centered intersect Uniform: {len(centered_top10 & uniform_top10)}/10')

print('\nMissingness-stratified difference (new vs corrected):')
for n in [7, 6, 5, 4]:
    subset = [r for r in results if r['n_avail'] == n]
    if subset:
        diffs = [r['centered_new'] - r['centered_old_corrected'] for r in subset if not np.isnan(r['centered_new']) and not np.isnan(r['centered_old_corrected'])]
        if diffs:
            print(f'  {n}/7 avail: n={len(subset)}, mean_diff={np.mean(diffs):.6f}, max={max(diffs):.6f}')