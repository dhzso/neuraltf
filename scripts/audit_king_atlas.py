import pandas as pd

df = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\data\king_atlas.tsv', sep='\t')

print('=== 1. TOTAL ROW COUNT ===')
print(len(df))
print()

print('=== 2. ALL UNIQUE compartment VALUES ===')
for c in sorted(df['compartment'].unique()):
    print('  -', c)
print()

print('=== 3. ALL UNIQUE subcluster VALUES ===')
for s in sorted(df['subcluster'].unique()):
    print('  -', s)
print()

print('=== 4. COUNT ROWS PER compartment ===')
for k, v in sorted(df['compartment'].value_counts().items()):
    print(f'  {k}: {v}')
print()

print('=== 5. COUNT ROWS PER neural subcluster ===')
neural_subs = [s for s in df['subcluster'].unique() if str(s).startswith('neural')]
print('Total neural subclusters:', len(neural_subs))
for s in sorted(neural_subs):
    count = (df['subcluster'] == s).sum()
    print(f'  {s}: {count} rows')
print()

print('=== 6. TOTAL UNIQUE v6_ids ===')
print(df['v6_id'].nunique())
print()

print('=== 7. NEURAL SUBCLUSTER FULL LIST ===')
print(neural_subs)
print()

print('=== 8. CELL_TYPE FOR NEURAL-RELATED ROWS ===')
ndf = df[df['subcluster'].astype(str).str.startswith('neural')]
seen = set()
for _, row in ndf[['subcluster', 'cell_type']].drop_duplicates().iterrows():
    ct = row['cell_type'] if str(row['cell_type']) != 'nan' else '(empty)'
    if (row['subcluster'], ct) not in seen:
        seen.add((row['subcluster'], ct))
        print(f'  {row["subcluster"]} -> {ct}')
print()

print('=== 9. SUBCLUSTER -> CELL_TYPE MAPPING ===')
groups = ndf.groupby('subcluster')['cell_type'].apply(
    lambda x: sorted(set([v for v in x if pd.notna(v) and str(v).strip() != '']))
)
for sub in sorted(groups.index):
    vals = ', '.join(groups[sub]) if groups[sub] else '(empty)'
    print(f'  {sub} => {vals}')
print()

print('=== 10. log2fc AND pval RANGES ===')
all_log = pd.to_numeric(df['log2fc'], errors='coerce').dropna()
all_pv = pd.to_numeric(df['pval'], errors='coerce').dropna()
print(f'OVERALL:    log2fc [{all_log.min():.4f}, {all_log.max():.4f}], pval [{all_pv.min():.3e}, {all_pv.max():.3e}]')

n_log = pd.to_numeric(ndf['log2fc'], errors='coerce').dropna()
n_pv = pd.to_numeric(ndf['pval'], errors='coerce').dropna()
nn = df[~df['subcluster'].astype(str).str.startswith('neural')]
nn_log = pd.to_numeric(nn['log2fc'], errors='coerce').dropna()
nn_pv = pd.to_numeric(nn['pval'], errors='coerce').dropna()

print(f'NEURAL:     log2fc [{n_log.min():.4f}, {n_log.max():.4f}] mean={n_log.mean():.4f}, pval [{n_pv.min():.3e}, {n_pv.max():.3e}]')
print(f'NON-NEURAL: log2fc [{nn_log.min():.4f}, {nn_log.max():.4f}] mean={nn_log.mean():.4f}, pval [{nn_pv.min():.3e}, {nn_pv.max():.3e}]')
print()

print('=== 11. GENES WITH log2fc > 2 IN AT LEAST ONE NEURAL SUBCLUSTER ===')
high = ndf[pd.to_numeric(ndf['log2fc'], errors='coerce') > 2]
print(f'Rows: {len(high)}')
print(f'Unique v6_ids: {high["v6_id"].nunique()}')
print(f'Unique gene_names: {high["gene_name"].nunique()}')
print('Gene names:')
for g in sorted(high['gene_name'].unique()):
    print(f'  - {g}')

print()
print('=== 12. X1 compartment - Neural-related subclusters ===')
x1_neural = df[(df['compartment'] == 'X1') & df['subcluster'].apply(lambda x: str(x).startswith('neural'))]
print(f'Rows in X1 with neural subclusters: {len(x1_neural)}')
print('Subclusters:')
for s in sorted(x1_neural['subcluster'].unique()):
    rows = x1_neural[x1_neural['subcluster'] == s]
    cts = sorted(set([v for v in rows['cell_type'] if pd.notna(v)]))
    print(f'  {s} -> cell_type: {cts if cts else "(empty)"}')

print()
print('=== 13. X1 Major Tissue Neural rows ===')
x1mt = df[(df['compartment'] == 'X1 Major Tissue') & (df['subcluster'] == 'Neural')]
print(f'Rows: {len(x1mt)}, Unique genes: {x1mt["v6_id"].nunique()}')
print('Genes:')
for _, r in x1mt.iterrows():
    print(f'  {r["gene_name"]} log2fc={r["log2fc"]}')

print()
print('=== 14. Full Fincher cluster mapping (fincher_cluster for neural rows) ===')
print('subcluster -> fincher_cluster (unique values):')
ficher_map = ndf.groupby('subcluster')['fincher_cluster'].apply(
    lambda x: sorted(set([v for v in x if pd.notna(v) and str(v).strip() != '']))
)
for sub in sorted(ficher_map.index):
    vals = ', '.join(ficher_map[sub]) if ficher_map[sub] else '(empty)'
    print(f'  {sub} -> {vals}')