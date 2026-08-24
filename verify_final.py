import pandas as pd

df_centered = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\dirichlet_overall_top10_byscore.csv')
df_uniform = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\dirichlet_uniform_overall_top10.csv')
df_fixed = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\top10_neural_tfs_prioritized.csv')

print('=== CORRECTED DIRICHLET VALIDATION SUMMARY ===')
print(f'Centered Dirichlet top-10: {len(pd.read_csv(r"D:\\Bioinformatics\\projects\\NeuralTF\\results\\dirichlet_overall_top10_byscore.csv"))}')
print(f'Uniform Dirichlet top-10: {len(pd.read_csv(r"D:\\Bioinformatics\\projects\\NeuralTF\\results\\dirichlet_uniform_overall_top10.csv"))}')
print(f'Fixed-weight top-10: {len(pd.read_csv(r"D:\\Bioinformatics\\projects\\NeuralTF\\results\\top10_neural_tfs_prioritized.csv"))}')

print()
print('Top candidate dd14115:')
df_centered = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\dirichlet_overall_top10_byscore.csv')
df_uniform = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\dirichlet_uniform_overall_top10.csv')
df_fixed = pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\top10_neural_tfs_prioritized.csv')

print(f'  Fixed composite: {df_fixed[df_fixed["gene_id_v6"]=="dd_Smed_v6_14115_0_1"]["composite_score"].values[0]:.4f}')
print(f'  Centered Dirichlet: {df_centered[df_centered["gene_id_v6"]=="dd_Smed_v6_14115_0_1"]["dirichlet_median_score"].values[0]:.4f}')
print(f'  Uniform Dirichlet: {df_uniform[df_uniform["gene_id_v6"]=="dd_Smed_v6_14115_0_1"]["uniform_median_score"].values[0]:.4f}')

# Check overlap
fixed_top10 = set(pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\top10_neural_tfs_prioritized.csv')['gene_id_v6'])
centered_top10 = set(pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\dirichlet_overall_top10_byscore.csv')['gene_id_v6'])
uniform_top10 = set(pd.read_csv(r'D:\Bioinformatics\projects\NeuralTF\results\dirichlet_uniform_overall_top10.csv')['gene_id_v6'])

print()
print(f'Fixed top-10: {len(fixed_top10)}')
print(f'Centered top-10: {len(centered_top10)}')
print(f'Uniform top-10: {len(uniform_top10)}')
print(f'Fixed ∩ Centered: {len(fixed_top10 & centered_top10)}/10')
print(f'Fixed ∩ Uniform: {len(fixed_top10 & uniform_top10)}/10')
print(f'Centered ∩ Uniform: {len(centered_top10 & uniform_top10)}/10')