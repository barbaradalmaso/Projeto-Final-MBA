############# TCC DE MBA EM DATA SCIENCE & ANALYTICS #############
###################### UNIVERSIDADE DE SÃO PAULO #################
###################### ESTUDANTE BARBARA DALMASO #################

# In[0.1]: Importar pacotes
# pip install numpy pandas matplotlib scanpy seaborn scipy anndata scikit-learn
# pip install git+https://github.com/theislab/scanpy

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import os
from scipy.sparse import csr_matrix
import anndata as ad
import scanpy.external as sce


# In[0.2]: Ajustar work-directory e fazer download dos dados
os.getcwd()
os.chdir("/Users/barbaradalmaso/Desktop/AMD-RPE-spatial/Dados/")

# Pré-processamento dos dados
# Download dos datasets
wet = sc.read_h5ad("/Users/barbaradalmaso/Desktop/AMD-RPE-spatial/Dados/GSE221042_RAW/wet_amd.h5ad") 
dry = sc.read_h5ad("/Users/barbaradalmaso/Desktop/AMD-RPE-spatial/Dados/GSE221042_RAW/dry_amd.h5ad") 
control = sc.read_h5ad("/Users/barbaradalmaso/Desktop/AMD-RPE-spatial/Dados/GSE221042_RAW/control.h5ad") 

# Transformar matrix em fator
wet.obs['matrix'] = wet.obs['matrix'].astype('category')
dry.obs['matrix'] = dry.obs['matrix'].astype('category')
control.obs['matrix'] = control.obs['matrix'].astype('category')

wet.obs['patient_code'] = wet.obs['matrix']
dry.obs['patient_code'] = dry.obs['matrix']
control.obs['patient_code'] = control.obs['matrix']

# Primeiro, como observei que a maior parte das celulas possuem baixo count de genes, decidi filtrar 
# as celulas com total count abaixo de 500 (em cada amostra)
control.obs['total_counts'] = control.X.sum(axis=1).A1 if isinstance(control.X, csr_matrix) else control.X.sum(axis=1)
control_filtered = control[(control.obs['total_counts'] > 1400),:]
print(f"Número de células antes da filtragem: {control.n_obs}")
print(f"Número de células após a filtragem: {control_filtered.n_obs}")

dry.obs['total_counts'] = dry.X.sum(axis=1).A1 if isinstance(dry.X, csr_matrix) else dry.X.sum(axis=1)
dry_filtered = dry[(dry.obs['total_counts'] > 1400), :]
print(f"Número de células antes da filtragem: {dry.n_obs}")
print(f"Número de células após a filtragem: {dry_filtered.n_obs}")

wet.obs['total_counts'] = wet.X.sum(axis=1).A1 if isinstance(wet.X, csr_matrix) else wet.X.sum(axis=1)
wet_filtered = wet[(wet.obs['total_counts'] > 1400), :]
print(f"Número de células antes da filtragem: {wet.n_obs}")
print(f"Número de células após a filtragem: {wet_filtered.n_obs}")

# Como tirei as celulas que nao possuem nenhuma contagem de expressao genica, o numero de celulas por amostra reduziu drasticamente.
# Agora eh possivel que eu junte todas as amostras de pacientes em um unico dataset, pois nao ficara mais um arquivo tao pesado.
# Mas para isso, necessito fazer metadados de qualidade para nao perder o agrupamento de pacientes.

# Adicionar coluna 'group' aos metadados de cada objeto
control_filtered.obs['group'] = 'control'
dry_filtered.obs['group'] = 'dry_amd'
wet_filtered.obs['group'] = 'wet_amd'

# Combinar os objetos AnnData
combined_data = ad.concat([control_filtered, dry_filtered, wet_filtered], label="group", keys=["control", "dry_amd", "wet_amd"])

# Normalizar
sc.pp.normalize_total(combined_data, target_sum=1000)
sc.pp.log1p(combined_data)

# Salvar o objeto combinado em um arquivo .h5ad
output_path = "/Users/barbaradalmaso/Desktop/AMD-RPE-spatial/Dados/combined_patient_sc.h5ad"
combined_data.write(output_path)

# In[0.3]: Download do dataset combinado e já normalizado
# Download dos dados
bulkdata = sc.read_h5ad("/Users/barbaradalmaso/Desktop/AMD-RPE-spatial/Dados/combined_patient_sc.h5ad") # Arquivo de dados scRNA-seq

# Filtrando genes com baixa expressão
sc.pp.highly_variable_genes(bulkdata, min_mean=0.5, max_mean=5, min_disp=0.5)

# In[1.1] Análise por PCA
# Padronizar os dados (z-score)
sc.pp.scale(bulkdata)

# Executar PCA
sc.tl.pca(bulkdata, n_comps=20)  # Ajuste o número de componentes conforme necessário

# Extrair as coordenadas de PC1 e PC2
pc1 = bulkdata.obsm['X_pca'][:, 0]
pc2 = bulkdata.obsm['X_pca'][:, 1]

group_colors = {
    'control': '#8ECEDA',  # Azul-esverdeado
    'dry_amd': '#F6D6C2',  # Bege claro
    'wet_amd': '#D47264'   # Vermelho terroso
}

# Plotar PCA com coloração personalizada
fig, ax = plt.subplots(figsize=(6, 4))
sc.pl.pca(bulkdata, color="group", show=False, ax=ax, palette=group_colors, s = 10)

# Adicionar título com valores de PC1 e PC2
pc1_var = bulkdata.uns['pca']['variance_ratio'][0] * 100  # Variância explicada por PC1
pc2_var = bulkdata.uns['pca']['variance_ratio'][1] * 100  # Variância explicada por PC2
ax.set_title(f'PCA - PC1: {pc1_var:.2f}% | PC2: {pc2_var:.2f}%', fontsize=16)

# Ajustar a opacidade (alpha=0.5)
for line in ax.lines:
    line.set_alpha(0.1)
    
ax.set_xlabel('PC1', fontsize=16)
ax.set_ylabel('PC2', fontsize=16)

# Mostrar o gráfico
plt.tight_layout()
plt.show()

# Visualizar a porcentagem de variância explicada por cada componente principal
sc.pl.pca_variance_ratio(bulkdata, log=True, show=True)

# Aplicar Harmony para corrigir batch effect
sce.pp.harmony_integrate(bulkdata, key='patient_code')  # Ajuste a chave para a variável correta de batch

# Substituir os PCs originais pelos PCs corrigidos
bulkdata.obsm['X_pca'] = bulkdata.obsm['X_pca_harmony']

# Rodar novamente o código de PCA para ver o PCA ajustado após a correção de batch

# In[1.2] Análise por UMAP
# Calcular a matriz de vizinhança (necessário para UMAP)
sc.pp.neighbors(bulkdata, n_neighbors=20, n_pcs=10)

# Executar UMAP
sc.tl.umap(bulkdata)  # Realiza a redução de dimensionalidade usando UMAP

# Plotar UMAP com coloração personalizada
fig, ax = plt.subplots(figsize=(6, 4))
sc.pl.umap(bulkdata, color="group", show=False, ax=ax, palette=group_colors, s=10)

# Adicionar título com valores de UMAP1 e UMAP2
ax.set_title('UMAP', fontsize=16)

# Ajustar a opacidade (alpha=0.5)
for line in ax.lines:
    line.set_alpha(0.9)

ax.set_xlabel('UMAP1', fontsize=16)
ax.set_ylabel('UMAP2', fontsize=16)

# Mostrar o gráfico
plt.tight_layout()
plt.show()

# In[1.3] Análise por t-SNE
# Calcular a matriz de vizinhança (necessário para t-SNE)
sc.pp.neighbors(bulkdata, n_neighbors=100, n_pcs=10)

# Executar t-SNE
sc.tl.tsne(bulkdata)  # Realiza a redução de dimensionalidade usando t-SNE

# Plotar t-SNE com coloração personalizada
fig, ax = plt.subplots(figsize=(6, 4))
sc.pl.tsne(bulkdata, color="group", show=False, ax=ax, palette=group_colors, s=7)

# Adicionar título com valores de t-SNE1 e t-SNE2
ax.set_title('t-SNE', fontsize=16)

# Ajustar a opacidade (alpha=0.5)
for line in ax.lines:
    line.set_alpha(0.5)

ax.set_xlabel('t-SNE1', fontsize=16)
ax.set_ylabel('t-SNE2', fontsize=16)

# Mostrar o gráfico
plt.tight_layout()
plt.show()

# In[1.4] Métricas para comparar os três métodos
# Métricas para comparar os três métodos
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from sklearn.preprocessing import StandardScaler
import numpy as np

scaler = StandardScaler()
bulkdata.obsm["X_pca"] = scaler.fit_transform(bulkdata.obsm["X_pca"])
bulkdata.obsm["X_umap"] = scaler.fit_transform(bulkdata.obsm["X_umap"])
bulkdata.obsm["X_tsne"] = scaler.fit_transform(bulkdata.obsm["X_tsne"])

# Converter grupos em valores numéricos diretamente
labels = bulkdata.obs["group"].astype("category").cat.codes.values

# Amostragem aleatória para reduzir consumo de memória (opcional)
sample_size = min(10000, bulkdata.shape[0])  # Limita a 10.000 pontos para otimizar

def compute_metrics(embedding, labels, sample_size=10000):
    idx = np.random.choice(len(labels), size=sample_size, replace=False)  # Seleção aleatória
    sample_embedding = embedding[idx]
    sample_labels = labels[idx]
    
    silhouette = silhouette_score(sample_embedding, sample_labels)
    davies_bouldin = davies_bouldin_score(sample_embedding, sample_labels)
    calinski_harabasz = calinski_harabasz_score(sample_embedding, sample_labels)
    
    return silhouette, davies_bouldin, calinski_harabasz

# PCA - 2 primeiras componentes principais
pca_silhouette, pca_davies, pca_calinski = compute_metrics(bulkdata.obsm["X_pca"], labels)
print(f"PCA - Silhouette Score: {pca_silhouette:.4f}, Davies-Bouldin: {pca_davies:.4f}, Calinski-Harabasz: {pca_calinski:.4f}")

# UMAP - Coordenadas da projeção UMAP
umap_silhouette, umap_davies, umap_calinski = compute_metrics(bulkdata.obsm["X_umap"], labels)
print(f"UMAP - Silhouette Score: {umap_silhouette:.4f}, Davies-Bouldin: {umap_davies:.4f}, Calinski-Harabasz: {umap_calinski:.4f}")

# t-SNE - Coordenadas da projeção t-SNE
tsne_silhouette, tsne_davies, tsne_calinski = compute_metrics(bulkdata.obsm["X_tsne"], labels)
print(f"t-SNE - Silhouette Score: {tsne_silhouette:.4f}, Davies-Bouldin: {tsne_davies:.4f}, Calinski-Harabasz: {tsne_calinski:.4f}")

# In[1.5] O silhouete score disse que MENTIexiste uma grande sobreposição dos dados,
# Silhouette Score - PCA: -0.0131
# Silhouette Score - UMAP: -0.0247
# Silhouette Score - t-SNE: -0.0210
# e isso está sendo observado qualitativamente. Vou colorir os gráficos pelos 
# diferentes clusters presentes no dataset, pra ver se essa dispersão é biológica
# Download dos dados
bulkdata = sc.read_h5ad("/Users/barbaradalmaso/Desktop/AMD-RPE-spatial/Dados/combined_patient_sc.h5ad") # Arquivo de dados scRNA-seq

# Filtrando genes com baixa expressão
sc.pp.highly_variable_genes(bulkdata, min_mean=0.5, max_mean=5, min_disp=0.5)

# 1. Calcular a matriz de vizinhança (necessário para análise de agrupamento)
sc.pp.neighbors(bulkdata, n_neighbors=10, n_pcs=10)

# 2. Rodar o método de agrupamento (louvain ou leiden)
# Usando o método 'leiden' (ou 'louvain', caso preferir)
sc.tl.leiden(bulkdata, resolution=1.0)  # Você pode ajustar a resolução conforme necessário

# 3. Agora, pode fazer a análise de genes diferencialmente expressos por grupo
sc.tl.rank_genes_groups(bulkdata, groupby="leiden", method="wilcoxon")  # 'leiden' é a coluna gerada com o agrupamento

# 4. Visualizar os genes mais diferencialmente expressos
sc.pl.rank_genes_groups(bulkdata, n_genes=20, sharey=False)

# 5. Obter os top 5 genes de cada cluster
markers = pd.DataFrame(bulkdata.uns["rank_genes_groups"]["names"]).head(5)
markers

# Refazer os gráficos de dimensionalidade, colorindo pelos clustersimport
# UMAP
sc.tl.umap(bulkdata)  # Realiza a redução de dimensionalidade usando UMAP
fig, ax = plt.subplots(figsize=(6, 4))
sc.pl.umap(bulkdata, color='leiden', show=False, ax=ax, palette='tab20', s=10)  # Colorir por clusters de Leiden
ax.set_title('UMAP - Clusters', fontsize=16)
ax.set_xlabel('UMAP1', fontsize=16)
ax.set_ylabel('UMAP2', fontsize=16)
plt.tight_layout()
plt.show()

# t-SNE
sc.tl.tsne(bulkdata)  # Realiza a redução de dimensionalidade usando t-SNE
fig, ax = plt.subplots(figsize=(6, 4))
sc.pl.tsne(bulkdata, color='leiden', show=False, ax=ax, palette='tab20', s=10)  # Colorir por clusters de Leiden
ax.set_title('t-SNE - Clusters', fontsize=16)
ax.set_xlabel('t-SNE1', fontsize=16)
ax.set_ylabel('t-SNE2', fontsize=16)
plt.tight_layout()
plt.show()

# PCA
sc.tl.pca(bulkdata)  # Realiza a redução de dimensionalidade usando PCA
fig, ax = plt.subplots(figsize=(6, 4))
sc.pl.pca(bulkdata, color='group', show=False, ax=ax, palette='tab20', s=10)  # Colorir por clusters de Leiden
ax.set_title('PCA - Clusters', fontsize=16)
ax.set_xlabel('PCA1', fontsize=16)
ax.set_ylabel('PCA2', fontsize=16)
plt.tight_layout()
plt.show()

# In[1.5] Fazer análise de métricas novamente, separando por célula
from sklearn.metrics import silhouette_score, davies_bouldin_score, calinski_harabasz_score
from sklearn.preprocessing import StandardScaler
import numpy as np

scaler = StandardScaler()
bulkdata.obsm["X_pca"] = scaler.fit_transform(bulkdata.obsm["X_pca"])
bulkdata.obsm["X_umap"] = scaler.fit_transform(bulkdata.obsm["X_umap"])
bulkdata.obsm["X_tsne"] = scaler.fit_transform(bulkdata.obsm["X_tsne"])

# Converter grupos em valores numéricos diretamente
labels = bulkdata.obs["leiden"].astype("category").cat.codes.values

# Amostragem aleatória para reduzir consumo de memória (opcional)
sample_size = min(10000, bulkdata.shape[0])  # Limita a 10.000 pontos para otimizar

def compute_metrics(embedding, labels, sample_size=10000):
    idx = np.random.choice(len(labels), size=sample_size, replace=False)  # Seleção aleatória
    sample_embedding = embedding[idx]
    sample_labels = labels[idx]
    
    silhouette = silhouette_score(sample_embedding, sample_labels)
    davies_bouldin = davies_bouldin_score(sample_embedding, sample_labels)
    calinski_harabasz = calinski_harabasz_score(sample_embedding, sample_labels)
    
    return silhouette, davies_bouldin, calinski_harabasz

# PCA - 2 primeiras componentes principais
pca_silhouette, pca_davies, pca_calinski = compute_metrics(bulkdata.obsm["X_pca"], labels)
print(f"PCA - Silhouette Score: {pca_silhouette:.4f}, Davies-Bouldin: {pca_davies:.4f}, Calinski-Harabasz: {pca_calinski:.4f}")

# UMAP - Coordenadas da projeção UMAP
umap_silhouette, umap_davies, umap_calinski = compute_metrics(bulkdata.obsm["X_umap"], labels)
print(f"UMAP - Silhouette Score: {umap_silhouette:.4f}, Davies-Bouldin: {umap_davies:.4f}, Calinski-Harabasz: {umap_calinski:.4f}")

# t-SNE - Coordenadas da projeção t-SNE
tsne_silhouette, tsne_davies, tsne_calinski = compute_metrics(bulkdata.obsm["X_tsne"], labels)
print(f"t-SNE - Silhouette Score: {tsne_silhouette:.4f}, Davies-Bouldin: {tsne_davies:.4f}, Calinski-Harabasz: {tsne_calinski:.4f}")

# In[2.0] Clusterização de fotorreceptores em diferentes condições

# Olhando manualmente os marcadores de cada grupo do Leiden, vou selecionar especificamente
# os que acredito serem fotorreceptores. Inicialmente, irei selecioná-los
selected_clusters = ['3', '4']  # Garantindo que os clusters são strings

# Filtrando células pertencentes aos clusters escolhidos
subset_adata = bulkdata[bulkdata.obs['leiden'].isin(selected_clusters)].copy()

# Substituindo os valores de 'leiden' pelos nomes dos fotorreceptores
subset_adata.obs['leiden'] = subset_adata.obs['leiden'].replace({'3': 'Bastonetes', '4': 'Bastonetes'})

# Conferindo distribuição dos clusters com as substituições
print(subset_adata.obs['leiden'].value_counts())

# In[2.05] Avaliação da clusterização de novo
scaler = StandardScaler()
subset_adata.obsm["X_umap"] = scaler.fit_transform(subset_adata.obsm["X_umap"])

# Converter grupos em valores numéricos diretamente
labels = subset_adata.obs["group"].astype("category").cat.codes.values

# Definir tamanho da amostra sem ultrapassar o total de células disponíveis
sample_size = min(10000, subset_adata.shape[0])

def compute_metrics(embedding, labels, sample_size):
    if sample_size < len(labels):  # Evita erro caso o subconjunto seja menor que a amostra definida
        idx = np.random.choice(len(labels), size=sample_size, replace=False)  # Seleção aleatória
        sample_embedding = embedding[idx]
        sample_labels = labels[idx]
    else:
        sample_embedding = embedding
        sample_labels = labels
    
    silhouette = silhouette_score(sample_embedding, sample_labels)
    davies_bouldin = davies_bouldin_score(sample_embedding, sample_labels)
    calinski_harabasz = calinski_harabasz_score(sample_embedding, sample_labels)
    
    return silhouette, davies_bouldin, calinski_harabasz

# UMAP - Coordenadas da projeção UMAP
umap_silhouette, umap_davies, umap_calinski = compute_metrics(subset_adata.obsm["X_umap"], labels, sample_size)
print(f"UMAP - Silhouette Score: {umap_silhouette:.4f}, Davies-Bouldin: {umap_davies:.4f}, Calinski-Harabasz: {umap_calinski:.4f}")


# In[2.1] Aplicar o método k-means de clusterização
# Selecionar o valor de k pelo Elbow
from sklearn.cluster import KMeans

inertia = []
K_range = range(2, 10)

for k in K_range:
    kmeans = KMeans(n_clusters=k, random_state=42)
    kmeans.fit(subset_adata.obsm['X_pca'])
    inertia.append(kmeans.inertia_)

plt.plot(K_range, inertia, marker='o')
plt.xlabel('Número de clusters')
plt.ylabel('Inércia (soma dos quadrados das distâncias)')
plt.title('Método do Cotovelo para Escolher k')
plt.show()

# Pelo Elbow, parece que o melhor é 4 para n de clusters
# Usando PCA para reduzir a dimensionalidade antes do KMeans
sc.tl.pca(subset_adata, svd_solver='arpack')

# Aplicando K-Means nos primeiros componentes principais
kmeans = KMeans(n_clusters=3, random_state=100)
subset_adata.obs['kmeans'] = kmeans.fit_predict(subset_adata.obsm['X_pca'])

# Definir os vizinhos
sc.pp.neighbors(subset_adata, n_neighbors=50, use_rep='X_pca')  # Ajuste n_neighbors conforme necessário
sc.tl.umap(subset_adata, min_dist=0.1) # Numeros menores deixam os pontos mais próximos
sc.tl.umap(subset_adata)

# Criar o gráfico UMAP
subset_adata.obs['kmeans'] = subset_adata.obs['kmeans'].astype(str)
fig, ax = plt.subplots(figsize=(5, 4))
# Plotar o gráfico UMAP colorindo pelos clusters K-means
sc.pl.umap(subset_adata, color='kmeans', ax=ax, show=False, s=25, palette='tab20')
# Personalizar o gráfico
ax.set_title("K-means Clusters")
ax.set_xlabel('UMAP 1')
ax.set_ylabel('UMAP 2')

# Exibir o gráfico
plt.show()

# t-sne
fig, ax = plt.subplots(figsize=(5, 4))
sc.pl.tsne(subset_adata, color='kmeans', show=False, ax=ax, palette='tab20', s=25)
ax.set_title('t-SNE - Clusters', fontsize=16)
ax.set_xlabel('t-SNE1', fontsize=16)
ax.set_ylabel('t-SNE2', fontsize=16)
plt.tight_layout()
plt.show()


# In[2.2] Fazer o gráfico com colorindo pelas células de cada grupo
# Criar o gráfico UMAP
fig, ax = plt.subplots(figsize=(5, 4))

# Plotar o gráfico UMAP colorindo pelos clusters K-means
sc.pl.umap(subset_adata, color='group', palette='tab20', ax=ax, show=False, s=25)

# Personalizar o gráfico
ax.set_title("Condição de saúde")
ax.set_xlabel('UMAP 1')
ax.set_ylabel('UMAP 2')

# Exibir o gráfico
plt.show()

fig, ax = plt.subplots(figsize=(5, 4))
sc.pl.tsne(subset_adata, color='group', show=False, ax=ax, palette='tab20', s=25)
ax.set_title('t-SNE - Clusters', fontsize=16)
ax.set_xlabel('t-SNE1', fontsize=16)
ax.set_ylabel('t-SNE2', fontsize=16)
plt.tight_layout()
plt.show()

# In[2.3] Fazer o gráfico com colorindo pelas células obtidas em Leiden
sc.tl.leiden(subset_adata, resolution=0.25)  # Você pode ajustar a resolução conforme necessário
sc.tl.rank_genes_groups(subset_adata, groupby="leiden", method="wilcoxon") 
# 4. Visualizar os genes mais diferencialmente expressos
sc.pl.rank_genes_groups(subset_adata, n_genes=20, sharey=False)

# 5. Obter os top 5 genes de cada cluster
markers = pd.DataFrame(subset_adata.uns["rank_genes_groups"]["names"]).head(5)
markers

# Refazer os gráficos de dimensionalidade, colorindo pelos clustersimport
# UMAP
sc.tl.umap(subset_adata)  # Realiza a redução de dimensionalidade usando UMAP
fig, ax = plt.subplots(figsize=(5, 4))
sc.pl.umap(subset_adata, color='leiden', show=False, ax=ax, palette='tab20', s=25)  # Colorir por clusters de Leiden
ax.set_title('Leiden - Clusters', fontsize=16)
ax.set_xlabel('UMAP1')
ax.set_ylabel('UMAP2')
plt.tight_layout()
plt.show()

fig, ax = plt.subplots(figsize=(5, 4))
sc.pl.tsne(subset_adata, color='leiden', show=False, ax=ax, palette='tab20', s=25)
ax.set_title('t-SNE - Clusters', fontsize=16)
ax.set_xlabel('t-SNE1', fontsize=16)
ax.set_ylabel('t-SNE2', fontsize=16)
plt.tight_layout()
plt.show()

import scanpy as sc
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.mixture import GaussianMixture

# 1. Reduzindo a dimensionalidade com PCA
sc.tl.pca(subset_adata, svd_solver='arpack')

# 2. Aplicando GMM
n_clusters = 3  # Defina o número de clusters desejado
gmm = GaussianMixture(n_components=n_clusters, random_state=42)
subset_adata.obs['gmm'] = gmm.fit_predict(subset_adata.obsm['X_pca']).astype(str)  # Convertendo para string para ser categórico

# 3. Definir vizinhos e calcular UMAP
sc.pp.neighbors(subset_adata, n_neighbors=50, use_rep='X_pca')  
sc.tl.umap(subset_adata, min_dist=0.3, random_state=42)  # Ajuste os parâmetros conforme necessário

# 4. Criar gráfico UMAP colorindo pelos clusters GMM
fig, ax = plt.subplots(figsize=(5,4))
sc.pl.umap(subset_adata, color='gmm', show=False, ax=ax, palette='tab20', s=25)  
ax.set_title('GMM - Clusters', fontsize=16)
ax.set_xlabel('UMAP1')
ax.set_ylabel('UMAP2')
plt.tight_layout()
plt.show()

# Agora que fiz as clusterizações, vou coletar os genes marcadores de cada grupo, para ver subgrupos
# Identificar genes diferenciais para os clusters K-means
sc.tl.rank_genes_groups(subset_adata, groupby="kmeans", method="wilcoxon")
markers_kmeans = pd.DataFrame(subset_adata.uns["rank_genes_groups"]["names"]).head(5)

# Identificar genes diferenciais para os clusters Leiden
sc.tl.rank_genes_groups(subset_adata, groupby="leiden", method="wilcoxon")
markers_leiden = pd.DataFrame(subset_adata.uns["rank_genes_groups"]["names"]).head(5)

# Identificar genes diferenciais para os clusters GMM
sc.tl.rank_genes_groups(subset_adata, groupby="gmm", method="wilcoxon")
markers_gmm = pd.DataFrame(subset_adata.uns["rank_genes_groups"]["names"]).head(5)

# Identificar genes diferenciais para os grupos de condição de saúde (Group)
sc.tl.rank_genes_groups(subset_adata, groupby="group", method="wilcoxon")
markers_group = pd.DataFrame(subset_adata.uns["rank_genes_groups"]["names"]).head(5)

print("🔹 Marcadores K-Means")
print(markers_kmeans)

print("\n🔹 Marcadores Leiden")
print(markers_leiden)

print("\n🔹 Marcadores GMM")
print(markers_gmm)

print("\n🔹 Marcadores Group")
print(markers_group)

# In[2.4] Agora,vamos fazer matrizes de confusão pra tentar capturar a qualidade dos agrupamentos
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import confusion_matrix

# Criar um dicionário de mapeamento
cluster_mapping = {0: 'dry_amd', 2: 'control', 1: 'wet_amd'}

# Aplicar a conversão em Leiden e K-means
subset_adata.obs['leiden_named'] = subset_adata.obs['leiden'].astype(int).map(cluster_mapping)
subset_adata.obs['kmeans_named'] = subset_adata.obs['kmeans'].astype(int).map(cluster_mapping)

# Criar DataFrame para análise
df = subset_adata.obs[['group', 'leiden_named', 'kmeans_named']].dropna()

# Definir classes
classes = ['dry_amd', 'control', 'wet_amd']

# Criar matrizes de confusão
conf_matrix_leiden = confusion_matrix(df['group'], df['leiden_named'], labels=classes)
conf_matrix_kmeans = confusion_matrix(df['group'], df['kmeans_named'], labels=classes)

# Função para plotar matriz de confusão
def plot_confusion_matrix(conf_matrix, title, xlabel, ylabel):
    plt.figure(figsize=(5, 4))
    sns.heatmap(conf_matrix, annot=True, fmt='d', cmap="Blues",
                xticklabels=classes, yticklabels=classes, cbar=True, linewidths=0)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.show()

# Plotar matrizes de confusão sem grades
plot_confusion_matrix(conf_matrix_leiden, "Matriz de confusão - Leiden", "Leiden Cluster", "Real")
plot_confusion_matrix(conf_matrix_kmeans, "Matriz de confusão - K-means", "K-means Cluster", "Real")

# Função para calcular métricas
def compute_metrics(conf_matrix, method_name):
    accuracy = np.trace(conf_matrix) / np.sum(conf_matrix)
    precision_per_class = np.diag(conf_matrix) / np.sum(conf_matrix, axis=0)
    recall_per_class = np.diag(conf_matrix) / np.sum(conf_matrix, axis=1)
    f1_per_class = 2 * (precision_per_class * recall_per_class) / (precision_per_class + recall_per_class)

    # Criar DataFrame para melhor visualização
    metrics_df = pd.DataFrame({
        "Classe": classes,
        "Precisão": precision_per_class,
        "Recall": recall_per_class,
        "F1-score": f1_per_class
    })

    print(f"\nMétricas para {method_name}:")
    print(f"Acurácia: {accuracy:.4f}\n")
    print(metrics_df.to_string(index=False))

# Calcular métricas para Leiden e K-means
compute_metrics(conf_matrix_leiden, "Leiden")
compute_metrics(conf_matrix_kmeans, "K-means")



# In[3.0] Agora vamos iniciar o treinamento do modelo pra identificar os grupos celulares
# Salvar o objeto com dados de fotorreceptores em um arquivo .h5ad
output_path = "/Users/barbaradalmaso/Desktop/AMD-RPE-spatial/Dados/fotorreceptores_sc.h5ad"
subset_adata.write(output_path)

# In[3.1] Coletar os dados que serão usados especificamente para o treinamento
dados = sc.read_h5ad("/Users/barbaradalmaso/Desktop/AMD-RPE-spatial/Dados/fotorreceptores_sc.h5ad")

# In[3.2] # Importar os pacotes necessários para análise de random forrests
import scanpy as sc
import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, accuracy_score
from sklearn.model_selection import GridSearchCV
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import confusion_matrix

# In[3.3] # Preciso organizar meus dados para ter uma variável X (com as features - no caso os genes)
# e uma variável Y (com a variável algo - que no caso serão as condições a serem previstas)
# 1. Separar os dados:
    # Matriz de expressão gênica
X = dados.X.toarray()  # Converter para numpy array se estiver como matriz esparsa

    # Variável alvo (grupo de cada célula)
y = dados.obs['group'].astype(str)  # Certificar que é string para classificação

    # Verificar a distribuição das classes
print(y.value_counts())  

# 2. Separar os dados em treino e teste - 80% deles serão treino, e 20% serão teste
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

print(f"Tamanho do treino: {X_train.shape}, Tamanho do teste: {X_test.shape}")

# 3. Buscar quais os melhores hiperparametros que encontram o melhor modelo
# Definir o modelo base
rf = RandomForestClassifier(random_state=42, class_weight='balanced')

# Definir a grade de parâmetros para teste
param_grid = {
    'n_estimators': [50, 100, 200],  # Testar diferentes quantidades de árvores
    'max_depth': [None, 10, 20],  # Permitir árvores mais profundas ou não
    'min_samples_split': [2, 5, 10],  # Evita overfitting
    'min_samples_leaf': [1, 2, 5],  # Número mínimo de amostras por folha
    'max_features': ['sqrt', 'log2']  # Define quantas features cada nó pode usar
}

# Criar o GridSearchCV para encontrar os melhores parâmetros
grid_search = GridSearchCV(rf, param_grid, cv=5, scoring='accuracy', n_jobs=-1, verbose=2)

# Rodar o grid search
grid_search.fit(X_train, y_train)

# Exibir os melhores hiperparâmetros encontrados
print("Melhores hiperparâmetros:", grid_search.best_params_)

# Treinar o modelo final com os melhores hiperparâmetros
best_rf = grid_search.best_estimator_

# Fazer previsões
y_pred = best_rf.predict(X_test)

# Avaliar o modelo ajustado
accuracy = accuracy_score(y_test, y_pred)
print(f"Acurácia do modelo otimizado: {accuracy:.2f}")

print(classification_report(y_test, y_pred))

# 4 Treinar e avaliar o modelo final Random Forest
# Fitting 5 folds for each of 162 candidates, totalling 810 fits
# Melhores hiperparâmetros: {'max_depth': None, 'max_features': 'sqrt', 'min_samples_leaf': 5, 'min_samples_split': 2, 'n_estimators': 200}
# Acurácia do modelo otimizado: 0.88
    # Criar o modelo Random Forest
rf_model_otimizado = RandomForestClassifier(
    max_depth=None,
    max_features='sqrt',
    min_samples_leaf=5,
    min_samples_split=2,
    n_estimators=200,
    random_state=42,
    n_jobs=-1,
    class_weight='balanced'
)
    # Treinar o modelo
rf_model_otimizado.fit(X_train, y_train)

    # Fazer previsões no conjunto de teste
y_pred_otimizado = rf_model_otimizado.predict(X_test)

    # Avaliar a acurácia
accuracy = accuracy_score(y_test, y_pred_otimizado)
print(f"Acurácia do modelo otimizado: {accuracy:.2f}")

# Gerar relatório de classificação
print(classification_report(y_test, y_pred_otimizado))


# 6. Analisar a Importância dos Genes
# Obter a importância das features
importances = rf_model_otimizado.feature_importances_

# Criar um dataframe com os nomes dos genes e sua importância
genes = dados.var_names  # Nome dos genes
feature_importance = pd.DataFrame({"Gene": genes, "Importância": importances})
feature_importance = feature_importance.sort_values(by="Importância", ascending=False)

# Mostrar os top 10 genes mais importantes
print(feature_importance.head(10))

# Plotar os 10 genes mais importantes
plt.figure(figsize=(6, 5))
plt.barh(feature_importance['Gene'][:10], feature_importance['Importância'][:10], color='royalblue')
sns.despine()
plt.xlabel("Importância")
plt.ylabel("Gene")
plt.title("Top 10 Genes mais importantes (Random Forest)")
plt.gca().invert_yaxis()  # Inverter para exibir o mais importante no topo
plt.show()

# Criar e plotar matriz de confusão
cm = confusion_matrix(y_test, y_pred_otimizado)
plt.figure(figsize=(6,5))
sns.heatmap(cm, annot=True, fmt="d", cmap="Blues", xticklabels=best_rf.classes_, yticklabels=best_rf.classes_)
plt.xlabel("Previsto")
plt.ylabel("Real")
plt.title("Matriz de Confusão - Random Forest")
plt.show()

# In[3.3] # Fazer o mesmo utilizando XGBoost
# Importar pacotes
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.metrics import classification_report, accuracy_score, confusion_matrix
from xgboost import XGBClassifier
import matplotlib.pyplot as plt
import seaborn as sns


# In[3.4] # Rodar o modelo
# 1. Separar os dados
X = dados.X.toarray()  # Converter para array numpy
y = dados.obs['group'].astype("category").cat.codes  # Converter para números

# Verificar a distribuição das classes
print(pd.Series(y).value_counts())

# 2. Separar treino e teste (80% treino, 20% teste)
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)
print(f"Tamanho do treino: {X_train.shape}, Tamanho do teste: {X_test.shape}")

# 3. Definir modelo base do XGBoost
xgb_model = XGBClassifier(
    objective="multi:softmax",  # Classificação multiclasse
    num_class=len(set(y)),  # Número de classes
    eval_metric="mlogloss",
    use_label_encoder=False,
    random_state=42 # Usa GPU se disponível (se não, remove essa linha)
)

# 4. Definir espaço de busca para otimização
param_grid = {
    'n_estimators': [50, 100],  
    'max_depth': [3, 6],  
    'learning_rate': [0.01, 0.1],  
    'subsample': [0.8],  
    'colsample_bytree': [0.8]
}

# 5. Criar e rodar o RandomizedSearchCV (mais rápido que GridSearchCV)
random_search = RandomizedSearchCV(
    xgb_model, param_distributions=param_grid, 
    n_iter=10, scoring="accuracy", n_jobs=-1, cv=3, verbose=2, random_state=42
)

random_search.fit(X_train, y_train)

# 6. Exibir os melhores hiperparâmetros
print("Melhores hiperparâmetros:", random_search.best_params_)

# 7. Treinar o modelo final com os melhores hiperparâmetros
best_xgb = random_search.best_estimator_
y_pred = best_xgb.predict(X_test)

# 8. Avaliação do modelo
accuracy = accuracy_score(y_test, y_pred)
print(f"Acurácia do modelo otimizado: {accuracy:.2f}")
print(classification_report(y_test, y_pred))

# 9. Criar matriz de confusão
plt.figure(figsize=(6, 5))
sns.heatmap(confusion_matrix(y_test, y_pred), annot=True, fmt='d', cmap="Blues", xticklabels=dados.obs['group'].cat.categories, yticklabels=dados.obs['group'].cat.categories)
plt.xlabel("Previsto")
plt.ylabel("Real")
plt.title("Matriz de Confusão - XGBoost")
plt.show()

# Fazer o gráfico com os genes
# Obter a importância dos genes no modelo treinado
importances = best_xgb.feature_importances_

# Criar um DataFrame com os nomes dos genes e a importância
genes = dados.var_names  # Nomes dos genes
feature_importance = pd.DataFrame({"Gene": genes, "Importância": importances})

# Ordenar pela importância (maior primeiro)
feature_importance = feature_importance.sort_values(by="Importância", ascending=False)

# Plotar os 10 genes mais importantes
plt.figure(figsize=(6, 5))
sns.barplot(x="Importância", y="Gene", data=feature_importance.head(10), color='royalblue')
sns.despine()
plt.xlabel("Importância")
plt.ylabel("Gene")
plt.title("Top 10 Genes Mais Importantes (XGBoost)")
plt.show()


# In[3.3] # Análise LightGBM
# Importar pacotes
import lightgbm as lgb
from sklearn.model_selection import train_test_split, RandomizedSearchCV
from sklearn.metrics import accuracy_score, classification_report, confusion_matrix
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# 1. Separar os dados
X = dados.X.toarray()  # Converter para array se for esparso
y = dados.obs['group'].astype('category').cat.codes  # Converter para categorias numéricas

# 2. Dividir em treino e teste
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42, stratify=y)

# 3. Definir o modelo base
lgb_model = lgb.LGBMClassifier(boosting_type='gbdt', objective='multiclass', random_state=42, n_jobs=-1)

# 4. Definir hiperparâmetros para busca randômica
param_dist = {
    'num_leaves': [20, 31, 40],
    'learning_rate': [0.01, 0.05, 0.1],
    'n_estimators': [50, 100, 200],
    'max_depth': [-1, 10, 20],
    'subsample': [0.7, 0.8, 0.9],
    'colsample_bytree': [0.7, 0.8, 1.0]
}

# 5. Rodar RandomizedSearchCV para otimizar
random_search = RandomizedSearchCV(lgb_model, param_distributions=param_dist, n_iter=20, cv=3, scoring='accuracy', verbose=1, n_jobs=-1, random_state=42)
random_search.fit(X_train, y_train)

# 6. Melhor modelo
tuned_lgb = random_search.best_estimator_
print("Melhores hiperparâmetros:", random_search.best_params_)

# 7. Avaliação no conjunto de teste
y_pred = tuned_lgb.predict(X_test)
accuracy = accuracy_score(y_test, y_pred)
print(f"Acurácia do modelo otimizado: {accuracy:.2f}")
print(classification_report(y_test, y_pred))

# 8. Matriz de confusão
plt.figure(figsize=(6,5))
sns.heatmap(confusion_matrix(y_test, y_pred), annot=True, fmt='d', cmap='Blues')
plt.xlabel("Previsto")
plt.ylabel("Real")
plt.title("Matriz de Confusão - LightGBM")
plt.show()

# 9. Importância dos genes
importances = tuned_lgb.feature_importances_
genes = dados.var_names
feature_importance = pd.DataFrame({'Gene': genes, 'Importância': importances}).sort_values(by='Importância', ascending=False)

# 10. Gráfico de importância
plt.figure(figsize=(6,5))
sns.barplot(x=feature_importance.Importância[:10], y=feature_importance.Gene[:10], color='royalblue')
plt.xlabel("Importância")
plt.ylabel("Gene")
plt.title("Top 10 Genes Mais Importantes - LightGBM")
sns.despine()
plt.show()
