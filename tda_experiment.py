import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
import kmapper as km
from sklearn.cluster import DBSCAN, KMeans, AgglomerativeClustering, OPTICS, SpectralClustering
from ripser import ripser
from persim import plot_diagrams

# -----------------------------
# loading and preprocessing our csv file
# -----------------------------
data_path = "drugs.csv"
df = pd.read_csv(data_path)

df.replace(r'^\s*$', np.nan, regex=True, inplace=True)
df.fillna(0, inplace=True)
#drop columns 
columns_to_drop = [
    'solubility predicted (DrugBank in mg/ml)',
    'Toxicity Rat DrugBank (mol/kg)',
    'Hepatoxicity Likelihood Score (PubChem) A,B,C,D,E = 1,2,3,4,5',
    '(HEPATO/Drug Induced Liver Injury Toxicity Grade) toxicity (Pubchem)',
    'pka (Strongest acidic) (Drugbank Predicted)',
    'pKa (Strongest basic) (DrugBank Predicted)',
    'aromatic rings (chembl)',
    'atom count (RDkit)',
    'Caco-2 factor (DrugBank)',
    'protein binding factor Min(not all will have; DrugBank) (DECIMAL)',
    'protein binding factor Max (not all will have; DrugBank) (Decimal)',
    'blood-brain barrier (DrugBank)',
    'signs for blood-brain barrier(+ = 1, - = 0) (DrugBank)']

df.drop(columns=columns_to_drop, inplace=True)

X = df.select_dtypes(include=[np.number])

X.columns = X.columns.astype(str)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

# -----------------------------
# trying out dimensionality reduction with tsnse
# -----------------------------
tsne = TSNE(n_components=2, perplexity=40, learning_rate=200, n_iter=1000, random_state=42)
tsne_results = tsne.fit_transform(X_scaled)

df['TSNE1'] = tsne_results[:, 0]
df['TSNE2'] = tsne_results[:, 1]

plt.figure(figsize=(10, 7))
if 'Disease' in df.columns:
    unique_categories = df['Disease'].unique()
    color_map = plt.cm.get_cmap('tab20', len(unique_categories))  
    category_map = {category: color_map(i) for i, category in enumerate(unique_categories)}
    for category, color in category_map.items():
        subset = df[df['Disease'] == category]
        plt.scatter(subset['TSNE1'], subset['TSNE2'], label=category, alpha=0.7, c=[color])
    plt.legend(title="Disease Categories")
else:
    plt.scatter(df['TSNE1'], df['TSNE2'], alpha=0.7)
plt.title("t-SNE Visualization")
plt.xlabel("Component 1")
plt.ylabel("Component 2")
plt.show()


# -----------------------------
# mapper algorithm testing and visualization
# -----------------------------
mapper = km.KeplerMapper(verbose=1)
projected_data = mapper.fit_transform(X_scaled, projection=TSNE(n_components=2, perplexity=40, learning_rate=200, n_iter=1000, random_state=42))
graph = mapper.map(projected_data, clusterer=DBSCAN(eps=100, min_samples=3))
mapper.visualize(graph, path_html="mapper.html", title="Mapper Visualization with DBSCAN")
print("Mapper visualization saved to 'mapper.html'")

# KMeans
graph = mapper.map(projected_data, clusterer=KMeans(n_clusters=5, random_state=42))
mapper.visualize(graph, path_html="mapper_kmeans.html", title="Mapper Visualization with KMeans")
print("Mapper visualization saved to 'mapper_kmeans.html'")

# Agglomerative Clustering
graph = mapper.map(projected_data, clusterer=AgglomerativeClustering(n_clusters=5))
mapper.visualize(graph, path_html="mapper_agglomerative.html", title="Mapper Visualization with Agglomerative Clustering")
print("Mapper visualization saved to 'mapper_agglomerative.html'")

# OPTICS
graph = mapper.map(projected_data, clusterer=OPTICS(min_samples=5, max_eps=100))
mapper.visualize(graph, path_html="mapper_optics.html", title="Mapper Visualization with OPTICS")
print("Mapper visualization saved to 'mapper_optics.html'")

# Spectral Clustering
graph = mapper.map(projected_data, clusterer=SpectralClustering(n_clusters=5, random_state=42))
mapper.visualize(graph, path_html="mapper_spectral.html", title="Mapper Visualization with Spectral Clustering")
print("Mapper visualization saved to 'mapper_spectral.html'")



# -----------------------------
# trying out different properties with the mapper and saving it in folder cluster_analysis
# -----------------------------

properties = [
    'MolecularWeight', 'XLogP', 'ExactMass', 'MonoisotopicMass', 'TPSA', 'HBondDonorCount',
    'HBondAcceptorCount', 'RotatableBondCount', 'HeavyAtomCount', 'AtomStereoCount', 'DefinedAtomStereoCount',
    'UndefinedAtomStereoCount', 'BondStereoCount', 'DefinedBondStereoCount', 'UndefinedBondStereoCount',
    'CovalentUnitCount', 'Volume3D', 'XStericQuadrupole3D', 'YStericQuadrupole3D', 'ZStericQuadrupole3D',
    'FeatureAcceptorCount3D', 'FeatureDonorCount3D', 'FeatureAnionCount3D', 'FeatureCationCount3D', 
    'FeatureRingCount3D', 'FeatureHydrophobeCount3D', 'ConformerModelRMSD3D', 'EffectiveRotorCount3D'
]

for property in properties:
    color_values = df[property].values.flatten()
    path_html = f"cluster_analysis/mapper_{property.lower()}.html"  
    title = f"TDA: {property} Analysis"

    mapper.visualize(graph, 
                     path_html=path_html,
                     title=title, 
                     custom_tooltips=df['Drug Name'].values, 
                     color_values=color_values,
                     color_function_name=[property])
    
    print(f"Mapper visualization saved to '{path_html}'")

# -----------------------------
# persistent homology testing and visualization
# -----------------------------
ph_results = ripser(X_scaled, maxdim=2)
diagrams = ph_results['dgms']
plot_diagrams(diagrams, show=True)

