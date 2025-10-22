import pandas as pd
import numpy as np
import base64
import kmapper as km
from sklearn.cluster import DBSCAN
from sklearn.manifold import TSNE
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import itertools
import warnings
warnings.filterwarnings('ignore')

N_CUBES_RANGE = [10, 13, 15, 20]
PERC_OVERLAP_RANGE = [0.2, 0.3, 0.4, 0.5]
DBSCAN_EPS_RANGE = [10, 15, 20, 25]
DBSCAN_MIN_SAMPLES_RANGE = [5]

TSNE_PERPLEXITY = 40
TSNE_LEARNING_RATE = 200
TSNE_N_ITER = 1000
RANDOM_STATE = 42
TSNE_METRIC = 'cosine'

df = pd.read_csv("fingering.csv")
fingerprints_base64 = df["Fingerprint2D"].dropna().tolist()

def decode_base64_fingerprint(fp_str):
    fp_bytes = base64.b64decode(fp_str)           
    bitstring = ''.join(f"{byte:08b}" for byte in fp_bytes)  
    return np.array([int(bit) for bit in bitstring], dtype=np.uint8)

X = np.array([decode_base64_fingerprint(fp) for fp in fingerprints_base64])
print("Decoded fingerprint matrix shape:", X.shape)

scaler = StandardScaler()
X_scaled = scaler.fit_transform(X)

print("Applying t-SNE dimensionality reduction...")
tsne = TSNE(
    n_components=2,
    perplexity=TSNE_PERPLEXITY,
    learning_rate=TSNE_LEARNING_RATE,
    max_iter=TSNE_N_ITER,
    random_state=RANDOM_STATE,
    metric=TSNE_METRIC
)
lens = tsne.fit_transform(X_scaled)
print("t-SNE projection shape:", lens.shape)

def optimize_parameters(X_scaled, lens, disease_labels):
    best_score = -1
    best_params = None
    results = []
    
    print("\nOptimizing parameters using silhouette score...")
    print("Testing parameter combinations...")
    
    total_combinations = len(N_CUBES_RANGE) * len(PERC_OVERLAP_RANGE) * len(DBSCAN_EPS_RANGE) * len(DBSCAN_MIN_SAMPLES_RANGE)
    current_combination = 0
    
    for n_cubes, perc_overlap, eps, min_samples in itertools.product(
        N_CUBES_RANGE, PERC_OVERLAP_RANGE, DBSCAN_EPS_RANGE, DBSCAN_MIN_SAMPLES_RANGE
    ):
        current_combination += 1
        print(f"Testing combination {current_combination}/{total_combinations}: "
              f"n_cubes={n_cubes}, overlap={perc_overlap}, eps={eps}, min_samples={min_samples}")
        
        try:
            temp_mapper = km.KeplerMapper(verbose=0)
            
            temp_graph = temp_mapper.map(
                lens,
                X_scaled,
                cover=km.Cover(n_cubes=n_cubes, perc_overlap=perc_overlap),
                clusterer=DBSCAN(eps=eps, min_samples=min_samples)
            )
            
            cluster_labels = []
            for i in range(len(X_scaled)):
                found_cluster = -1
                for node_id, members in temp_graph['nodes'].items():
                    if i in members:
                        found_cluster = int(node_id.split('_')[0].replace('cube', ''))
                        break
                cluster_labels.append(found_cluster)
            
            cluster_labels = np.array(cluster_labels)
            
            unique_clusters = np.unique(cluster_labels)
            n_clusters = len(unique_clusters[unique_clusters != -1])
            
            if n_clusters > 1 and n_clusters < len(X_scaled):
                sil_score = silhouette_score(X_scaled, cluster_labels)
                
                results.append({
                    'n_cubes': n_cubes,
                    'perc_overlap': perc_overlap,
                    'eps': eps,
                    'min_samples': min_samples,
                    'silhouette_score': sil_score,
                    'n_clusters': n_clusters,
                    'n_nodes': len(temp_graph['nodes']),
                    'n_edges': len(temp_graph['links'])
                })
                
                print(f"  Silhouette score: {sil_score:.4f}, Clusters: {n_clusters}, Nodes: {len(temp_graph['nodes'])}")
                
                if sil_score > best_score:
                    best_score = sil_score
                    best_params = (n_cubes, perc_overlap, eps, min_samples)
                    print(f"  *** New best score: {best_score:.4f} ***")
            else:
                print(f"  Skipped (invalid clustering: {n_clusters} clusters)")
                
        except Exception as e:
            print(f"  Error: {e}")
            continue
    
    print(f"\nOptimization complete!")
    print(f"Best parameters: n_cubes={best_params[0]}, overlap={best_params[1]}, eps={best_params[2]}, min_samples={best_params[3]}")
    print(f"Best silhouette score: {best_score:.4f}")
    
    return best_params, best_score, results

disease_labels = df.loc[df["Fingerprint2D"].notna(), "Disease"]
drug_names = df.loc[df["Fingerprint2D"].notna(), "Drug Name"]

best_params, best_score, optimization_results = optimize_parameters(X_scaled, lens, disease_labels)
N_CUBES, PERC_OVERLAP, DBSCAN_EPS, DBSCAN_MIN_SAMPLES = best_params

results_df = pd.DataFrame(optimization_results)
results_df = results_df.sort_values('silhouette_score', ascending=False)
results_df.to_csv('parameter_optimization_results.csv', index=False)
print(f"Optimization results saved to parameter_optimization_results.csv")

print(f"\nUsing optimized parameters:")
print(f"N_CUBES = {N_CUBES}")
print(f"PERC_OVERLAP = {PERC_OVERLAP}")
print(f"DBSCAN_EPS = {DBSCAN_EPS}")
print(f"DBSCAN_MIN_SAMPLES = {DBSCAN_MIN_SAMPLES}")

mapper = km.KeplerMapper(verbose=1) 

disease_to_int = {d: i for i, d in enumerate(disease_labels.unique())}
color_values = disease_labels.map(disease_to_int).values

graph = mapper.map(
    lens,
    X_scaled,
    cover=km.Cover(n_cubes=N_CUBES, perc_overlap=PERC_OVERLAP),  
    clusterer=DBSCAN(eps=DBSCAN_EPS, min_samples=DBSCAN_MIN_SAMPLES)
)

mapper.visualize(
    graph,
    path_html="mapper_fingerprint_tsne.html",
    title="Molecular Fingerprint Mapper Graph (t-SNE, Colored by Disease)",
    custom_tooltips=drug_names.values,
    color_values=color_values,
    color_function_name="Disease"
)

print("Mapper graph saved as mapper_fingerprint_tsne.html")
