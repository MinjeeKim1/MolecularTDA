# Install UMAP if needed
!pip install umap-learn

# Imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import umap.umap_ as umap
from sklearn.cluster import KMeans
from sklearn.preprocessing import RobustScaler
from matplotlib.patches import Ellipse

# Load data
df = pd.read_csv('/content/finaldata62425.csv')

# --- Update this column name if needed ---
label_col = "Disease"  # or "DrugClass" or whatever label column you have

# Select numeric data
data_numeric = df.select_dtypes(include=[np.number]).copy()
data_filled = data_numeric.fillna(data_numeric.mean())

# RobustScaler to reduce outlier influence
scaler = RobustScaler()
data_scaled = scaler.fit_transform(data_filled)

# UMAP
reducer = umap.UMAP(n_neighbors=15, min_dist=0.1, random_state=42)
embedding = reducer.fit_transform(data_scaled)

# KMeans clustering (k=3)
kmeans = KMeans(n_clusters=3, random_state=42)
clusters = kmeans.fit_predict(embedding)

# Prepare plot
plt.figure(figsize=(10, 8))

# Disease label colors
unique_labels = df[label_col].dropna().unique()
colors = plt.cm.get_cmap('tab10', len(unique_labels))

# Plot points by disease label
for i, label in enumerate(unique_labels):
    mask = df[label_col] == label
    plt.scatter(
        embedding[mask, 0], embedding[mask, 1],
        s=20, alpha=0.7, label=label,
        color=colors(i)
    )

# Add dotted cluster outlines (KMeans)
for i in range(3):
    points = embedding[clusters == i]
    center = points.mean(axis=0)
    stds = np.std(points, axis=0)
    ellipse = Ellipse(
        xy=center,
        width=stds[0]*4, height=stds[1]*4,
        edgecolor='black', facecolor='none',
        linestyle='dotted', linewidth=2
    )
    plt.gca().add_patch(ellipse)

# Final plot
plt.title("UMAP of Drugs Colored by Disease, with Dotted Cluster Outlines (k=3)")
plt.xlabel("UMAP 1")
plt.ylabel("UMAP 2")
plt.legend(title="Disease", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)
plt.tight_layout()
plt.show()
