# lib310 Python package [![CircleCI](https://dl.circleci.com/status-badge/img/gh/310-ai/lib310/tree/master.svg?style=shield)](https://dl.circleci.com/status-badge/redirect/gh/310-ai/lib310/tree/master)

![Alt Text](media/logo.png)

## Sample Usages

### Protein Functional Annotation

```python
# 1. import lib310 
import lib310

# 2. Get Spike SARS2 related proteins from database
seqs = lib310.db.fetch(
    name="SPIKE_SARS2",
    feature='sequence',
    limit=500
)

# 3. Instantiate a pre-trained GO Annotation machine learning model (e.g. TALE)
goa = lib310.ml.goa(model="TALE", v="512_756_4layers")

# 4. Predict!
annotations = goa.predict(seqs)

# 5. (Optional) Get sequences' embeddings for visualization purposes or downstream analysis 
embeddings = goa.get_embeddings(seqs, layer=-1)
```

### Protein Generation

```python
# 1. import lib310 
import lib310

# 2. Instantiate a pre-trained Generative Machine Learning model (e.g. GPT3)
lm = lib310.ml.AutoRegressiveLM(model="GPT3", v="512_756_4layers")

# 3. Predict!
generated_sequences = lm.run(num_samples=1024)

# 4. Downstream Analysis...
clusters = lib310.tools.cluster(generated_sequences, method='kcluster')

# 5. Visualization
lib310.plot.umap(generated_sequences, clusters)
```



