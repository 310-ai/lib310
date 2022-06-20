# B10 Python package

![Alt Text](media/logo.)

## Sample Usages

### Protein Functional Annotation

```python
# 1. import b10 
from b10 import dataAPI as d, MLAPI as m

# 2. Get Spike SARS2 related proteins from database
seqs = d.fetch(
    name="SPIKE_SARS2",
    feature='sequence',
    limit=500
)

# 3. Instantiate a pre-trained GO Annotation machine learning model (e.g. TALE)
goa = m.goa(model="TALE", v="512_756_4layers")

# 4. Predict!
annotations = goa.predict(seqs)

# 5. (Optional) Get sequences' embeddings for visualization purposes or downstream analysis 
embeddings = goa.get_embeddings(seqs, layer=-1)
```

### Protein Generation

```python
# 1. import b10 
import b10

# 2. Instantiate a pre-trained Generative Machine Learning model (e.g. GPT3)
goa = b10.MLAPI.plm(model="GPT3", v="512_756_4layers")

# 4. Predict!
generated_sequences = goa.generate(num_samples=1024)

# 5. Downstream Analysis...
clusters = b10.tl.cluster(generated_sequences, method='kcluster')
```



