import sys
try:
    del sys.modules['lib310']
    del sys.modules['google']
    del sys.modules['google.cloud']
    del sys.modules['google.logging']
    del sys.modules['numpy']
    del sys.modules['scanpy']
    del sys.modules['transformers']
    del sys.modules['torch']
    del sys.modules['pandas']
    del sys.modules['seaborn']
    del sys.modules['matplotlib']
    del sys.modules['plotly']
    del sys.modules['graphviz']

except Exception as e:
    print(e)
    print("Module was not cached")
import lib310
