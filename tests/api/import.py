def test_b10_import():
    try:
        from b10 import dataAPI as d, MLAPI as m
    except ImportError as e:
        raise e
