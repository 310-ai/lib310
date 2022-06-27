def test_b10_import():
    try:
        from b10 import db
        from b10 import MLAPI
        from b10 import pl
        from b10 import tl
    except ImportError as e:
        raise e
