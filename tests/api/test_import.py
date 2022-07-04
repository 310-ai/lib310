def test_lib310_import():
    try:
        from lib310 import data
        from lib310 import db
        from lib310 import ml
        from lib310 import plot
        from lib310 import tools
    except ImportError as e:
        raise e
