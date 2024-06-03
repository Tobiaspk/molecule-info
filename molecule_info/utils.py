import importlib.resources

def get_package_location():
    try:
        # support python 3.9
        return importlib.resources.files("molecule_info")
    except:
        # earlier python versions
        import os
        return os.path.dirname(os.path.dirname(__file__))