import importlib.resources

def get_package_location():
    return importlib.resources.files("molecule_info")