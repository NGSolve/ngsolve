def dynamic_metadata(field, settings):

    if field != "dependencies":
        msg = "Only the 'dependencies' field is supported"
        raise ValueError(msg)

    if settings:
        msg = "No inline configuration is supported"
        raise ValueError(msg)

    import tomllib
    import importlib.metadata

    config = tomllib.load(open("pyproject.toml", "rb"))
    dependencies = config["project"].get("dependencies", [])

    ngs_version = importlib.metadata.version("ngsolve")
    dependencies = dependencies + [f"ngsolve=={ngs_version}"]
    return dependencies


def get_requires_for_dynamic_metadata(_settings):
    return []
