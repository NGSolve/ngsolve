def mydir(x=None):
    if x==None:
        return []
    else:
        return [i for i in dir(x) if not '__' in i]

def startConsole():
    import readline # optional, will allow Up/Down/History in the console
    import code
    vars = globals()
    vars.update(locals())
    shell = code.InteractiveConsole(vars)
    shell.interact()


startConsole()


