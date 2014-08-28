print ("hello from init.py")

def mydir(x=None):
    if x==None:
        return []
    else:
        return [i for i in dir(x) if not '__' in i]

def startConsole():
    import rlcompleter, readline # optional, will allow Up/Down/History in the console
    import code
    readline.parse_and_bind("tab:complete") # autocomplete
    vars = globals()
    vars.update(locals())
    shell = code.InteractiveConsole(vars)
    shell.interact()


startConsole()


