print ("Hello from init.py")
import ngfem_cf

def mydir(x=None):
    if x==None:
        return []
    else:
        return [i for i in dir(x) if not '__' in i]

def execfile(fname):
    exec(open(fname).read())

def startConsole():
    import code
    try:
        import readline
        import rlcompleter
        readline.parse_and_bind("tab:complete") # autocomplete
    except:
        try:
            import pyreadline as readline
            import rlcompleter
            readline.parse_and_bind("tab:complete") # autocomplete
        except:
            print('readline not found')
    vars = globals()
    vars.update(locals())
    shell = code.InteractiveConsole(vars)
    shell.interact()


startConsole()


