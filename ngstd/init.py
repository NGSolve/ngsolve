print ("Hello from init.py")
import code
import ngfem_cf
import ngmpi
import sys

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

class InteractiveMPIConsole(code.InteractiveConsole):
    def runsource(self, source, filename="<input>", symbol="single"):
        ngmpi.SendCommand('ngs_py '+source)
        code.InteractiveConsole.runsource(self, source, filename, symbol)
        self.Barrier()
    def interact(self):
        self.write("MPI Shell\n")
        self.write("================\n")
        self.write("Use pprint(str) to print with MPI ranks\n\n")
        self.runsource("from ngmpi import *")
        self.runsource("pprint=lambda p='':print('Process ' + str(ngmpi.Rank()) + '\\n'+str(p)+'\\n')")
        self.locals.update(locals())
        self.locals.update(globals())
        sys.ps1 = "MPI >>> "
        sys.ps2 = "MPI ... "
        code.InteractiveConsole.interact(self,'')
    def Barrier(self):
        source = 'ngmpi.Barrier()'
        ngmpi.SendCommand('ngs_py '+source)
        code.InteractiveConsole.runsource(self, source)


def MpiShell():
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
    mshell = InteractiveMPIConsole(vars)
    mshell.interact()

startConsole()


