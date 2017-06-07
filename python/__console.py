import code
# import ngfemcf
# import ngmpi

def mydir(x=None):
    if x==None:
        return []
    else:
        return [i for i in dir(x) if not '__' in i]

def execfile(fname):
    exec(open(fname).read())

class InteractiveMPIConsole(code.InteractiveConsole):
# Function copied from /usr/lib/python3.4/code.py line 38 
    def runsource(self, source, filename="<input>", symbol="single"):
        try:
            compiled_code = self.compile(source, filename, symbol)
        except (OverflowError, SyntaxError, ValueError):
            # Case 1
            self.showsyntaxerror(filename)
            return False

        if compiled_code is None:
            # Case 2
            return True

        # Case 3 -- first send code to other mpi processes
        ngmpi.SendCommand('ngs_py '+source)
        # then run it on master
        code.InteractiveConsole.runcode(self, compiled_code)

        # Avoid the prompt to show up before other processes' output
        self.Barrier()
        return False

    def interact(self):
        import sys
        self.write("MPI Shell\n")
        self.write("================\n")
        self.write("Use pprint(str) to print with MPI ranks\n\n")
        self.runsource("from ngmpi import *\n")
        self.runsource("pprint=lambda p='':print('Process ' + str(ngmpi.Rank()) + '\\n'+str(p)+'\\n')\n")
        self.locals.update(locals())
        self.locals.update(globals())
        sys.ps1 = "MPI >>> "
        sys.ps2 = "MPI ... "
        code.InteractiveConsole.interact(self,'')
        sys.ps1 = ">>> "
        sys.ps2 = "... "
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
    vars2 = globals()
    vars2.update(locals())
    mshell = InteractiveMPIConsole(vars2)
    mshell.interact()

def startConsole(vars2):
    try:
        import readline
        import rlcompleter
        if 'libedit' in readline.__doc__:
            readline.parse_and_bind("bind ^I rl_complete")
        else:
            readline.parse_and_bind("tab: complete")
    except:
        try:
            import pyreadline as readline
            import rlcompleter
            readline.parse_and_bind("tab:complete") # autocomplete
        except:
            print('readline not found')
    shell = code.InteractiveConsole(vars2)
    shell.push('from netgen import *')
    shell.push('from ngsolve import *')
    shell.interact(banner="NGS-Python console is up ...")
