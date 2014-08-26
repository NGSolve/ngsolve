# load a pde, and exec in the py-console

# expression templates from python
globals().update(run_module('bla',globals()))
globals().update(run_module('la',globals()))

vecu = pde.gridfunctions["u"].Vector()
vecf = pde.linearforms["f"].Vector()
mata = pde.bilinearforms["a"].Matrix()

vecu.size
mata.height
mata.width


vecv = vecu.CreateVector()
vecw = vecu.CreateVector()

# no freedofs yet
inva = mata.Inverse()

vecv.data = inva * vecf

vecw.data = vecu - vecv

print ("diff^2 : ", InnerProduct (vecw, vecw))
