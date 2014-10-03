# load a pde, and exec in the py-console


pde.Solve()

vecu = pde.gridfunctions["u"].Vector()
vecf = pde.linearforms["f"].Vector()
mata = pde.bilinearforms["a"].Matrix()

v = pde.spaces["v"]

vecu.size
mata.height
mata.width

vecv = vecu.CreateVector()
vecw = vecu.CreateVector()

# atn: don't use eliminate_internal
inva = mata.Inverse(v.FreeDofs())

vecv.data = inva * vecf
vecw.data = vecu - vecv

print ("diff^2 : ", la.InnerProduct (vecw, vecw))

