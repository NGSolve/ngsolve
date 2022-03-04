import numpy as np

#returns volume fraction of negative part of triangle
def GetVolumeFraction(psi0, psi1, psi2, EPS):

    if psi0 < -EPS and psi1 < -EPS and psi2 < -EPS:
        return 1
    elif psi0 > EPS and psi1 > EPS and psi2 > EPS:
        return 0
    
    ### if p0 is cut...
    if abs(psi0)<EPS and abs(psi1)>=EPS and abs(psi2)>=EPS: #psi0 is on vertex, psi1 and psi2 not
        if np.sign(psi1)!=np.sign(psi2):        # edge p1p2 is cut
            s = 1 * (-psi1)/(psi2 - psi1)
            if psi1 > 0:
                s = 1-s
        else:   #element is only cut in p0 (tangent)
            s = -0.5 * np.sign(psi1) + 0.5      # if psi1 <0 (-->psi2<0) then s=1, else (i.e. psi1, psi2>0) then s=0)
        return s
    elif abs(psi0)<EPS and abs(psi1)<EPS and abs(psi2)>=EPS: #edge p0p1 is on interface
        if psi2>0:
            s=0
        else:
            s=1
        return s
    elif abs(psi0)<EPS and abs(psi1)>=EPS and abs(psi2)<EPS: #edge p0p2 is on interface
        if psi1>0:
            s=0
        else:
            s=1
        return s
            
    ### if p1 is cut...
    if abs(psi1)<EPS and abs(psi2)>=EPS and abs(psi0)>=EPS: #psi1 is on vertex, psi2 and psi0 not
        if np.sign(psi2)!=np.sign(psi0):        # edge p2p0 is cut
            s = 1 * (-psi2)/(psi0 - psi2)
            if psi2 > 0:
                s = 1-s
        else:   #element is only cut in p1 (tangent)
            s = -0.5 * np.sign(psi2) + 0.5      # if psi2 <0 (-->psi0<0) then s=1, else (i.e. psi2, psi0>0) then s=0)
        return s
    elif abs(psi1)<EPS and abs(psi2)<EPS and abs(psi0)>=EPS: #edge p1p2 is on interface
        if psi0>0:
            s=0
        else:
            s=1
        return s
    elif abs(psi1)<EPS and abs(psi2)>=EPS and abs(psi0)<EPS: #edge p1p0 is on interface
        if psi2>0:
            s=0
        else:
            s=1
        return s
            
    ### if p2 is cut...
    if abs(psi2)<EPS and abs(psi0)>=EPS and abs(psi1)>=EPS: #psi2 is on vertex, psi0 and psi1 not
        if np.sign(psi0)!=np.sign(psi1):        # edge p0p1 is cut
            s = 1 * (-psi0)/(psi1 - psi0)
            if psi0 > 0:
                s = 1-s
        else:   #element is only cut in p1 (tangent)
            s = -0.5 * np.sign(psi0) + 0.5      # if psi0 <0 (-->psi1<0) then s=1, else (i.e. psi0, psi1>0) then s=0)
        return s
    elif abs(psi2)<EPS and abs(psi0)<EPS and abs(psi1)>=EPS: #edge p2p0 is on interface
        if psi1>0:
            s=0
        else:
            s=1
        return s
    elif abs(psi2)<EPS and abs(psi0)>=EPS and abs(psi1)<EPS: #edge p1p0 is on interface
        if psi0>0:
            s=0
        else:
            s=1
        return s
         
    if abs(psi0)>=EPS and abs(psi1)>=EPS and abs(psi2)>=EPS:    
        if np.sign(psi1) == np.sign(psi2) and np.sign(psi0) != np.sign(psi1):
            s = psi0*psi0 / ((psi1 - psi0)*(psi2 - psi0))
            if psi0 > 0:
                s = 1-s
        elif np.sign(psi2) == np.sign(psi0) and np.sign(psi1) != np.sign(psi2):
            s = psi1*psi1 / ( (psi2-psi1)*(psi0-psi1) )
            if psi1 > 0:
                s = 1-s
        elif np.sign(psi0) == np.sign(psi1) and np.sign(psi2) != np.sign(psi0):
            s = psi2*psi2 / ( (psi0-psi2)*(psi1-psi2) )
            if psi2 > 0:
                s = 1-s
        else:
            print(psi0, psi1, psi2)
            print("Error. This should not happen.")
    else:       #all three values are below EPS
        print(psi0, psi1, psi2)
        print("Error. This shoul not happen. Think about re-initializing your level set function as a signed distance function.")
    return s
        
        
#val1 where levelset_p1 < 0; val2 where levelset_p1 > 0
def InterpolateLevelSetToElems(levelset_p1, val1, val2, func_p0, mesh, EPS ):
    s  = 0
    nCutElems = 0
    for i,el in enumerate(mesh.ngmesh.Elements2D()):

        psi0 = levelset_p1.vec[el.vertices[0].nr-1]
        psi1 = levelset_p1.vec[el.vertices[1].nr-1]
        psi2 = levelset_p1.vec[el.vertices[2].nr-1]
        
        if psi0 < 0 and psi1 < 0 and psi2 < 0:
            func_p0.vec[i] = val1
        elif psi0>0 and psi1>0 and psi2 > 0:
            func_p0.vec[i] = val2
        else:
            s = GetVolumeFraction(psi0, psi1, psi2, EPS)  #s... volume fraction of negative part of triangle
            nCutElems += 1
            func_p0.vec[i] = val2 + s*(val1-val2) #== s*val1 + (1-s)*val2
