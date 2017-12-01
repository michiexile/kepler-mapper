import dionysus
import scipy as sp
import numpy as np
from numpy.random import uniform

class Homology:
    def __init__(self, maxdim=2):
        self.maxdim = maxdim

    def barcode(self, mapper):
        closure = dionysus.closure([dionysus.Simplex(s) for s in mapper['simplices']], self.maxdim)
        cpx = dionysus.Filtration(closure)
        persistence = dionysus.homology_persistence(cpx)
        diagram = dionysus.init_diagrams(persistence, cpx)
        return [[d for d in dgm] for dgm in diagram[:self.maxdim]]

    def cycles(self, mapper):
        closure = dionysus.closure([dionysus.Simplex(s) for s in mapper['simplices']], self.maxdim)
        cpx = dionysus.Filtration(closure)

        # Add cone point to force homology to finite length; Dionysus only gives out cycles of finite intervals
        spxs = [dionysus.Simplex([-1])] + [c.join(-1) for c in cpx]
        for spx in spxs:
            spx.data = 1
            cpx.append(spx)

        persistence = dionysus.homology_persistence(cpx)
        diagram = dionysus.init_diagrams(persistence, cpx)
        representatives = [[persistence[persistence.pair(interval.data)]
                            for interval in dgm
                            if persistence.pair(interval.data) != persistence.unpaired]
                           for dgm in diagram[:self.maxdim]]
        naked_chains = [[{frozenset(cpx[s.index]): s.element for s in rep} for rep in dim] for dim in representatives]
        return naked_chains

class Persistence:
    def __init__(self, maxdim=2, maxeps=None):
        self.maxdim = maxdim
        self.maxeps = maxeps

    def barcode(self, points):
        if self.maxeps is None:
            eps = points.std(axis=0).max()
        else:
            eps = self.maxeps
        cpx = dionysus.fill_rips(points, self.maxdim, self.maxeps)
        persistence = dionysus.homology_persistence(cpx)
        diagram = dionysus.init_diagrams(persistence, cpx)
        return [[d for d in dgm] for dgm in diagram[:self.maxdim]]


    def cycles(self, points, top=None):
        if self.maxeps is None:
            eps = points.std(axis=0).max()
        else:
            eps = self.maxeps
        cpx = dionysus.fill_rips(points, self.maxdim, eps)

        # Add cone point to force homology to finite length; Dionysus only gives out cycles of finite intervals
        spxs = [dionysus.Simplex([-1])] + [c.join(-1) for c in cpx]
        for spx in spxs:
            spx.data = 10*eps
            cpx.append(spx)

        persistence = dionysus.homology_persistence(cpx)
        diagram = dionysus.init_diagrams(persistence, cpx)

        if top is None:
            representatives = [[persistence[persistence.pair(interval.data)]
                                for interval in dgm
                                if persistence.pair(interval.data) != persistence.unpaired]
                               for dgm in diagram[:self.maxdim]]
            naked_chains = [[{frozenset(cpx[s.index]): s.element for s in rep} for rep in dim] for dim in representatives]
            return naked_chains
        else:
            intervals = sorted([d for dgm in diagram[:self.maxdim] for d in dgm], key=lambda d: d.death-d.birth, reverse=True)
            representatives = [persistence[persistence.pair(interval.data)]
                                for interval in intervals[:top+1]
                                if persistence.pair(interval.data) != persistence.unpaired]
            naked_chains = [{frozenset(cpx[s.index]): s.element for s in rep} for rep in representatives]
            return naked_chains


