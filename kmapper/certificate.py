import scipy as sp
import numpy as np
from numpy.random import multivariate_normal, uniform
from sklearn import decomposition
from sklearn import cluster
from collections import defaultdict
import networkx as nx

from .persistence import Persistence
from .nerve import GraphNerve

def second_interval(barcode):
    return sorted([c - b for (a, (b, c)) in barcode])[-2]

def simulation_test(data, N=100, invariant=second_interval):
    persistence = Persistence(maxdim = data.shape[1], maxeps = 4*max(data.std(axis=0)))
    simulationP = [invariant([(i, (interval.birth,interval.death)) for (i, dgm) in enumerate(persistence.barcode(data)) for interval in dgm])]
    for simulation in [
                multivariate_normal(mean=data.mean(axis=0), cov=np.cov(data.T), size=data.shape[0])
                for _ in range(N - 1)]:
        simulationP.append(invariant([(i, (interval.birth,interval.death)) for (i, dgm) in enumerate(persistence.barcode(simulation)) for interval in dgm]))
    return 1 - np.where(np.argsort(simulationP) == 0)[0] / N


def certificate(graph, data, minsize=4, significance = 0.05, simulations=40, factor=0.5, nerve=GraphNerve()):
    """
    Will certify goodness of a cover for a Mapper complex, and if the cover is not good will suggest a refinement.
    """
    simplices = graph['simplices'].copy()
    nodes = graph['nodes'].copy()
    links = graph['links'].copy()
    splits = {n: None for n in nodes.keys()}
    queue = sorted(simplices, key=len, reverse=True)
    done = []
    persistence = Persistence(maxdim=data.shape[1], maxeps=max(data.std(axis=0)))
    while len(queue) > 0:
        simplex = queue.pop(0)
        if any([splits[s] is not None for s in simplex]):
            # This simplex has already been split up and handled; continue
            continue

        idxs = set(nodes[simplex[0]]).intersection(*[nodes[s] for s in simplex[1:]])
        if len(idxs) < minsize: continue
        points = data[sorted(idxs), :]
        if simulation_test(points, N = simulations) < significance:
            # time to split something!
            cycles = persistence.cycles(points, top=5)
            vertices = [(len(k[0]), sorted(k[0].union(*k[1:]))) for n in cycles for k in [list(n.keys())]]
            dim, largest_cycle = sorted(vertices, key=lambda dk: len(dk[1]))[-1]
            kmeans = cluster.KMeans(n_clusters=dim+1)
            kmeans.fit(data[largest_cycle, :])
            changes = simplex.copy()
            for coverset in changes:
                # split each top-level set
                barycentric = kmeans.transform(data[nodes[coverset]])
                membership = ((1+factor)*barycentric).T - factor*(barycentric.sum(axis=1))
                membership = membership.T > 0

                ## next up: actually perform the split! we now have required cluster memberships
                # Add references in splits
                splits[coverset] = [coverset+"_split{}".format(i) for i in range(kmeans.n_clusters)]
                for label in splits[coverset]: splits[label] = None

                # Adjust the nodes collection
                newnodes = []
                for i, label in enumerate(splits[coverset]):
                    nodes[label] = [nodes[coverset][j] for j in range(membership.shape[0]) if membership[j,i]]
                    newnodes.append(label)

                # Adjust the simplices and links
                newnodes.extend([[t for t in spx if t != coverset] for spx in simplices if coverset in spx])
                newnodes.remove([])

                newlinks, newsimplices = nerve({n: nodes[n] for n in newnodes})

                for spx in simplices:
                    if coverset in spx:
                        simplices.remove(spx)
                simplices.extend([s for s in newsimplices if s not in simplices])
                if coverset in links:
                    del(links[coverset])
                for k in links:
                    if coverset in links[k]:
                        links[k].remove(coverset)
                for k in newlinks:
                    links[k].extend([s for s in newlinks[k] if s not in links[k]])

                # Adjust the queue
                for spx in queue:
                    if coverset in spx:
                        queue.remove(spx)
                for spx in newsimplices:
                    if spx not in queue and spx not in done:
                        queue.append(spx)
        else:
            done.append(simplex)





def code():
    # check whether all clusters and all cluster intersections are acyclic
    # [using a persistent homology software library... gudhi!]
    # gather up all non-acyclic cases in a list of tasks to handle later

    # go from deep to shallow intersections: a deep intersection having homology will carry out to shallow ones
    splits = []
    for (id1, id2) in network.edges():
        data = list(set(id1).intersection(id2))
        if len(data) <= 2:
            continue
        if simulation_test(X[data, :]) < .05:
            splits.append((id1, id2))

    for data in network.nodes():
        if len(data) <= 1:
            continue
        if simulation_test(X[data, :]) < .05:
            splits.append(data)

    # for each non-acyclic cluster, split the deepest acyclic intersection into several pieces.
    # Then use some nice heuristic to propagate the split out to clusters
    #
    # * local PCA to create a global additional filter function?
    # * local PCA to create a local split?
    # * local k-means to create a local split?
    # * use representative to ensure split? [gudhi doesn't provide chains]
    #
    # For propagation:
    # * use projection onto local PCA coordinates to split things
    # * use local centroids to create voronoi-partition of supersets
    while splits:
        print("Should be doing stuff...")
        print(splits)
        splits = []
        # step 1: use `get_cycles` below to extract cycle supports
        # step 2: for a cycle in dimension j, do (j+2)-means clustering
        # step 3: soft nearest neighbor overlapping splits on each containing
        #         cover element

    # repeat until everything is sufficiently acyclic

    return network

def whiten(data):
    return (data - data.mean(axis=0)) / data.std(axis=0)


