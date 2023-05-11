import numpy as np
import scipy as sci
import itertools
from collections import defaultdict
from functools import partial
from scipy.optimize import fmin, fmin_cg, fmin_powell, fmin_bfgs
from scipy.spatial import ConvexHull, Delaunay

def circum_points(center, r, num_pts=100):

    indices = np.arange(0, num_pts, dtype=float) + 0.5
    theta = np.pi * (1 + 5**0.5) * indices
    phi = np.arccos(1 - 2*indices/num_pts)

    if len(center) == 2:
        x, y = r * np.cos(theta) + center[0], r * np.sin(theta) + center[1]
        circle = np.column_stack((x,y))
        return circle
    elif len(center) == 3:

        x, y, z = r * np.cos(theta) * np.sin(phi) + center[0], r * np.sin(theta) * np.sin(phi) + center[1], r * np.cos(phi) + center[2]
        sphere = np.column_stack((x,y,z))
        return sphere
    
    else:
        return None

def buffer_hull(vertices, r=None, num_pts=10):

    if r is None:
        center = np.mean(vertices, axis=0)
        r = np.amax(np.linalg.norm(vertices-center, axis=1))
        #print("r=" + str(r) + " is set for the buffer hull")

    buffer_points = []
    for vertice in vertices:
        buffer_points.append(circum_points(vertice, r, num_pts))
    buffer_points = np.vstack(buffer_points)

    return ConvexHull(buffer_points)

def bisector_between(point1=None, point2=None):

    middle_pnt = (point1 + point2) / 2
    normal_vec = (point2 - point1) / np.linalg.norm(point2 - point1)

    if normal_vec.dot(point1.T) <= normal_vec.dot(middle_pnt):
        A = normal_vec
        b = normal_vec.dot(middle_pnt)

    else:
        A = -normal_vec
        b = -normal_vec.dot(middle_pnt)
    
    return A, np.array([b])

def intersection_of(vertices1=None, vertices2=None, tol=1e-10):

    if any(in_hull(vertices1, vertices2, tol=tol)) or any(in_hull(vertices2, vertices1, tol=tol)):
        A1, b1, Aeqn_1, beqn_1 = V2Ab(vertices1, tol=tol)
        A2, b2, Aeqn_2, beqn_2 = V2Ab(vertices2, tol=tol)

        Ab1 = np.concatenate((A1,b1), axis=1)
        Ab2 = np.concatenate((A2,b2), axis=1)

        _, indices = unique_rows(np.concatenate((Ab1, -Ab2)), tol=tol, return_index=True)
        if False in indices:
            V = np.array([])

        else:
            A = np.concatenate((A1, A2))
            b = np.concatenate((b1, b2))

            V = Ab2V(A, b, tol=tol)

            if V.size == 0 or False in in_hull(np.mean(V, axis=0), vertices1, tol=tol) or False in in_hull(np.mean(V, axis=0), vertices2, tol=tol):
                V = np.array([])
            elif False in in_hull(V, vertices1, tol=tol) or False in in_hull(V, vertices2, tol=tol):
                A = A + np.random.normal(0, 1, A.shape) * tol
                b = b + np.random.normal(0, 1, b.shape) * tol
                V = Ab2V(A, b, tol=tol)
                print("function intersection_of got some unexpected results")
                print("please be careful of the final result, it could have some numerical noise beyond the tolerance")

    else:
        V = np.array([])

    return V

def rownormalize(A=None, b=None):

    row_norms = np.linalg.norm(A, ord=2, axis=1, keepdims=True)

    A = A / row_norms
    b = b / row_norms

    return A, b

def in_hull(points, hull_or_vertices, tol=1e-10):

    if isinstance(hull_or_vertices, ConvexHull):
        hull = hull_or_vertices
    
    else:
        hull = ConvexHull(hull_or_vertices)
    
    points = np.atleast_2d(points)
    return np.all(np.add(np.dot(points, hull.equations[:,:-1].T), hull.equations[:,-1]) < tol, axis=1)

def V2ineqnAb(V=None, tol=1e-10):

    V = np.atleast_2d(V)
    V[abs(V) < tol] = 0

    s = ConvexHull(V).simplices
    centroid = np.mean(V[np.unique(s),:], axis=0)[np.newaxis]
    V = V - centroid

    A = np.full((s.shape[0], V.shape[1]), np.nan)

    dim = V.shape[1]
    I = np.ones(s.shape[1])

    i = 0
    for j in range(len(s)):
        F = V[s[j,:],:]

        if np.linalg.matrix_rank(F, tol) == dim:

            A[i,:], _, _, _ = np.linalg.lstsq(F, I, rcond=None)
            i = i + 1

    A = A[range(i),:]
    b = np.ones((A.shape[0], 1))

    b = b + A.dot(centroid.T)

    [A, b] = rownormalize(A, b)

    Ab = np.concatenate((A,b), axis=1)
    _, indices = unique_rows(Ab, tol=tol, return_index=True)

    A = A[indices,:]
    b = b[indices]

    A[abs(A) < tol] = 0
    b[abs(b) < tol] = 0

    return A, b

def V2Ab(V=None, tol=1e-10):

    V = np.atleast_2d(V)
    V[abs(V) < tol] = 0

    p = V[0,:][np.newaxis].T
    X = V.T - p

    Q, R, P = sci.linalg.qr(X, mode='economic', pivoting=True)

    diagr = abs(np.diag(R))
    if np.count_nonzero(diagr):
        rank = np.count_nonzero(diagr > tol * diagr[0])
        Rsub = R[range(rank)][:,P.argsort()].T

        if rank > 1:
            A, b = V2ineqnAb(Rsub, tol)

        else:
            if rank == 1:
                A = np.array([[1], [-1]])
                b = np.array([max(Rsub), -min(Rsub)])

        A = A.dot(Q[:,range(rank)].T)
        b = b + A.dot(p)

        if rank < V.shape[1]:
            Aeqn = Q[:,rank:].T
            beqn = Aeqn.dot(p)
        else:
            Aeqn = np.array([])
            beqn = np.array([])
    else:
        A = np.array([])
        b = np.array([])
        Aeqn = np.eye(V.shape[1])
        beqn = p

    A[abs(A) < tol] = 0
    b[abs(b) < tol] = 0
    Aeqn[abs(Aeqn) < tol] = 0
    beqn[abs(beqn) < tol] = 0

    return A, b, Aeqn, beqn

def unique_rows(V, tol=1e-10, return_index=False):
    
    remove = np.zeros(len(V), dtype=bool)
    for i in range(len(V)):
        equals = np.all(np.isclose(V[i,:].astype(float), V[(i+1):,:].astype(float), atol=tol), axis=1)
        remove[(i+1):] = np.logical_or(remove[(i+1):], equals)
    if return_index:
        return V[np.logical_not(remove)], np.logical_not(remove)
    else:
        return V[np.logical_not(remove)]

def f(x=None, param=None, tol=1e-10):

    A = param[0]
    b = param[1]

    d = A.dot(np.atleast_2d(x).T) - b

    d[d > tol] = d[d > tol] + 1

    d = max(np.concatenate((np.atleast_2d(0), d)))

    return d

def Ab2V(A=None, b=None, tol=1e-10):

    A[abs(A) < tol] = 0
    b[abs(b) < tol] = 0

    guess_point, _, _, _ = np.linalg.lstsq(A, b, rcond=None)

    if not all(abs(A.dot(guess_point) - b) < tol):
        interior_point, fopt, _, _, warnflag = fmin_cg(lambda x=None: f(x,(A,b)), x0=guess_point, full_output=True, disp=False)
        interior_point, fopt, _, _, warnflag = fmin(lambda x=None: f(x,(A,b)), x0=interior_point, xtol=tol, ftol=tol, maxfun=1/tol, full_output=True, disp=False)
        interior_point, fopt, _, _, _, _, warnflag = fmin_bfgs(lambda x=None: f(x,(A,b)), x0=interior_point, full_output=True, disp=False)
        interior_point, fopt, _, _, _, warnflag = fmin_powell(lambda x=None: f(x,(A,b)), x0=interior_point, xtol=tol, ftol=tol, maxfun=1/tol, full_output=True, disp=False)
         
        i = 0
        max_noise = 1000
        while not fopt < tol and i < max_noise:
            i = i + 1
            A_copy = A.copy()
            b_copy = b.copy()
            A = A + np.random.normal(0, 1, A.shape) * tol
            b = b + np.random.normal(0, 1, b.shape) * tol

            interior_point = 0.5 * (np.atleast_2d(interior_point).T + guess_point)
            interior_point, fopt, _, _, warnflag = fmin_cg(lambda x=None: f(x,(A,b)), x0=interior_point, full_output=True, disp=False)
            interior_point, fopt, _, _, warnflag = fmin(lambda x=None: f(x,(A,b)), x0=interior_point, xtol=tol, ftol=tol, maxfun=1/tol, full_output=True, disp=False)
            interior_point, fopt, _, _, _, _, warnflag = fmin_bfgs(lambda x=None: f(x,(A,b)), x0=interior_point, full_output=True, disp=False)
            interior_point, fopt, _, _, _, warnflag = fmin_powell(lambda x=None: f(x,(A,b)), x0=interior_point, xtol=tol, ftol=tol, maxfun=1/tol, full_output=True, disp=False)

        if warnflag > 0 or fopt > tol:
            V = np.array([])
            print("function Ab2V failed to find a point clearly inside the region defined by halfspaces for")
            print("A=" + str(A_copy))
            print("b=" + str(b_copy))

        else:
            b = b - A.dot(np.atleast_2d(interior_point).T)
            D = A / np.tile(b,np.array([1, A.shape[1]]))
            
            try:
                s = ConvexHull(D).simplices
                G = np.zeros((s.shape[0], D.shape[1]))
            except:
                V = np.array([])
            else:
                for i in range(len(s)):
                    F = D[s[i,:],:]
                    g, _, _, _ = np.linalg.lstsq(F,np.ones((F.shape[0], 1)), rcond=None)
                    G[i,:] = g.T

                V = G + np.tile(interior_point.T,np.array([G.shape[0], 1]))
                V = unique_rows(V, tol=tol)
                V[abs(V) < tol] = 0
        
    return V

def bounded_voronoi(points=None, boundary_vertices=None, tol=1e-10):

    if False in in_hull(points, boundary_vertices, tol=tol):
        print("some points are outside the boundary")
        print("you should use function soft_bounded_voronoi instead of bounded_voronoi")
        print(in_hull(points, boundary_vertices, tol=tol))

    s = ConvexHull(boundary_vertices).simplices
    V_boundary = boundary_vertices[s,:].reshape((s.size, s.shape[1]))
    A_boundary, b_boundary, _, _ = V2Ab(V_boundary, tol=tol)
    
    tri = Delaunay(points)
    neighbors = defaultdict(set)
    for p in tri.vertices:
        for i, j in itertools.combinations(p,2):
            neighbors[i].add(j)
            neighbors[j].add(i)
    
    A_bisector = defaultdict(lambda: defaultdict(partial(np.ndarray, 0)))
    b_bisector = defaultdict(lambda: defaultdict(partial(np.ndarray, 0)))
    for i in range(len(points)):
        for j in neighbors[i]:
            A_bisector[i][j], b_bisector[i][j] = bisector_between(points[i], points[j])

    A_voronoi = defaultdict(partial(np.ndarray, 0))
    b_voronoi = defaultdict(partial(np.ndarray, 0))
    for i in range(len(points)):
        A_voronoi[i] = np.concatenate((np.stack(list(A_bisector[i].values())), A_boundary))
        b_voronoi[i] = np.concatenate((np.stack(list(b_bisector[i].values())), b_boundary))

    V_voronoi = defaultdict(partial(np.ndarray, 0))
    voronoi_cells = defaultdict(partial(ConvexHull, 0))
    for i in range(len(points)):
        V_voronoi[i] = Ab2V(A_voronoi[i], b_voronoi[i], tol=tol)
        try:
            voronoi_cells[tuple(points[i])] = ConvexHull(V_voronoi[i])
        except:
            hull = ConvexHull(V_voronoi[i], qhull_options='QJ')
            hull_vertices = hull.points[hull.vertices]
            hull_vertices = intersection_of(hull_vertices, boundary_vertices, tol=tol)
            if hull_vertices.size != 0:
                voronoi_cells[tuple(points[i])] = ConvexHull(hull_vertices)
                print("qhull_options='QJ is used for point " + str(points[i]))

    return voronoi_cells

def soft_bounded_voronoi(points=None, boundary_vertices=None, r=None, tol=1e-10):

    buffer = buffer_hull(boundary_vertices, r, num_pts=10)
    points = points[in_hull(points, buffer, tol=tol)]

    if all(in_hull(points, boundary_vertices, tol=tol)) == True:
        voronoi_cells = bounded_voronoi(points, boundary_vertices, tol=1e-10)
        return voronoi_cells

    boundary = ConvexHull(boundary_vertices)
    cells = bounded_voronoi(points, points[ConvexHull(points).vertices], tol=tol)

    V_sum = 0.0
    voronoi_cells = defaultdict(partial(ConvexHull, 0))

    for point, cell in cells.items():
        cell_vertices = cell.points[cell.vertices]
        V = intersection_of(cell_vertices, boundary_vertices, tol=tol)

        try:
            voronoi_cells[point] = ConvexHull(V)
            V_sum += ConvexHull(V).volume
        except:
            pass
    
    if not np.isclose(V_sum, boundary.volume, atol=tol):

        points = np.array(list(voronoi_cells.keys()))
        if all(in_hull(points, boundary_vertices, tol=tol)) == True:
            voronoi_cells = bounded_voronoi(points, boundary_vertices, tol=1e-10)
            return voronoi_cells

        print("The final result has some numerical noise beyond the tolerance" + str(tol))
        print("Function soft_bounded_voronoi is attempting to resolve this")

        cells = bounded_voronoi(points, points[ConvexHull(points).vertices])

        V_sum = 0.0
        voronoi_cells = defaultdict(partial(ConvexHull, 0))

        for point, cell in cells.items():
            cell_vertices = cell.points[cell.vertices]
            V = intersection_of(cell_vertices, boundary_vertices, tol=tol)

            try:
                voronoi_cells[point] = ConvexHull(V)
                V_sum += ConvexHull(V).volume
            except:
                pass
        
        if np.isclose(V_sum, boundary.volume, atol=tol):
            print("Succeeded, the numerical noise now is" + str(boundary.volume - V_sum) + "<" + str(tol))
        else:
            print("Failed, the numerical noise now is " + str(boundary.volume - V_sum) + ">" + str(tol))
            print("The final results usually are still useful, but please be careful of it and report the bug to Siyu (sc2090@cam.ac.uk)")
    
    return voronoi_cells
