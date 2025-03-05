

import spherogram
from spherogram import links
import numpy as np
import sympy
import ast

from sympy import symbols, Matrix 
from datetime import datetime


from scipy.optimize import least_squares
from functools import lru_cache


# Generate timestamp
timestamp = datetime.now().strftime("%Y%m%d_%H%M")

# Define the log file name with timestamp
log_file = f"log_{timestamp}.txt"

# Save the log file in the current folder
with open(log_file, "a", encoding="utf-8") as f:
    f.write("Log file created with timestamp.")




######################### Equations from PSL representation paper #########################


# Save printed results
def print1(message=' '):
    #print(message)
    print(message)
    with open(log_file, "a") as file:
        file.write(f"{message}\n")     

# Save printed results
def print2(message=' '):
    #print(message)
    print((message), end="")

    with open(log_file, "a") as file:
        file.write(f"{message}\n")     



BLACK = "Black"
WHITE = "White"

OPPOSITE = {WHITE: BLACK, BLACK: WHITE}


# MUST CHECK if edge goes left to right as viewed by region or right to left. Currently assumes left to right. 
# This determines g. If a strand moves from an underpass to an over pass, g = -1. If over -> under, g=+1. If they are paralell, g=0. 

def g(edge):
    if edge.strand_index in [1, 3] and edge.opposite().strand_index in [0, 2]:
        return -1
    elif edge.strand_index in [0, 2] and edge.opposite().strand_index in [1, 3]:
        return 1
    else:
        return 0



def find_next_edge(faces, current_edge):
    for face in faces:
        if current_edge in face:
            current_index = face.index(current_edge)
            next_index = (current_index + 1) % len(face)  # Wrap around for circular navigation
            return face[next_index]
    return None



def solve(equations, crossing_vars, edge_vars, meridian_var, crossing_labels, edge_labels, meridian_label, num_iterations=30):
    """
    Find a real-valued solution to the given system of equations.

    Uses Newton--Raphson iteration to solve the system.
    """
    # Create lambdified versions of the equations and their Jacobian matrix
    f_temp = sympy.lambdify(crossing_vars + edge_vars + meridian_var, equations, "numpy")
    Jac = sympy.Matrix(equations).jacobian(crossing_vars + edge_vars + meridian_var)
    Df_temp = sympy.lambdify(crossing_vars + edge_vars + meridian_var, Jac, "numpy")

    # Function to compute f(x)
    f = lambda x: f_temp(*x)

    # Function to compute Jacobian Df(x)
    Df = lambda x: Df_temp(*x)

    # Initial guess (convert to real numbers explicitly)
    labels = np.array(crossing_labels + edge_labels + meridian_label, dtype=np.float64)

    # Newton--Raphson iteration
    for _ in range(num_iterations):
        f_x = f(labels)
        Df_x = Df(labels)

        # QR decomposition for numerical stability
        q, r = np.linalg.qr(Df_x)
        incx = np.linalg.solve(r, -np.matmul(q.T, f_x))  # Restrict to real arithmetic
        labels = np.around(labels + incx, 6)

        # If f_x is sufficiently close to zero, exit early
        if np.linalg.norm(f_x, ord=2) < 1e-6:
            break

    # Extract solutions for crossing, edge, and meridian variables
    crossing_solutions = labels[: len(crossing_labels)]
    edge_solutions = labels[len(crossing_labels) : len(crossing_labels) + len(edge_labels)]
    meridian_solution = [labels[-1]]

    return crossing_solutions, edge_solutions, meridian_solution


def analyze_PSL(L):
    # Collect edges.  Each edge is a (CrossingStrand, CrossingEntryPoint) pair, where the latter is the 'head' of the oriented edge.
    # Note that comparing a (CrossingStrand, CrossingEntryPoint) pair to a (CrossingStrand, CrossingStrand) pair (as are given in the faces) still works properly since CrossingEntryPoint inherits from CrossingStrand.
    
    faces = [[x.opposite() for x in face] for face in L.faces()]  # This orientation matters since it defines clockwise in each region.
    
    
    # Edge variables and starting labels. 
    " 1.First checkerboard color all the faces."
    coloring = [None] * len(faces)
    coloring[0] = BLACK
    add_index = [0]  # A stack for DFT.
    while add_index:
        face_i = add_index.pop()
        for nbr_i, nbr in enumerate(faces):
            if coloring[nbr_i] is None and any(e.opposite() in faces[face_i] for e in nbr):
                coloring[nbr_i] = OPPOSITE[coloring[face_i]]
                add_index.append(nbr_i)


    "2.Next we iterate through all edges by iterating through black faces and assigning indices and names to the corresponding edge u_i and the reversed edge v_i (bordering a white region)" 
    black_faces = [face for face, color in zip(faces, coloring) if color == BLACK]
    edge_indices = {side: index                     # Create a dictionay that assigns a unique index to every edge and its opposite.( u_i and v_i have the same index i )
                    for index, side in enumerate(
                        oedge 
                        for face in black_faces 
                        for edge in face 
                        for oedge in [edge, edge.opposite()]
                        )
                    }
    
    equations = []


    # Create a local name for the set of crossings for consistency
    crossings = L.crossings

    # Meridian matrix
    m = symbols('m')
    M = Matrix([[m, 1], [0 , 1/m]]) 

    # Set meridian element initial value
    meridian_var = [m]
    


    # Create matrix labels and variables of matrices. 
    edge_matrix_label = [f"{prefix}{i+1}" for i in range(2 * len(crossings)) for prefix in ["P_l", "P_r"]] # P_r1, P_l1, P_r2, P_l2, ..... len(edges) = 2*len(crossings) 
    edge_element_vars = sympy.symbols([f"{prefix}{i+1}" for i in range(2 * len(crossings)) for prefix in ["U_l", "V_l", "U_r", "V_r"]]) # U_l1, V_l1, U_r1, V_r1 ....
    print1(edge_element_vars)

    # Create edge matrices P_l & P_r
    edge_matrices = {} 
    for face in black_faces:
        for edge in face:
            i = edge_indices[edge] 
            i_opposite = edge_indices[edge.opposite()]

            edge_matrices[edge_matrix_label[i]] = Matrix([      # P_l
            [edge_element_vars[2*i+1], edge_element_vars[2*i]], 
            [0, 1/edge_element_vars[2*i+1]]
            ])

            edge_matrices[edge_matrix_label[i_opposite]] = Matrix([      # P_r
            [edge_element_vars[2*i_opposite+1], edge_element_vars[2*i_opposite]], 
            [0, 1/edge_element_vars[2*i_opposite]]
            ])
            

   

    # Crossing variables & crossing matirx labels
    crossing_vars = sympy.symbols([f"w{i+1}" for i in range(len(crossings))]) # w1, w2, ....
    crossing_matrix_label = [f"{prefix}{i+1}" for i in range(len(crossings)) for prefix in ["W"]] # W1, W2, ... 
    print1(crossing_vars)

    
    

   
    # Create crossing matrices
    crossing_matrices = {}
    for i in range(len(crossings)):
        crossing_matrices[crossing_matrix_label[i]] = Matrix([[0, crossing_vars[i]], [1, 0]])


    
    " 3. Go through all faces and create region equations P_1 * W_1 ..... = xI "
    print1()
    print1("@ Region equations")
    for face_i, face in enumerate(faces):
        if len(face) == 1:
            raise RuntimeError("Monogon")
        
        elif len(face) == 2:
            print1(f"Bigon case: Face number {face_i}: {face}")
            e1, e2 = face
            i_b1 = edge_indices[e1]
            i_b2 = edge_indices[e2]
            
            equations.append(edge_element_vars[2*i_b1]) # U_l1 = 0
            exp = edge_element_vars[2*i_b1]
            print2(exp)
            print1('= 0')

            equations.append(edge_element_vars[2*i_b2]) # U_l2 = 0
            exp = edge_element_vars[2*i_b2]
            print2(exp)
            print1('= 0')

            equations.append(edge_element_vars[2*i_b1 + 1] - 1/edge_element_vars[2*i_b1 + 1] ) # V_l1 - 1 / V_l1
            exp = (edge_element_vars[2*i_b1 + 1] - 1/edge_element_vars[2*i_b1 + 1])
            print2(exp)
            print1('= 0')

            equations.append(edge_element_vars[2*i_b2 + 1] - 1/edge_element_vars[2*i_b2 + 1] )  # V_l2 - 1 / V_l2
            exp = (edge_element_vars[2*i_b2 + 1] - 1/edge_element_vars[2*i_b2 + 1])
            print2(exp)
            print1('= 0')

            equations.append(crossing_vars[e1[0].label] - crossing_vars[e2[0].label])
            exp = (crossing_vars[e1[0].label] - crossing_vars[e2[0].label])
            print2(exp)
            print1('= 0')
        
        else:  # len(face) >= 2:  # Generate equations for faces with more than 1 side.
            print1()
            print1((f"Face number {face_i}: {face}"))
            Region_matrix_eq = Matrix([[1, 0], [0, 1]])
            k= 1 # Pr
            if coloring[face_i] == BLACK: 
                k = 0 # Pl
            
            for edge in face:
                next_edge = find_next_edge(faces, edge)
                Region_matrix_eq *= edge_matrices[edge_matrix_label[edge_indices[edge]]] # *Pl or *Pr
                Region_matrix_eq *= crossing_matrices[crossing_matrix_label[next_edge[0].label]]  # *W // next_edge[0].label is the label of the crossing at the start point of the next edge
            
            equations.append(Region_matrix_eq[0, 0] - Region_matrix_eq[1, 1]) 
            equations.append(Region_matrix_eq[0, 1])
            equations.append(Region_matrix_eq[1, 0])    

            equation_region_temp = []
            equation_region_temp.append(Region_matrix_eq[0, 0] - Region_matrix_eq[1, 1])
            equation_region_temp.append(Region_matrix_eq[0, 1]) 
            equation_region_temp.append(Region_matrix_eq[1, 0])
            for eq in equation_region_temp:
                print2(eq)
                print1('= 0')



    " 1. Create commute equations PM = MP <=> (m-1/m)u_i = (v_i - 1/v_i) "
    print1()
    print1("@ Commute equations")
    equation_temp=[]
    for i in range(len(edge_element_vars)//4) :
        equation_temp.append( (m - 1/ m)* edge_element_vars[4*i] - ( edge_element_vars[4*i+1] - 1/edge_element_vars[4*i+1] )  ) # U_li = edge_element_vars[4*i] , v_li = edge_element_vars[4*i+1]
        exp = (m - 1/ m)* edge_element_vars[4*i] - ( edge_element_vars[4*i+1] - 1/edge_element_vars[4*i+1])
        print2(exp)
        print1('= 0')
        equation_temp.append( (m - 1/ m)* edge_element_vars[4*i+2] - ( edge_element_vars[4*i+3] - 1/edge_element_vars[4*i+3] )  ) # U_ri = edge_element_vars[4*i+2] , v_ri = edge_element_vars[4*i+3]
        exp = (m - 1/ m)* edge_element_vars[4*i+2] - ( edge_element_vars[4*i+3] - 1/edge_element_vars[4*i+3])
        print2(exp)
        print1('= 0')
    equations.extend(equation_temp)



    " 2. Go through each black face and create the edge matrix equations P_l = M^g P_r "     
    print1()
    print1("@ Edge equations")  
    for face in black_faces:
        for edge in face:
            i = edge_indices[edge]
            i_op = edge_indices[edge.opposite()] 

            # Create edge equations
            equations.append(edge_element_vars[2*i+1] - 1/ edge_element_vars[2*i_op+1] - m ** g(edge)) # Vl - 1/Vr = m^g
            exp = edge_element_vars[2*i+1] * 1/ edge_element_vars[2*i_op+1] - m ** g(edge)
            print2(exp)
            print1('= 0')

            equations.append(edge_element_vars[2*i] * edge_element_vars[2*i_op+1] - edge_element_vars[2*i_op] * edge_element_vars[2*i+1]   - g(edge)) # Ul * Vr - Ur * Vl = g
            exp = edge_element_vars[2*i] * edge_element_vars[2*i_op+1] - edge_element_vars[2*i_op] * edge_element_vars[2*i+1]   - g(edge) 
            print2(exp)
            print1('= 0')


    

    
    
    " 1. Set meridian initial value " 
    meridian_initial = [ 0.5 ] 

    " 2. Set crossing initial values "
    crossing_initial = [ 1  if crossing.sign == 1 else  -1 for crossing in crossings]

    " 3. Set edge initinal values "
    edge_initial = [0] * 8 * len(crossings)
    for face in black_faces:
        for edge in face:
            i = edge_indices[edge]
            i_opp = edge_indices[edge.opposite()]
            edge_initial[2*i] = 0.5 
            edge_initial[2*i + 1] = 0.5  # initial edge element values for black side
            edge_initial[2*i_opp] = -0.5 
            edge_initial[2*i_opp + 1] = -0.5  # initial edge element values for white side 


    # Solve equations
    crossing_solutions, edge_solutions, meridian_solution = solve(equations, crossing_vars, edge_element_vars, meridian_var, crossing_initial, edge_initial, meridian_initial)
    
    

  
    # Output time
    print1(("Faces:"))
    for face, color in zip(faces, coloring):
        print1((f"{color} face: " + " -> ".join(str(edge[0][0]) for edge in face) + " ->"))

    print1(("'Ul, Vl' labels correspond to the black sides of the edges."))
    print1(("'Ur, Vr' labels to the white sides."))
    print1(("'w' labels correspond crossings."))

    print1(("Calculated values:"))
    print1(("1.Meridian value"))
    for m in meridian_solution:
        print1((f"m = {m}"))

    print1(("2.Crossing values"))
    for index, value in enumerate(crossing_solutions):
        print1((f"w{index + 1} = {value}"))

    print1(("3.Edge values"))
    for face in black_faces:
        print1((f"Black face: {face}"))
        for edge in face:
            i = edge_indices[edge]
            i_op = edge_indices[edge.opposite()]
            print1((f"Ul{i//2 + 1} = {edge_solutions[i * 2]}"))
            print1((f"Vl{i//2 + 1} = {edge_solutions[i * 2 + 1]}"))
            print1((f"Ur{i_op//2 + 1} = {edge_solutions[i_op * 2 ]}"))
            print1((f"Vr{i_op//2 + 1} = {edge_solutions[i_op * 2 + 1 ]}"))
            print1((f"\t(From edge {edge[0]} -> {edge.opposite()[0]})"))

   
        
######################### Equations and solutions from alternative approach to hyperbolic structure paper #########################
       







def reverse_edge_D(edge):
    return (edge[0].opposite(),edge[1].opposite())

def opposite_color_D(color):
    if color == 'white':
        return 'black'
    if color == 'black':
        return 'white'

def pars_D(code):
    "Given an PD code as a string, return a PD code as a list of tuples"
    return ast.literal_eval(code)

def get_adj_faces_D(face_i,faces):
    "ALEX NOTES: This finds adjacent faces with the presumption that they share an edge of opposite orientation."
    """
    Args:
        face_i -- index in faces of the face to find adjacent faces to
        faces -- list of all faces for the link
    """
    return [f for f in faces if any(reverse_edge_D(e) in faces[face_i] for e in f)]

def compute_k_D(edge):
    #MUST CHECK if edge goes left to right as viewed by region or right to left.
    #currently assumes left to right
    if(edge[0].strand_index in [1,3] and edge[1].strand_index in [0,2]):
        return -1
    elif(edge[0].strand_index in [0,2] and edge[1].strand_index in [1,3]):
        return 1
    else:
        return 0

def compute_kappa_D(edge_1,edge_2,oriented_edges):
    "ALEX NOTES: oriented edges is the set of edges positively oriented. This is just an XNOR for existance in oriented_edges."
    if edge_1 in oriented_edges and edge_2 in oriented_edges:
        return 1
    if edge_1 not in oriented_edges and edge_2 not in oriented_edges:
        return 1
    return -1

@lru_cache()
def get_f_n_D(shape_params,n):
    if n < 3:
        raise ValueError("parameter n must be at least 3")
    elif n == 3:
        return [1 - s for s in shape_params] # this is assigning the first entry in shape_params to correspond to zeta_2 in the notation of the paper
    elif n == 4:
        return [1 - shape_params[i] - shape_params[(i+1) % len(shape_params)] for i in range(len(shape_params))]
    else:
        f_twoback = get_f_n_D(shape_params,n-2)
        f_oneback = get_f_n_D(shape_params,n-1)
        shape_n = shape_params[2:] + shape_params[:2]
        return [f_oneback[i] - shape_n[i]*f_twoback[i] for i in range(len(shape_params))]

def analyze_alt_eq(L):
    # Collect edges.  Each edge is a (CrossingStrand, CrossingEntryPoint) pair, where the latter is the 'head' of the oriented edge.
    # Note that comparing a (CrossingStrand, CrossingEntryPoint) pair to a (CrossingStrand, CrossingStrand) pair (as are given in the faces) still works properly since CrossingEntryPoint inherits from CrossingStrand.
    oriented_edges = [(x.opposite(),x) for x in L.crossing_entries()]
    # every edge appears once as-is and once reversed in the faces
    faces = [[(x.opposite(),x) for x in face] for face in L.faces()] #hopefully right orientation
    # We create a local name for the set of crossings for consistency
    crossings = L.crossings

    # First color all the regions.  The colors are referenced by the list 'coloring' which stores the colors as strings "black" or "white"
    coloring = ['']*len(faces)
    coloring[0] = 'black'
    to_add = [faces.index(adj) for adj in get_adj_faces_D(0,faces)]
    for to_add_face in to_add:
        coloring[to_add_face] = 'white'
    while(to_add):
        face_i = to_add.pop()
        for nbr in get_adj_faces_D(face_i,faces):
            nbr_i = faces.index(nbr)
            if not coloring[nbr_i]:
                coloring[nbr_i] = opposite_color_D(coloring[face_i])
                to_add.append(nbr_i)

    # Next we iterate through all edges by iterating through black faces and assigning indices and names to the corresponding edge
    # At the same time, we assign indices and names to the reversed edge (bordering a white region)

    edge_indices = {} # used to index the edge in arrays or lists (such as 'edge_names' below)
    edge_names = [] # human readable names

    count = 0
    for face_i,face in enumerate(faces):
        if coloring[face_i] == 'white':
            continue
        for edge in face:
            edge_indices[edge] = count*2
            edge_indices[reverse_edge_D(edge)] = count*2+1
            edge_names.append('u' + str(count+1))
            edge_names.append('v' + str(count+1))
            count += 1

    #ALEX NOTES: constructor then zeroing out a sq matrix of complex numbers."
    labels = np.zeros(len(crossings) + 2*len(oriented_edges), dtype=np.complex128) #initial value estimates for all variables
    crossing_labels = labels[:len(crossings)] # reference crossing labels by slice
    edge_labels = labels[len(crossings):] # reference edge labels by slice
    edge_vars = sympy.symbols(edge_names)
    crossing_vars = sympy.symbols(['w' + str(k) for k in range(len(crossing_labels))])

    #labels = np.zeros(len(crossings+2*len(oriented_edges)), dtype=np.complex_)
    equations = []
    

    for crossing in L.crossings: # initial crossing label values
        if crossing.sign == 1: #-0.5
            crossing_labels[crossing.label] = .5+.5j  # .5j
        else:
            crossing_labels[crossing.label] = -.5-.5j  # -.5j
    # Now lets go through each face and create the equations
    count=0
    for face_i,face in enumerate(faces):
        if len(face) < 2:
            return
        if coloring[face_i] == 'black': # Add edge label equations (for the two sides of each edge) only if the color is black
            for edge in face:
               edge_labels[edge_indices[edge]] = -.5-.5j # initial edge label value for black side
               edge_labels[edge_indices[reverse_edge_D(edge)]] = .5-.5j #initial edge label value for white side
               k = compute_k_D(edge) #Set this to +-1 appropriately as per paper
               equations.append(edge_vars[edge_indices[edge]]-edge_vars[edge_indices[reverse_edge_D(edge)]]-k)

        if len(face) > 2: # generate equations for faces with more than 2 sides
            shape_params = []
            kl = compute_kappa_D(face[0],face[-1],oriented_edges)
            el = edge_vars[edge_indices[face[-1]]]
            e1 = edge_vars[edge_indices[face[0]]]
            w1 = crossing_vars[face[0][0][0].label]
            shape_par = kl*w1/el/e1
            shape_params.append(shape_par)
            for i in range(len(face) - 1):
                e_i = edge_vars[edge_indices[face[i]]]
                e_ii = edge_vars[edge_indices[face[i+1]]]
                w_i = crossing_vars[face[i+1][0][0].label]
                k_i = compute_kappa_D(face[i],face[i+1],oriented_edges)
                shape_par = k_i*w_i/e_i/e_ii
                shape_params.append(shape_par)
            edge_prod = 1
            for e in face:
                edge_prod *= edge_vars[edge_indices[e]]
            ##one=[]  #code error
            ##for i in range(len(face)):
            ##    one.append(str(edge_vars[edge_indices[face[i]]]))
            f_n= get_f_n_D(tuple(shape_params),len(shape_params)) # get list of f_n functions. The will have fractions
            #multiply both sides of the equation by denominators
            ##for i in range(len(f_n)):
            ##    n,d=sympy.fraction(sympy.together(f_n[i]))
            ##    f_n[i]=sympy.simplify(sympy.Mul(d,f_n[i]))
            #f_n=[sympy.simplify(sympy.Mul(sympy.sympify("*".join(one)),f_n_rational[i])) for i in range(len(f_n_rational))]
            equations = equations + f_n[:3]
    for face_i,face in enumerate(faces): # special initial value assignments for bigons
        if len(face) == 2:
            for edge in face:
                equations.append(edge_vars[edge_indices[edge]]) # enforce edge_label == 0 for sides of 2 sided regions
                edge_labels[edge_indices[edge]] = 0 # Reset edge label to 0
            equations.append(crossing_vars[face[0][0][0].label]-crossing_vars[face[1][0][0].label]) #ALEX: added crossing relation for bigons.
    print1("We run Newton-Raphson method to solve the following system of equations:")
    for eq in equations:
        print1("{} = 0".format(eq))
    f_temp = sympy.lambdify(crossing_vars+edge_vars,equations,'numpy')
    Jac = sympy.Matrix([[sympy.diff(eq,var) for var in crossing_vars + edge_vars] for eq in equations])
    Df_temp = sympy.lambdify(crossing_vars + edge_vars, Jac,'numpy')
    f = lambda x : f_temp(*x)
    Df = lambda x: Df_temp(*x)
    newton_iter = 30   #number of ierations in Newton-Raphson method
    for _ in range(newton_iter):
        f_x = f(labels)
        Df_x = Df(labels)
        q,r = np.linalg.qr(Df_x)
        incx = np.linalg.solve(r,-np.matmul(q.conj().T,f_x))
        labels = np.around(labels+incx,6)
    
      
    







def analyze_alt_sol(L):
    # Collect edges.  Each edge is a (CrossingStrand, CrossingEntryPoint) pair, where the latter is the 'head' of the oriented edge.
    # Note that comparing a (CrossingStrand, CrossingEntryPoint) pair to a (CrossingStrand, CrossingStrand) pair (as are given in the faces) still works properly since CrossingEntryPoint inherits from CrossingStrand.
    oriented_edges = [(x.opposite(),x) for x in L.crossing_entries()]
    # every edge appears once as-is and once reversed in the faces
    faces = [[(x.opposite(),x) for x in face] for face in L.faces()] #hopefully right orientation
    # We create a local name for the set of crossings for consistency
    crossings = L.crossings

    # First color all the regions.  The colors are referenced by the list 'coloring' which stores the colors as strings "black" or "white"
    coloring = ['']*len(faces)
    coloring[0] = 'black'
    to_add = [faces.index(adj) for adj in get_adj_faces_D(0,faces)]
    for to_add_face in to_add:
        coloring[to_add_face] = 'white'
    while(to_add):
        face_i = to_add.pop()
        for nbr in get_adj_faces_D(face_i,faces):
            nbr_i = faces.index(nbr)
            if not coloring[nbr_i]:
                coloring[nbr_i] = opposite_color_D(coloring[face_i])
                to_add.append(nbr_i)

    # Next we iterate through all edges by iterating through black faces and assigning indices and names to the corresponding edge
    # At the same time, we assign indices and names to the reversed edge (bordering a white region)

    edge_indices = {} # used to index the edge in arrays or lists (such as 'edge_names' below)
    edge_names = [] # human readable names

    count = 0
    for face_i,face in enumerate(faces):
        if coloring[face_i] == 'white':
            continue
        for edge in face:
            edge_indices[edge] = count*2
            edge_indices[reverse_edge_D(edge)] = count*2+1
            edge_names.append('u' + str(count+1))
            edge_names.append('v' + str(count+1))
            count += 1

    #ALEX NOTES: constructor then zeroing out a sq matrix of complex numbers."
    labels = np.zeros(len(crossings) + 2*len(oriented_edges), dtype=np.complex128) #initial value estimates for all variables
    crossing_labels = labels[:len(crossings)] # reference crossing labels by slice
    edge_labels = labels[len(crossings):] # reference edge labels by slice
    edge_vars = sympy.symbols(edge_names)
    crossing_vars = sympy.symbols(['w' + str(k) for k in range(len(crossing_labels))])

    #labels = np.zeros(len(crossings+2*len(oriented_edges)), dtype=np.complex_)
    equations = []
    

    for crossing in L.crossings: # initial crossing label values
        if crossing.sign == 1: #-0.5
            crossing_labels[crossing.label] = .5+.5j  # .5j
        else:
            crossing_labels[crossing.label] = -.5-.5j  # -.5j
    # Now lets go through each face and create the equations
    count=0
    for face_i,face in enumerate(faces):
        if len(face) < 2:
            return
        if coloring[face_i] == 'black': # Add edge label equations (for the two sides of each edge) only if the color is black
            for edge in face:
               edge_labels[edge_indices[edge]] = -.5-.5j # initial edge label value for black side
               edge_labels[edge_indices[reverse_edge_D(edge)]] = .5-.5j #initial edge label value for white side
               k = compute_k_D(edge) #Set this to +-1 appropriately as per paper
               equations.append(edge_vars[edge_indices[edge]]-edge_vars[edge_indices[reverse_edge_D(edge)]]-k)

        if len(face) > 2: # generate equations for faces with more than 2 sides
            shape_params = []
            kl = compute_kappa_D(face[0],face[-1],oriented_edges)
            el = edge_vars[edge_indices[face[-1]]]
            e1 = edge_vars[edge_indices[face[0]]]
            w1 = crossing_vars[face[0][0][0].label]
            shape_par = kl*w1/el/e1
            shape_params.append(shape_par)
            for i in range(len(face) - 1):
                e_i = edge_vars[edge_indices[face[i]]]
                e_ii = edge_vars[edge_indices[face[i+1]]]
                w_i = crossing_vars[face[i+1][0][0].label]
                k_i = compute_kappa_D(face[i],face[i+1],oriented_edges)
                shape_par = k_i*w_i/e_i/e_ii
                shape_params.append(shape_par)
            edge_prod = 1
            for e in face:
                edge_prod *= edge_vars[edge_indices[e]]
            ##one=[]  #code error
            ##for i in range(len(face)):
            ##    one.append(str(edge_vars[edge_indices[face[i]]]))
            f_n= get_f_n_D(tuple(shape_params),len(shape_params)) # get list of f_n functions. The will have fractions
            #multiply both sides of the equation by denominators
            ##for i in range(len(f_n)):
            ##    n,d=sympy.fraction(sympy.together(f_n[i]))
            ##    f_n[i]=sympy.simplify(sympy.Mul(d,f_n[i]))
            #f_n=[sympy.simplify(sympy.Mul(sympy.sympify("*".join(one)),f_n_rational[i])) for i in range(len(f_n_rational))]
            equations = equations + f_n[:3]
    for face_i,face in enumerate(faces): # special initial value assignments for bigons
        if len(face) == 2:
            for edge in face:
                equations.append(edge_vars[edge_indices[edge]]) # enforce edge_label == 0 for sides of 2 sided regions
                edge_labels[edge_indices[edge]] = 0 # Reset edge label to 0
            equations.append(crossing_vars[face[0][0][0].label]-crossing_vars[face[1][0][0].label]) #ALEX: added crossing relation for bigons.
    
    

    #for eq in equations:
    #    print1("{} = 0".format(eq))

    
    f_temp = sympy.lambdify(crossing_vars+edge_vars,equations,'numpy')
    Jac = sympy.Matrix([[sympy.diff(eq,var) for var in crossing_vars + edge_vars] for eq in equations])
    Df_temp = sympy.lambdify(crossing_vars + edge_vars, Jac,'numpy')
    f = lambda x : f_temp(*x)
    Df = lambda x: Df_temp(*x)
    newton_iter = 30   #number of ierations in Newton-Raphson method
    for _ in range(newton_iter):
        f_x = f(labels)
        Df_x = Df(labels)
        q,r = np.linalg.qr(Df_x)
        incx = np.linalg.solve(r,-np.matmul(q.conj().T,f_x))
        labels = np.around(labels+incx,6)
    #output time
    print1("\nThe face with vertices {} is black.  Other faces can be colored in a checkerboard fashion.\n'u' labels correspond to the black sides of the edges, and 'v' labels to the white sides.".format([edge[1][0] for edge in faces[0]]))
    for i in range(len(crossings)):
        print1("w{} = {}".format(i,labels[i]))
    count = 0
    print1('\nCalculated values:')
    for face_i,face in enumerate(faces):
        print1("\n\nFace Number: {}".format(face_i))
        for edge in face:
            print1("crossing {} to crossing {}".format(edge[0][0],edge[1][0]))
        if coloring[face_i] == 'white':
            print1("Face Color: White")
            continue
        print1("Face Color: Black\n")
        for edge in face:
            print1('Edge from crossing {} to crossing {}'.format(edge[0][0],edge[1][0]))
            print1('v{} = {}'.format(str(count+1),labels[len(crossings)+count*2+1]))
            count += 1





def main():
    if __name__ == "__main__":
        print1("Input DT (Dowker-Thislethwaite) code for your link diagram. See for example, https://en.wikipedia.org/wiki/Dowker%E2%80%93Thistlethwaite_notation. If the link is alternating, PD (Planar Diagram) code can be used as well.")
        print1("Example input for a figure-eight knot diagram: 'DT: [(4,6,8,2)]'. Example input for Borromean Rings link diagram: 'DT: [(6,8),(10,12),(4,2)]'.")
        print1("")
        print1("Input a link code:")

    # temp_inp = input()
    temp_inp = input()
    if temp_inp.startswith("DT: "):
        L = links.Link(temp_inp)
        print1(f"Analyzing link with DT code {L.DT_code()}.")
    else:
        L = links.Link(ast.literal_eval(temp_inp))
        print1(f"Analyzing link with PD code {L.PD_code()}")
    
    
    while True:
        print1("\nMenu:")
        print1("1. Find equations for the complete hyperbolic structure of a knot/link.")
        print1("2. Find values of edge and crossing labels describing the complete hyperbolic structure of a knot/link. (Based on the work by Thislethwaite and Tsvietkova: https://arxiv.org/abs/1108.0510.)")
        print1("3. Find equations for the canonical component of the PSL(2, C)-character variety of a knot.(Based on the work by Petersen and Tsvietkova.) Note: for 1-3, the given knot/link diagram should be taut.")
        print1("4. Exit the program.")

        choice = input("Enter your choice (1, 2, 3, or 4): ")
        
        if choice == '1':
            analyze_alt_eq(L)
        elif choice == '2':
            analyze_alt_sol(L)
        elif choice == '3':
            analyze_PSL(L)
        elif choice == '4':
            print1("Exiting the program.")
            break
        else:
            print1("Invalid input. Please enter 1, 2, 3, or 4.")

main()
