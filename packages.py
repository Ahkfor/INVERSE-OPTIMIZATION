import gurobipy as gp
from gurobipy import GRB
from helper_functions import *



def RP(OS_parameters, X_hat, S):
    '''
    Function to solve the RP problem
    :param OS_parameters: include the constraint matrix and the objective matrix
    :param X_hat: an input data set
    :param S: a list of binary values representing an objective set
    :return: the optimal value of the RP problem
    '''

    # Supress output
    env = gp.Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.start()

    A_mat = OS_parameters[0]
    obj = OS_parameters[1]

    # Get the number of objectives
    (num_obj, num_var) = obj.shape

    # Get the number of constraints and variables
    (num_constraint, num_var) = A_mat.shape

    # Get number of x solutions in the list
    num_sol = len(X_hat)

    RP = gp.Model('mip1', env = env)

    # Add the variables
    # esp
    esp = RP.addMVar(num_sol, name='esp')
    # x
    x = RP.addMVar((num_sol, num_var), name = 'x', lb = 0.1, ub = 2)
    # Objective
    RP.setObjective(sum(esp[i] for i in range(num_sol)), GRB.MINIMIZE)

    # The duality gap constraints
    for k in range(num_sol):
        for i in range(num_obj):
            if S[i] != 0:
                RP.addConstr(esp[k]*(sum(obj[i, j]*X_hat[k][j] for j in range(num_var)))>= sum(obj[i, j]*x[k, j] for j in range(num_var)))


    # The Ax <= b constraints
    for k in range(num_sol):
        for i in range(num_constraint):
            RP.addConstr((sum(A_mat[i, j]*x[k, j] for j in range(num_var))) >= 2)

    RP.optimize()

    # Get the return values
    solution = RP.getAttr('x')
    esp = solution[0: num_sol]
    return esp


def algorithm5(K, X_hat, k_prime, theta, tal):
    '''
    Function to implement algorithm 5
    :param OS_parameters: parameters of the forward problem
    :param X_hat: input data set
    :param k_prime: prespecified objective
    :param theta: number of objectives to be chosen
    :return: an objective set S
    '''

    # Get number of samples
    N = len(X_hat)

    # Initiate S
    S = [0]*K
    S[k_prime] = 1

    eta_star = float('inf')
    eta = float('inf')
    i = 1
    while i <= theta:
        for k in range(K):
            if S[k] == 1:
                continue
            # Store all eps values
            eps_li = []
            # S = S union k
            S[k] = 1
            # eps_li is a list of eps values
            eps_li = RP(X_hat, S)
            eps_li = inverse_list(eps_li)
            if max(eps_li) <= eta:
                eta = max(eps_li)
                k_star = k
            S[k] = 0
        S[k_star] = 1
        eta_star = eta
        i += 1
    return S, eta_star


def algorithm6(K, theta, tal, X_hat_li):
    '''
    This function will generate initial I^q sets
    :param K: number of objectives
    :param theta: number of objectives to be chosen
    :param tal: relaxation level
    :param X_hat_li: data set containing Q samples
    :return: I_mat, a binary matrix; Z, a binary matrix where each column is an objective set
    '''
    Q = len(X_hat_li)
    I_li = []
    for j in range(Q):
        I_li.append([])
    Z = []
    i = 1
    for q in range(Q):
        for k in range(K):
            S, eta_star = algorithm5(K, X_hat_li[q], k, theta, tal)
            if eta_star <= tal:
                if S not in Z:
                    Z.append(S)
                    I_li[q].append(i)
                    i += 1
                else:
                    I_li[q].append(Z.index(S))

    # Convert I_li to a binary matrix
    I = len(Z)
    I_mat = convert_binary(I_li, I)
    Z = np.array(Z).transpose()

    return I_mat, Z


def SOS_LIN(I_mat, Z):
    '''
    This function solves SOS-LIN problem
    :param I_li: a list of lists of integers
    :return:
    '''

    # Supress output
    env = gp.Env(empty=True)
    env.setParam("OutputFlag", 0)
    env.start()

    _, K = Z.shape
    I, Q = I_mat.shape
    print(I_mat)


    # Build model
    SOS = gp.Model()
    y = SOS.addMVar((Q, I), lb=0, ub=1)

    # # Set Objective
    # chosen_li = []
    # for q in range(Q):
    #     chosen = [0]*K
    #     for j in range(len(I_li[q])):
    #         chosen = np.add(chosen, y[j, q]*Z[I_li[q][j], :])
    #     chosen_li.append(chosen)
    # chosen_li = np.array(chosen_li)

    print(stability(y@Z))
    SOS.setObjective(stability(y@Z))

    # Add constraints
    for q in range(Q):
        SOS.addConstr(sum(y[j, q] for j in range(len(I))) == 1)
        for j in range(len(I)):
            if j not in I_li[q]:
                SOS.addConstr(y[j, q] == 0)

    SOS.optimize()
    # Retrieve y solution
    y_sol = SOS.getAttr('x')
    y_sol = y_sol.reshape((len(I), Q))

    return y_sol


# def algorithm7(OS_parameters, K, theta, tal, X_hat_li):
#     I_li, Z = algorithm6(K, theta, tal, X_hat_li)
#     Q = len(X_hat_li)
#     for q in range(Q):
#         for q1 in range(Q):
#             if q != q1:
#                 for i in range(len(I_li[q1])):
#                     eps_li = RP(OS_parameters, X_hat_li[q], Z[i, :])
#                     eps_li = inverse_list(eps_li)
#                     if max(eps_li) <= tal:
#                         I_li[q].append(i)
