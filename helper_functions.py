import numpy as np




def inverse_list(eps):
    '''
    This is a helper function that convert each element of the list to its inverse
    :param eps: a list of numbers
    :return:
    '''
    for i in range(len(eps)):
        eps[i] = 1/eps[i]
    return eps


def stability(Z):
    '''
    This function implements the phi function
    :param Z: a matrix (can have fractional values)
    :return: stability of Z
    '''
    frequency = np.mean(Z, axis=0)
    theta = sum(Z[0, :])
    s_2 = []
    Q, K = Z.shape
    for i in range(K):
        s_2.append(Q*frequency[i]*(1-frequency[i])/(Q-1))
    phi = 1 - sum(s_2)/(theta*(1-theta/K))
    return phi


def convert_binary(I_li, I):
    '''
    This is a helper function that convert a list of list of integers to binary matrix
    :param I_li: list of lists of integers
    :param I: an integer
    :return: an I*Q binary matrix, each column represent choices of objective sets for a sample
    '''


    Q = len(I_li)
    binary_matrix = np.zeros((I, Q))
    for q in range(Q):
        for j in range(len(I_li[q])):
            binary_matrix[I_li[q][j], q] = 1

    return binary_matrix
