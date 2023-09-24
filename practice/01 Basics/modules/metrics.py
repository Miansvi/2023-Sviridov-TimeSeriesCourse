import numpy as np


def ED_distance(ts1, ts2):
    """
    Calculate the Euclidean distance.

    Parameters
    ----------
    ts1 : numpy.ndarray
        The first time series.

    ts2 : numpy.ndarray
        The second time series.

    Returns
    -------
    ed_dist : float
        Euclidean distance between ts1 and ts2.
    """
    
    ed_dist = np.sqrt(np.sum(np.power(ts1 - ts2, 2)))

    # INSERT YOUR CODE
    
    return ed_dist


def norm_ED_distance(ts1, ts2):
    """
    Calculate the normalized Euclidean distance.

    Parameters
    ----------
    ts1 : numpy.ndarray
        The first time series.

    ts2 : numpy.ndarray
        The second time series.

    Returns
    -------
    norm_ed_dist : float
        The normalized Euclidean distance between ts1 and ts2.
    """

    norm_ed_dist = 0

    # INSERT YOUR CODE 

    return norm_ed_dist


def DTW_distance(ts1, ts2, r=None):
    """
    Calculate DTW distance.

    Parameters
    ----------
    ts1 : numpy.ndarray
        The first time series.

    ts2 : numpy.ndarray
        The second time series.

    r : float
        Warping window size.
    
    Returns
    -------
    dtw_dist : float
        DTW distance between ts1 and ts2.
    """

    dtw_dist = 0
    N, M = sum(ts1.shape), sum(ts2.shape)
    dist_mat=np.zeros((N,M))
    for i in range(N):
      for j in range(M):
        dist_mat[i,j] = (ts1[i]- ts2[j])**2

    N,M=dist_mat.shape
    D_mat = np.zeros((N+1,M+1))
    for i in range(1,N+1):
        D_mat[i,0]=np.inf
    for i in range(1,M+1):
        D_mat[0,i]=np.inf

    for i in range(1,N+1):
      for j in range(1,M+1):
        D_mat[i][j] = dist_mat[i-1][j-1]+min(D_mat[i-1][j],D_mat[i,j-1],D_mat[i-1][j-1])
    return  D_mat[N][M]