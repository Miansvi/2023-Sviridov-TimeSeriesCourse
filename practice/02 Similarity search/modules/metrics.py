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

def calc_sd(ts, mu):
  sd = 0
  for i in ts:
    sd += i**2 - mu**2
  
  return np.sqrt((1.0/len(ts)) * sd)



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
    n = len(ts1)

    mu_ts1 = 1.0/n * np.sum(ts1)
    mu_ts2 = 1.0/n * np.sum(ts2)

    sd_ts1 = calc_sd(ts1, mu_ts1)
    sd_ts2 = calc_sd(ts2, mu_ts2)

    norm_ed_dist = np.sqrt(2 * n * (1 - (np.dot(ts1, ts2) - n * mu_ts1 * mu_ts2) / (n * sd_ts1 * sd_ts2)))

    # INSERT YOUR CODE 

    return norm_ed_dist


def DTW_distance(ts1, ts2, r=0):
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
    N, M = len(ts1), len(ts2)
    
    dist_mat=np.zeros((N,M))
    for i in range(N):
      for j in range(M):
        dist_mat[i,j] = (ts1[i]- ts2[j])**2

    D_mat = np.zeros((N+1,M+1))
    D_mat[:, :] = np.inf
    D_mat[0, 0] = 0

    r_int = int(np.floor(M*r))
    for i in range(1,N+1):
      left = max(1, i - r_int) if r > 0 else 1
      right = (min(M, i + r_int) + 1) if r > 0 else M + 1
      for j in range(left,right):
        D_mat[i][j] = dist_mat[i-1][j-1]+min(D_mat[i-1][j],D_mat[i,j-1],D_mat[i-1][j-1])
    return  D_mat[N][M]