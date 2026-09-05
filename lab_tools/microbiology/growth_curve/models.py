"""Reparameterized logistic and Gompertz growth models (Zwietering)."""
import numpy as np


def logistic(t, A, mu, lag):
    """
    Zwietering logistic growth model.

    y(t) = A / (1 + exp[(4·μ/A)·(lag - t) + 2])

    Parameters
    ----------
    t : array_like
        Time points.
    A : float
        Asymptotic maximum (carrying capacity).
    mu : float
        Maximum specific growth rate.
    lag : float
        Lag time.

    Returns
    -------
    np.ndarray
        Predicted OD values.
    """
    inner = (4.0 * mu / A) * (lag - t) + 2.0
    return A / (1.0 + np.exp(inner))


def gompertz(t, A, mu, lag):
    """
    Zwietering Gompertz growth model.

    y(t) = A · exp{-exp[(μ·e/A)·(lag - t) + 1]}

    Parameters
    ----------
    t : array_like
        Time points.
    A : float
        Asymptotic maximum.
    mu : float
        Maximum specific growth rate.
    lag : float
        Lag time.

    Returns
    -------
    np.ndarray
        Predicted OD values.
    """
    inner = (mu * np.exp(1.0) / A) * (lag - t) + 1.0
    return A * np.exp(-np.exp(inner))
