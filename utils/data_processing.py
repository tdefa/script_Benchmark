



import numpy as np


def sctransform_from_parameters(fit_params, array_of_vect):

    """
    inspire from the init files of ssam, compute sc transform from th...
    :param fit_params: directly from the r script (gene, 3)
    :param array_of_vect: array shape (nb_of_vect, nb_of_genes)
    :return: the normalize version of array_of_vect
    """

    fit_params = np.array(fit_params).T
    nvec = array_of_vect.shape[0]
    regressor_data = np.ones([nvec, 2])
    regressor_data[:, 1] = np.log10(np.sum(array_of_vect, axis=1))

    mu = np.exp(np.dot(regressor_data, fit_params[1:, :]))
    with np.errstate(divide='ignore', invalid='ignore'):
        res = (array_of_vect - mu) / np.sqrt(mu + mu ** 2 / fit_params[0, :])

    return res


### input preparation





