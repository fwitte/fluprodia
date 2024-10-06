import numpy as np


def _create_pattern_from_base(base_pattern, minimum, maximum):
    """Create a pattern spanning over the magnitudes between minimum and
    maximum given a specific base pattern for every magnitude

    Parameters
    ----------
    base_pattern : numpy.ndarray
        base pattern per magnitude
    minimum : float
        minimum of the range
    maximum : float
        maximum of the range

    Returns
    -------
    numpy.ndarray
        pattern repeating the base pattern for every magnitude
    """
    magnitude_min = np.floor(np.log10(minimum))
    magnitude_max = np.ceil(np.log10(maximum))
    magnitudes = (10 ** np.arange(magnitude_min, magnitude_max + 1))[..., None]
    return (np.array(base_pattern) * magnitudes).flatten()


def _generate_pattern_candidate(base_pattern, minimum, maximum, include_endpoints):
    """Create possible patterns for a base pattern

    Parameters
    ----------
    base_pattern : numpy.ndarray
        pattern repeating the base pattern for every magnitude
    minimum : float
        minimum of the range
    maximum : float
        maximum of the range
    include_endpoints : bool
        Make the range start below or at the minimum and extend to the maximum,
        or above

    Returns
    -------
    _type_
        _description_
    """

    pattern = _create_pattern_from_base(base_pattern, minimum, maximum)

    last_below_range = np.argwhere(pattern <= minimum + 1e-6)
    if len(last_below_range) == 0:
        last_below_range = None
    else:
        last_below_range = last_below_range[-1][0]
    first_above_range = np.argwhere(pattern >= maximum - 1e-6)

    if len(first_above_range) == 0:
        first_above_range = None
    else:
        first_above_range = first_above_range[0][0] + 1

    if not include_endpoints:
        if last_below_range is not None:
            last_below_range += 1
        else:
            last_below_range = 0

        if first_above_range is not None:
            first_above_range -= 1
        else:
            first_above_range = len(pattern)
    return pattern[last_below_range:first_above_range]


def _log_range(minimum, maximum, include_endpoints=True, target_num_points=8):
    """Create a roughly log spaced range with nice numbers.

    Parameters
    ----------
    minimum : float
        Minimum value, not necessarily included in range
    maximum : float
        Maximum value, not necessarily included in range
    include_endpoints : bool, optional
        Make the range start below or at the minimum and extend to the maximum,
        or above, by default True
    target_num_points : int, optional
        Target number of points in the range, by default 8

    Returns
    -------
    numpy.ndarray
        Roughly log spaced range
    """
    patterns = [
        [1],
        [1, 3],
        [1, 2, 5],
        [1, 2, 3, 5],
        [1, 2, 4, 6, 8],
        [1, 1.5, 2, 3, 5, 8],
        [1, 2, 3, 4, 5, 6, 7, 8, 9]
    ]
    min_distance = 1000

    for i, pattern in enumerate(patterns[::-1]):

        pattern = _generate_pattern_candidate(pattern, minimum, maximum, include_endpoints)
        distance = len(pattern) - target_num_points
        if abs(distance) < min_distance:
            if distance == 0:
                break
            min_distance = distance
        elif distance < 0:
            i = i - 1
            break

    return _generate_pattern_candidate(patterns[::-1][i], minimum, maximum, include_endpoints)


def _linear_range_excluding_endpoints(minimum, maximum, stepsizes):
    """Create a linear range which does not contain minimum and maximum.

    Parameters
    ----------
    minimum : float
        Minimum value, not necessarily included in range
    maximum : float
        Maximum value, not necessarily included in range
    stepsizes : numpy.ndarray
        Possible base step sizes.

    Returns
    -------
    tuple
        tuple of possible minima and maxima given the base step sizes.
    """
    minima = minimum + stepsizes - minimum % stepsizes
    maxima = maximum - maximum % stepsizes
    return minima, maxima


def _linear_range_including_endpoints(minimum, maximum, stepsizes):
    """Create a linear range which does contain minimum and maximum.

    Parameters
    ----------
    minimum : float
        Minimum value, not necessarily included in range
    maximum : float
        Maximum value, not necessarily included in range
    stepsizes : numpy.ndarray
        Possible base step sizes.

    Returns
    -------
    tuple
        tuple of possible minima and maxima given the base step sizes.
    """
    minima = minimum - minimum % stepsizes
    maxima = maximum + stepsizes - maximum % stepsizes
    return minima, maxima


def _linear_range(minimum, maximum, include_endpoints=True, target_num_steps=8):
    """Create a linearly spaced range with nice numbers.

    Parameters
    ----------
    minimum : float
        Minimum value, not necessarily included in range
    maximum : float
        Maximum value, not necessarily included in range
    include_endpoints : bool, optional
        Make the range start below or at the minimum and extend to the maximum,
        or above, by default True
    target_num_points : int, optional
        Target number of points in the range, by default 8

    Returns
    -------
    numpy.ndarray
        Linearly spaced range
    """
    base_steps = np.array([1, 2, 2.5, 3, 4, 5, 7.5])

    step_magnitude = np.floor(np.log10(maximum - minimum))
    stepsizes = 10 ** (step_magnitude - 1) * base_steps
    stepsizes2 = 10 ** step_magnitude * base_steps
    stepsizes = np.append(stepsizes, stepsizes2)

    if include_endpoints:
        minima, maxima = _linear_range_including_endpoints(minimum, maximum, stepsizes)
    else:
        minima, maxima = _linear_range_excluding_endpoints(minimum, maximum, stepsizes)

    mask = maximum % stepsizes != 0
    maxima[mask] += stepsizes[mask]

    min_distance = 100
    for i in range(len(stepsizes)):
        pattern = np.arange(minima[i], maxima[i], stepsizes[i])
        distance = len(pattern) - target_num_steps
        if abs(distance) < min_distance:
            if distance == 0:
                break
            min_distance = distance
        elif distance < 0:
            i = i - 1
            break

    return np.arange(minima[i], maxima[i], stepsizes[i])


def _isolines_log(val_min, val_max):
    """Generate default logarithmic isolines.

    Parameters
    ----------
    val_min : float
        Minimum value for isoline range.

    val_max : float
        Maximum value for isoline range.

    Returns
    -------
    arr : ndarray
        numpy array with logarithmically spaced values starting from the
        minimum value going to the maximum value in steps of :code:`1ek`,
        :code:`2ek` and :code:`5ek`.
    """
    arr = [val_min]
    digits = int(np.floor(np.log10(val_min)))
    while arr[-1] < val_max:
        arr += [1 * 10 ** digits]
        arr += [2 * 10 ** digits]
        arr += [5 * 10 ** digits]
        digits += 1

    arr = np.unique(np.asarray(arr + [val_max]))
    return arr[(arr >= val_min) & (arr <= val_max)]


def _beautiful_unit_string(unit):
    r"""Convert unit fractions to latex.

    Parameters
    ----------
    unit : str
        Value of unit for input, e.g. :code:`m^3/kg`.

    Returns
    -------
    unit : str
        Value of unit for output, e.g. :code:`$\frac{m^3}{kg}$`.
    """
    if '/' in unit:
        numerator = unit.split('/')[0]
        denominator = unit.split('/')[1]
        unit = '$\\frac{' + numerator + '}{' + denominator + '}$'

    return unit

