import numpy as np
import IsoSpecPy
import IsoSpecPy.Distributions
import math
import cmath
# import masserstein

from scipy.optimize import minimize
from pyteomics import mzml

formula = {
    'averagine': np.array([4.92453424, 7.77240127, 1.35557053, 1.46005576, 0.03567328]),  # mass: 110.4728 Da
    'averagine_normed': np.array([0.523082530, 0.825582100, 0.143988290, 0.155086679, 0.003789205]),
    'variance_changer': np.array([-0.34306864, -0.03723128, -0.49608047, -0.51840443, 0.60504442]),  # vector that change variance without changing average mass
    'variance': np.array([0.0107433995, 0.0001171534, 0.0036072342, 0.0086044585, 0.1700188292]),  # mass: 0.09292879 Da
    'variance_normed': np.array([0.0629693102, 0.0006866606, 0.0211427537, 0.0504325301, 0.9965158997]),
    'mass_avg': np.array([12.010824, 1.007941, 14.006705, 15.999409, 32.064888]),
    'mass_avg_normed': np.array([0.29789673, 0.02499932, 0.34739928, 0.39682302, 0.79528473]),
    'mass_monoisotopic': np.array([12.0, 1.00782503227, 14.0030740042, 15.9949146202, 31.9720711741])
}

######################################
### THEORETICAL SPECTRA PREDICTION ###
######################################

def final_rounding(protein, zeta, delta):
    '''
    Gives linear model prediction rounded to grid defined by zeta and delta.
    '''
    lm_pred = initial_prediction(protein)

    l = np.floor((lm_pred - delta) / zeta) * zeta + delta
    if lm_pred - l < l + zeta - lm_pred:
        alm_pred = l
    else:
        alm_pred = l + zeta
    return alm_pred - 1.198223e-07 * lm_pred


def initial_prediction(protein,
                      intercept=-0.1455702755046254237570,
                      mass_coef=0.9997869891271233822039,
                      var_coef=-0.5981736335758824907316):
    '''
    Monoisotopic peak mass prediction (initial) by linear model.
    '''
    mass = protein.empiric_average_mass()
    var = protein.empiric_variance()
    return intercept + mass_coef * mass + var_coef * var


def mzMLread(file):
    '''
    Read mzML file as dictionary of lists of pairs mass-intensity.
    '''
    reader = mzml.read(file)
    output = dict()
    for i in reader:
        output[i['index']] = [
            np.array(j) for j in [*zip(i['m/z array'], i['intensity array'])]
        ]
    return output


def optimal_delta2(iso, zeta):
    '''
    Compute optimal delta parameter (grid shift) for given theoretical spectrum and zeta parameter.
    '''
    masses = iso.np_masses()

    circle = []
    for i in range(len(masses)):
        circle.append((math.e**(2 * math.pi * 1j * masses[i] / zeta)))
    mean = sum(l[0] * l[1] for l in zip(circle, iso.np_probs()))
    delta = zeta / (2 * math.pi) * cmath.log(mean).imag
    return delta + zeta * (delta < 0)


def optimal_zeta(protein, prob_to_cover=.99):
    '''
    Compute zeta that minimizes variance.
    '''
    iso = IsoSpecPy.IsoTotalProb(prob_to_cover, protString(protein))
    return minimize(variance_for_zeta, 1.000514, args=(iso)).x[0]


def protString(formula):
    '''
    Convert np.array (CHNOS) chemical formula to string 'C...H...N... ...'.
    '''
    formula = list(map(int, formula))
    result = ''
    for i in range(5):
        result += 'CHNOS'[i] + str(formula[i])
    return result


def variance_for_zeta(zeta, iso):
    '''
    Compute variance of spectrum on circle with circumference equal to zeta in complex space.
    '''
    masses = iso.np_masses()
    probs = iso.np_probs()
    c1 = zeta / (2 * math.pi * 1j)
    c2 = 1 / c1

    circle = []
    for i in range(len(masses)):
        circle.append((math.e**(c2 * masses[i])))
    mean_peak = sum(l[0] * l[1] for l in zip(circle, probs)) # average peak in complex space
    mean_teta = math.e**(-cmath.log(mean_peak).imag * 1j) # how much circle is rotated from zero
    circle = [l * mean_teta for l in circle] # rotate to zero

    projection = [(c1 * cmath.log(l)).real for l in circle]
    mean = sum(l[0] * l[1] for l in zip(projection, probs))
    variance = sum(
        [probs[i] * (projection[i] - mean)**2 for i in range(len(masses))])
    return variance


def zeta_linear(mass_avg):
    '''
    Compute zeta parameter by use of linear model for given protein average mass.
    '''
    return 1.0023554156 + 6.9584683113e-10 * mass_avg


#######################################
### SIMULATED SPECTRUM CONSTRUCTION ###
#######################################

def get_averagine(input_mass, grid_example=None, variance_shift=0):
    '''
    Returns spectrum of averagine with changed variance. Parameter variance_shift determine, how much shift chemical formula into direction of vector that changes vawiance without changeing average mass. If grid_example is given, then the spectrum is shifted so, that the highest peak to the left of average mass but no further than zeta [Da] (if possible) is in grid_example place.
    '''
    in_mass = input_mass/110.4728
    chemical_formula = np.round(in_mass * formula['averagine'] + variance_shift * formula['variance_changer'])
    if chemical_formula[4] < 0: # in case of lack of more sulfur, we subtract equivalent in carbons
        chemical_formula[0] += np.round(chemical_formula[4] * 2.6696)
    chemical_formula[chemical_formula < 0] = 0

    mass = np.dot(chemical_formula, formula['mass_avg'])
    chemical_formula[1] += np.round(input_mass - mass)
    if chemical_formula[1] < 0:
        return None

    spectrum = IsoSpecPy.IsoTotalProb(0.99, protString(chemical_formula))
    spectrum.normalize()
    if grid_example is not None:
        spectrum.add_mass(grid_example - get_central_mass(spectrum))

    return spectrum


def get_central_mass(experimental, zeta=1.002381):
    '''
    Returns highest peak that is to the left from average mass and no further than zeta from it.
    If there is no peak in the interval, returns highest peak over whole spectrum.
    '''
    mass_avg = experimental.empiric_average_mass()
    masses = [*experimental.masses]
    intensities = [*experimental.probs]

    s = zip(masses, intensities)
    t = [j for j in s if j[0] > mass_avg-zeta and j[0] < mass_avg]
    if len(t) > 0:
        return max(t, key = lambda x: x[1])[0]
    else:
        return max(s, key = lambda x: x[1])[0]


def get_most_abu(spectrum):
    '''
    Returns most abundant peak of a given spectrum.
    '''
    masses = [*spectrum.masses]
    intensities = [*spectrum.probs]
    return masses[np.argmax(intensities)]


def monoisotopic_mass_prediction(experimental_masses, experimental_intensities, charge=None, return_simulated_spectrum=False):
    '''
    Function than do whole envemind monoisotopic mass prediction. Takes masses and intensities of experimental spectrum wor which prediction have to be run (experimental_masses and experimental_intensities parameters). If spectrum have m/z instead of daltons, then charge should be given. Function takes only fragment of spectrum (short range in Da or m/z) with single substance (already deconvoluted, if needed). If parameter return_simulated_spectrum is turned to True, function returns spectrum that fitted to experimental best.
    '''
    if charge is not None:
        experimental_masses = [(m - formula['mass_avg'][1]) * charge for m in experimental_masses]

    experimental_spectrum = IsoSpecPy.IsoDistribution(masses = experimental_masses, probs = experimental_intensities)
    experimental_spectrum.normalize()

    experimental_MostAbu = get_most_abu(experimental_spectrum)
    zeta = zeta_linear(experimental_MostAbu)

    simulated = optimize_position(experimental_spectrum, experimental_MostAbu, zeta)

    if return_simulated_spectrum == True:
        return simulated

    delta = optimal_delta2(simulated, zeta)
    return final_rounding(simulated, zeta, delta)


def optimize_position(experimental, experimental_MostAbu, zeta):
    '''
    Function returns spectrum that will be used to prediction of monoisotopic mass. Function optimizes for which shift averagine with optimized variance is most similar to experimental spectrum (in measure defined in spectraLikeness function). To ommit repeating calculations, most abudant peak mass have to be given (experimental_MostAbu parameter) and zeta.
    '''
    spectrum_best, score_best = optimize_variance(experimental, experimental_MostAbu, experimental_MostAbu)
    shift_best = 0
    def shift_score(shift):
        return optimize_variance(experimental, experimental_MostAbu, experimental_MostAbu+shift*zeta)

    spectrum_left, left = shift_score(-1)
    spectrum_right, right = shift_score(1)
    if left < right:
        shift = -2
        score = left
        spectrum = spectrum_left
    else:
        shift = 2
        score = right
        spectrum = spectrum_right

    if score_best > score:
        score_best = score
        spectrum_best = spectrum
        while True:
            spectrum, score = shift_score(shift)
            if score_best > score:
                score_best = score
                spectrum_best = spectrum
                shift += np.sign(shift)
            else:
                break

    return spectrum_best


def optimize_variance(experimental, experimental_avg, center_place, variance_parameter=1.6527):
    '''
    Function optimizes how much of vector that change variance have to be added to averagine, to obtain spectrum most similar to experimental one. Returns best spectrum and its score (in measure defined in spectraLikeness function). Averagine is calculated with given average mass (experimental_avg parameter), and centered in center_place (like described in get_averagine function). Parameter variance_parameter determines size of step in optimization proces, and should not be changed without important reason.
    '''
    simulated = get_averagine(experimental_avg, center_place)
    score_best = spectraLikeness(experimental, simulated)
    shift_best = 0

    minus = spectraLikeness(experimental, get_averagine(experimental_avg, center_place, -variance_parameter))
    plus = spectraLikeness(experimental, get_averagine(experimental_avg, center_place, variance_parameter))
    if minus < plus:
        shift = -2
        score = minus
    else:
        shift = 2
        score = plus

    if score_best >= score:
        score_best = score
        shift_best = shift - np.sign(shift)
        while True:
            score = spectraLikeness(experimental, get_averagine(experimental_avg, center_place, shift*variance_parameter))
            if score_best >= score:
                score_best = score
                shift_best = shift
                shift += np.sign(shift)
            else:
                break

    return get_averagine(experimental_avg, center_place, shift_best*variance_parameter), score_best


def spectraLikeness(experimental, simulated, method='wasserstein'):
    '''
    Function that collects measures that can be used to compare spectra (comparison between experimental and simulated). Envemind uses default method and it can be changed only here. Function for future development.
    '''
    if simulated is None:
        return np.Inf

    if method == 'wasserstein':
        return experimental.wassersteinDistance(simulated)
    if method == 'masserstein':
        simulated = masserstein.Spectrum(confs = [*zip(simulated.np_masses(), simulated.np_probs())])
        experimental = masserstein.Spectrum(confs = [*zip(experimental.np_masses(), experimental.np_probs())])
        score = masserstein.estimate_proportions(experimental, [simulated], MTD=0.05, progress=False)['proportions'][0]
        return 1/score
    if method == 'massersteinB':
        simulated = masserstein.Spectrum(confs = [*zip(simulated.np_masses(), simulated.np_probs())])
        experimental = masserstein.Spectrum(confs = [*zip(experimental.np_masses(), experimental.np_probs())])
        score = masserstein.estimate_proportions2(experimental, [simulated], MTD=0.05, MTD_th=0.05, noise='in_both_alg2', progress=False)['proportions'][0]
        return 1/score
