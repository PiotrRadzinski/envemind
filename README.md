# envemind
An algorithm for monoisotopic mass prediction for high-resolution mass spectra.

# Installation
To install our package run
```
git clone git@github.com:PiotrRadzinski/envemind.git
```
```
pip install envemind --user
```

# How to use
You can simply run prediction by use of `monoisotopic_mass_prediction(masses, intensities)`, where masses (in daltons) and intensities are two numpy arrays. If spectrum's x-axis contains m/z, you should also provide charge of the spectrum as third parameter e.g. `monoisotopic_mass_prediction(masses, intensities, 8)`. Remember, that the spectrum given to the function should contain single molecule (deconvoluted), and you should provide only a short range of masses where isotopic envelope appear.

### Example
```
import envemind as ev

spectra_dict = ev.mzMLread('data/Multinomial_600.mzML')
S = spectra_dict[0]

masses = [j[0] for j in S if j[0] > 2119 and j[0] < 2121]
intensities = [j[1] for j in S if j[0] > 2119 and j[0] < 2121]

ev.monoisotopic_mass_prediction(masses, intensities, 8)
```

# Citing
soon
