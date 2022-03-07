# envemind
An algorithm for monoisotopic mass prediction for high-resolution mass spectra.

# Installation
To install our package you need following packages first (available on PIP): `numpy`, `scipy`, `IsoSpecPy` and `pyteomics`. Then, just run
```
git clone git@github.com:PiotrRadzinski/envemind.git
```
and
```
pip install envemind --user
```

# How to use
You can simply run prediction by use of `monoisotopic_mass_prediction(masses, intensities)`, where masses and intensities (in daltons) are two numpy arrays. If intensities are m/z, you should also provide charge of the spectrum as third parameter e.g. `monoisotopic_mass_prediction(masses, intensities, 8)`. Remember, that the spectrum given to the function should contain single molecule (deconvoluted), and you should provide only a short range of intensities where isotopic envelope appear.

# Citing
soon
