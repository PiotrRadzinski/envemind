# envemind
An algorithm for monoisotopic mass prediction for high-resolution mass spectra.

# Installation
```
git clone 
```



# How to use
You can simply run prediction by use of `monoisotopic_mass_prediction(masses, intensities)`, where masses and intensities (in daltons) are two numpy arrays. If intensities are m/z, you should also provide charge of the spectrum as third parameter e.g. `monoisotopic_mass_prediction(masses, intensities, 8)`. Remember, that the spectrum given to the function should contain single molecule (deconvoluted), and you should provide only a short range of intensities where isotopic envelope appear.  

# Citing
soon
