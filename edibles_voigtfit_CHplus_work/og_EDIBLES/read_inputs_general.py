#================================ USER input ================================ 
# which star you are looking?
star = "HD ...."

#wave and flux for spectra
wave = []
flux = []




# ------------------ continuum fitting parameters ------------------
# absorption range is the region where continuum is not present i.e the region where the absorption lines are present. So this rigion won't be fitted for continuum.
absorption_range = ()

# Degree of Chebyshev polynomial for continuum fitting
degree = 3 
#-------------------------------------------------------------------


# Initial guess parameters for the molecules to fit the spectra  (b: gaussian width, N: column density, v_rad: radial velocity)
species_params =  {}
'''
species_params = {
                0: {  # species index
                'lambda': [list of wavelengths],
                'f': [list of f values],
                'gamma': [list of gamma values],
                'b': [list of b values for each component],
                'N': [list of N values for each component],
                'v_rad': [list of radial velocities for each component]
                },
                1: { ... },  # next species
                ...
                }
'''


















############### THIS IS FOR PLOTTING -- NO NEED TO CHANGE ANYTHING BY THE USER ###########
from astrovoigtfit import astro_voigt_fit, fit_continuum
import matplotlib.pyplot as plt

# before running lets check if the continuum fitted properly or not 
normalized_flux, continuum, poly = fit_continuum(wavelength=wave, flux=flux, absorption_range= absorption_range, degree=degree, return_std=False)


plt.figure(figsize=(8, 8))
plt.subplot(2, 1, 1)
plt.plot(wave, flux, 'b-', label='Original Flux')
plt.plot(wave, continuum, 'r-', label='Continuum Fit')
plt.axvspan(absorption_range[0], absorption_range[1], color='gray', alpha=0.3, label='Absorption Region')
plt.xlabel('Wave')
plt.ylabel('Flux')
plt.legend()
plt.title('Continuum Fitting')

# plotting continuu, Normalized flux
plt.subplot(2, 1, 2)
plt.plot(wave, normalized_flux, 'g-', label='Normalized Flux')
plt.axhline(1.0, color='k', linestyle='--', label='Continuum = 1.0')
plt.axvspan(absorption_range[0], absorption_range[1], color='gray', alpha=0.3)
plt.xlabel('Wavelength')
plt.ylabel('Normalized Flux')
plt.legend()
plt.title('Normalized Spectrum')
plt.show()


user_analysis = input("Do you want to run the fitting model? (yes/no): ").strip().lower()
    
if user_analysis == 'yes' or user_analysis == 'y':
    print("running the fitting model")
    fitresult = astro_voigt_fit(wavegrid=wave, ydata=normalized_flux, species_params=species_params,v_resolution=v_resolution, n_step=25, std_dev=std_dev)
    fitresult.params.pretty_print()
    print("chi-square value ", fitresult.chisqr)
    print("reduced chi-square value ", fitresult.redchi)
    print("FITTING RESULT :", fitresult.success)
    
    # plt.subplot(3, 1, 3)
    plt.plot(wave,fitresult.best_fit,color ='purple',label ="fit")
    plt.plot(wave, normalized_flux,color ='gray',label ='data',alpha = 0.7)
    plt.xlabel("Wavelength ($\AA$)")
    plt.ylabel("Normalised flux")
    plt.title("Multi cloud single line model fit for CH+",color = 'darkgreen')
    plt.grid()
    plt.legend()
    plt.show()
    
else: 
    print("Please rerun using your desired parameters.")