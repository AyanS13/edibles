#================================ USER input ================================ 

# which star you are looking?
star = "HD 145502"

#interested molecules for fitting?
# molecule = ['13CH+_4032','12CH+_4032']  # List of molecules to analyze
molecule = ['12CH+_4032']

# file number from the EDIBLES DR4 data set 
file_no = 1

# wavelength range for the fit, i.e. the range between which you want to fit the lines  (*** it's a highly sensitive parameter, so choose wisely)
wave_range = [4231.65, 4233.15]



# absorption range is the region where continuum is not present i.e the region where the absorption lines are present. So this rigion won't be fitted for continuum.
absorption_range = (4232.1, 4232.5)


# Initial guess parameters for the molecules to fit the spectra  (b: gaussian width, N: column density, v_rad: radial velocity)
# species_params = {
#     0: {  # 1st molecule 13CH+
#         'b': [3, 3],  
#         'N': [1e12/30, 1e13/30],
#         'v_rad': [-12, -8]
#     },
#     1: {  # 2nd molecule 12CH+
#         'b': [3, 4],
#         'N': [1e13, 1e13],
#         'v_rad': [-12,-8]
#     }
# }

species_params = {
    0: {  # 1st molecule only 12CH+
        'b': [4, 3],  
        'N': [1e12, 1e12],
        'v_rad': [-15, -12]
    }
}








############### THIS IS FOR PLOTTING -- NO NEED TO CHANGE ANYTHING BY THE USER ###########
from main_run import astrovoigtfit_run, continuum_fit
import matplotlib.pyplot as plt

# before running lets check if the continuum fitted properly or not 
wave,continuum_normalized_flux ,flux, continuum = continuum_fit(star, molecule,file_no, wave_range, absorption_range)

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
plt.plot(wave, continuum_normalized_flux, 'g-', label='Normalized Flux')
plt.axhline(1.0, color='k', linestyle='--', label='Continuum = 1.0')
plt.axvspan(absorption_range[0], absorption_range[1], color='gray', alpha=0.3)
plt.xlabel('Wavelength')
plt.ylabel('Normalized Flux')
plt.legend()
plt.title('Normalized Spectrum')
plt.show()

# IF THE USER IS SATISFIED WITH THE NORMALIZED FLUX THEN PROCEED TO THE FINAL STEP THAT IS FITTING THE SPECTRA
user_analysis = input("Do you want to run the fitting model? (yes/no): ").strip().lower()
    
if user_analysis == 'yes' or user_analysis == 'y':
    print("running the fitting model")
    astrovoigtfit_run(star,molecule, wave_range,species_params,absorption_range,file_no)
else: 
    print("Please rerun using your desired parameters.")

