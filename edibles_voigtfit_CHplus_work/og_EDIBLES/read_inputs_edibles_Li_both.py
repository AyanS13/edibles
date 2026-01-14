from main_run import astrovoigtfit_run, continuum_fit
import matplotlib.pyplot as plt

#user inputsn
star =     "HD 147889" #HD 149404"

#interested molecules
molecule = ['7Li_6707','6Li_6707']  # List of molecules to analyze

# file number of the observations
file_no = 3

# wavelength range for the fit  *** it's a highly sensitive parameter, so choose wisely
wave_range = [6707,6708.5]  


# absorption range is the region where continuum is not present i.e the region where the absorption lines are present
absorption_range = (6707.45, 6708)

# guessed species parameters for the molecules to fit the spectra 
species_params = {
    0: {  # 1st molecule 
        'b': [1.5],  
        'N': [10**(10.1651)],
        'v_rad':[-8.3]
    },
    1: {  # 2nd molecule 
        'b': [1.5],
        'N': [10**(9.0765)],
        'v_rad': [-8.3]
    }
}


degree = 3

# before running lets check if the continuum fitted properly or not 
wave,continuum_normalized_flux ,flux, continuum = continuum_fit(star, molecule,file_no, wave_range, absorption_range, degree)

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


user_analysis = input("Do you want to run the fitting model? (yes/no): ").strip().lower()
    
if user_analysis == 'yes' or user_analysis == 'y':
    print("running the fitting model")
    astrovoigtfit_run(star,molecule, wave_range,species_params,absorption_range,file_no)
else: 
    print("Please rerun using your desired parameters.")