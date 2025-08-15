####################################################################
# _author_ = "Debayan Das"
# __doc__ = "Largest ever Lithium survey in ISM"
####################################################################



#======================== PACKAGES =================================
import numpy as np
import VoigtFit
from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from scipy.interpolate import CubicSpline
#===================================================================



## Target
target ='HD 113904' 

# #======================== GETTING THE SPECTRUM ======================

filename = '/HD113904/BLUE_437/HD113904_w437_blue_20141221_O15.fits'   #test[0]           #4

# wrange = [6707,6708.5]      #1
# wrange = [6707.06,6708.1]   #2
# wrange = [6706.85,6708.7]   #3
wrange = [4231.75, 4233.25] #4

print(filename)
sp = EdiblesSpectrum(filename)
sp.getSpectrum(xmin=wrange[0], xmax=wrange[1])
wave = sp.bary_wave
flux = sp.bary_flux

# Filter the spectrum within the specified wavelength ranges
idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
wave = wave[idx]
flux = flux[idx]
# Normalize the flux
flux = flux / np.median(flux)
#===================================================================


# import numpy as np
# import matplotlib.pyplot as plt
# from numpy.polynomial.chebyshev import Chebyshev
# wavelength = wave



'''
This is for continuum fitting using chebyshev polynomial.
'''
# # mask_wavelength =[6707.4, 6708]                    #1
# # mask_wavelength =[6707.38, 6707.87]                #2
# # mask_wavelength =[6707.38, 6707.9]                 #3
# mask_wavelength =[6707.38, 6708.05]                 #4


# mask = (wavelength < mask_wavelength[0]) | (wavelength > mask_wavelength[1])
# norm_wavelength = 2 * (wavelength - wavelength.min()) / (wavelength.ptp()) - 1

# degree = 2 
# cheb_fit = Chebyshev.fit(norm_wavelength[mask], flux[mask], degree)

anchor_points_wave = np.array([4231.75, 4232, 4232.25, 4232.5, 4232.75, 4233, 4233.25]) 
anchor_points_flux = np.array([1.0015, 1.0015, 1.0015, 1.0015, 1.0013, 1.0005, 1.00])


# Fit cubic splines to the anchor points
cs = CubicSpline(anchor_points_wave, anchor_points_flux) # so the cubic spline routine is not really taking all the points into account. only anchors

# Generate the continuum model
continuum = cs(wave)

# fitted_continuum = cheb_fit(norm_wavelength)
normalized_flux = flux / continuum

#======================== FITTING THE SPECTRUM ======================

z_DLA =  0
# wl =wave
# spec =flux #data

wl = wave
spec = normalized_flux
err = np.full(len(wl), 0.0015)
dataset = VoigtFit.DataSet(z_DLA)
dataset.add_data(wl=wl, flux=spec, res=3,normalized=True)
dataset.set_name(f'{target}_plot')
dataset.verbose = True
res = 3
# dataset.velspan = 50
dataset.cheb_order = -1

# -- Add the data loaded from thecd 
dataset.add_data(wl, spec, res)

# -- Define absorption lines:
dataset.add_many_lines(['12CHI_1', '13CHI_2'])

# -- Add components for each ion:

#                      ion    z         b   logN
# dataset.add_component_velocity('12CHI', 13, 3, 13)   # towards left 
# dataset.copy_components(from_ion='12CHI',to_ion='13CHI',tie_b=True)

# dataset.add_component_velocity('12CHI', -24, 3, 13)   # towards left 
# dataset.copy_components(from_ion='12CHI',to_ion='13CHI',tie_b=True)

# dataset.add_component_velocity('12CHI', -22, 3, 13)   # towards left 
# # dataset.copy_components(from_ion='12CHI',to_ion='13CHI',tie_b=True)


# dataset.add_component_velocity('12CHI', -25.3, 3, 13)   # towards left 
# dataset.copy_components(from_ion='12CHI',to_ion='13CHI',tie_b=True)

# # # # # dataset.prepare_dataset(norm=False, mask=False)

dataset.add_component_velocity('12CHI', 1, 3, 13)   # towards left 
# dataset.copy_components(from_ion='12CHI',to_ion='13CHI',tie_b=True)

# dataset.add_component_velocity('12CHI', -10, 3, 13)   # towards left 
# # dataset.copy_components(from_ion='12CHI',to_ion='13CHI',tie_b=True)


# dataset.add_component_velocity('12CHI', 15, 3, 13)   # towards left 
# # dataset.copy_components(from_ion='12CHI',to_ion='13CHI',tie_b=True)

# dataset.add_component_velocity('12CHI', 10, 3, 13)   # towards left 
# dataset.copy_components(from_ion='12CHI',to_ion='13CHI',tie_b=True)


# dataset.add_component_velocity('12CHI', -25.3, 3, 13)   # towards left 
# dataset.copy_components(from_ion='12CHI',to_ion='13CHI',tie_b=True)

# dataset.add_component_velocity('12CHI', 11, 3, 12)   # towards left 
# dataset.copy_components(from_ion='12CHI',to_ion='13CHI',tie_b=True)

dataset.prepare_dataset()

#---
import sys
import io

# Create a string stream to capture the output
captured_output = io.StringIO()

# Redirect standard output to the string stream
sys.stdout = captured_output

# Run the fit:
popt, chi2 = dataset.fit(verbose=True) # USING THIS GIVES THIS ERROR: TypeError: Improper input: func input vector length N=4 must not exceed func output vector length M=0
dataset.plot_fit(filename=dataset.name,loc = 'center', fontsize=8,legend=False,xunit='wl')
# print(dataset.best_fit['logN0_7LiI'].value)
# print(dataset.best_fit['logN1_7LiI'].value)
# print(dataset.best_fit['logN0_6LiI'].value)
# print(dataset.best_fit['logN1_6LiI'].value)
# print("12CH+/13CH+ ratio = ", 10**(dataset.best_fit['logN_12CHI_1'].value)/10**(dataset.best_fit['logN_13CHI_2'].value))
# print("7Li/6Li ration = ", 10**(dataset.best_fit['logN1_7LiI'].value)/10**(dataset.best_fit['logN1_6LiI'].value))

print("chi-squared value",chi2)

# Reset standard output back to the console
sys.stdout = sys.__stdout__

# Now you can access the captured output as a string
terminal_output = captured_output.getvalue()
#---

# dataset.save_fit_regions(filename=dataset.name)




# popt.params.pretty_print()
print("chi-squared value",chi2)
# print("12CH+/13CH+ ratio = ", 10**(dataset.best_fit['logN_12CHI_1'].value)/10**(dataset.best_fit['logN_13CHI_2'].value))
VoigtFit.io.output.plot_all_lines(filename=dataset.name,dataset=dataset)
dataset.save_fit_regions(filename=dataset.name)

import matplotlib.pyplot as plt
from PyPDF2 import PdfReader, PdfWriter
from io import BytesIO
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import letter

def append_best_fit_to_pdf(existing_pdf, output_pdf, best_fit_text):
    # Step 1: Read the existing PDF
    reader = PdfReader(existing_pdf)
    writer = PdfWriter()

    for page_num, page in enumerate(reader.pages):
        packet = BytesIO()

        # Step 2: Create a new canvas to overlay text
        can = canvas.Canvas(packet, pagesize=letter)
        can.setFont("Times-Roman", 7)
        # Positioning the text in the lower half of the page
        x, y = 10, 233  # Adjust (x, y) as needed for placement
        can.drawString(x, y, "Best Fit Parameters:")
        
        # Adding each parameter in a new line
        for i, line in enumerate(best_fit_text.splitlines()):
            can.drawString(x, y - (7 * (i + 1)), line)

        can.save()

        # Step 3: Merge overlay with the original page
        packet.seek(0)
        overlay = PdfReader(packet)
        page.merge_page(overlay.pages[0])

        # Add the modified page to the writer
        writer.add_page(page)

    # Step 4: Save the updated PDF to the output path
    with open(output_pdf, "wb") as output_file:
        writer.write(output_file)

existing_pdf_path = f'{target}_plot.pdf'


# Pass the captured output to `best_fit_text`
best_fit_text = terminal_output

append_best_fit_to_pdf(existing_pdf_path, existing_pdf_path, best_fit_text)


#------------------

# from matplotlib.backends.backend_pdf import PdfPages

# # Your plot data
# # Assuming 'wavelength', 'flux', 'fitted_continuum', and 'normalized_flux' are already defined

# # Create a PDF file
# continuum_pdf_path = f'{target}_continuum.pdf'
# with PdfPages(continuum_pdf_path) as pdf:
#     # If the first page has content, just add a new page for the plot
#     plt.figure(figsize=(20, 10))
#     plt.plot(wavelength, flux, label='Observed Flux', color='blue')
#     # plt.plot(wavelength, fitted_continuum, label='Fitted Continuum', color='orange', linestyle='--')
#     plt.plot(wavelength, normalized_flux, label='Normalized Flux', color='green')
#     # plt.axvspan(mask_wavelength[0], mask_wavelength[1], color='red', alpha=0.3, label='Masked Region')
#     plt.xlabel('Wavelength (Å)')
#     plt.ylabel('Flux')
#     plt.legend()
#     plt.title('Continuum Fitting and Normalization')

#     # Save the plot to the PDF (this will be the second page)
#     pdf.savefig()
#     plt.close()

























