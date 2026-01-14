import tkinter as tk
from tkinter import ttk, messagebox, scrolledtext
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import sys
from io import StringIO

class SpectralFittingGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Spectral Fitting Tool - EDIBLES DR4")
        self.root.geometry("1400x900")
        self.root.configure(bg='#2b2b2b')
        
        # Style configuration
        style = ttk.Style()
        style.theme_use('clam')
        style.configure('TFrame', background='#2b2b2b')
        style.configure('TLabel', background='#2b2b2b', foreground='#ffffff', font=('Arial', 10))
        style.configure('TButton', font=('Arial', 10, 'bold'), padding=10)
        style.configure('Header.TLabel', font=('Arial', 12, 'bold'), foreground='#4a9eff')
        
        # Main container
        main_container = ttk.Frame(root, padding="10")
        main_container.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Configure grid weights
        root.columnconfigure(0, weight=1)
        root.rowconfigure(0, weight=1)
        main_container.columnconfigure(1, weight=1)
        main_container.rowconfigure(0, weight=1)
        
        # Left panel for inputs
        self.create_input_panel(main_container)
        
        # Right panel for plots
        self.create_plot_panel(main_container)
        
        # Bottom console
        self.create_console(main_container)
        
    def create_input_panel(self, parent):
        input_frame = ttk.Frame(parent, padding="10")
        input_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S), padx=(0, 10))
        
        # Scrollable frame for inputs
        canvas = tk.Canvas(input_frame, bg='#2b2b2b', highlightthickness=0, width=400)
        scrollbar = ttk.Scrollbar(input_frame, orient="vertical", command=canvas.yview)
        scrollable_frame = ttk.Frame(canvas)
        
        scrollable_frame.bind(
            "<Configure>",
            lambda e: canvas.configure(scrollregion=canvas.bbox("all"))
        )
        
        canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
        canvas.configure(yscrollcommand=scrollbar.set)
        
        canvas.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        
        # Title
        title = ttk.Label(scrollable_frame, text="SPECTRAL FITTING PARAMETERS", style='Header.TLabel')
        title.grid(row=0, column=0, columnspan=2, pady=(0, 20))
        
        row = 1
        
        # Star name
        ttk.Label(scrollable_frame, text="Star Name:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.star_var = tk.StringVar(value="HD 183143")
        ttk.Entry(scrollable_frame, textvariable=self.star_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        # File number
        ttk.Label(scrollable_frame, text="File Number:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.file_no_var = tk.IntVar(value=0)
        ttk.Entry(scrollable_frame, textvariable=self.file_no_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        # Molecules
        ttk.Label(scrollable_frame, text="Molecules (comma-separated):").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.molecules_var = tk.StringVar(value="13CH+_4032,12CH+_4032")
        ttk.Entry(scrollable_frame, textvariable=self.molecules_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        # Wavelength range section
        ttk.Label(scrollable_frame, text="WAVELENGTH RANGE", style='Header.TLabel').grid(row=row, column=0, columnspan=2, pady=(20, 10))
        row += 1
        
        ttk.Label(scrollable_frame, text="Min Wavelength:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.wave_min_var = tk.DoubleVar(value=4231.5)
        ttk.Entry(scrollable_frame, textvariable=self.wave_min_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        ttk.Label(scrollable_frame, text="Max Wavelength:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.wave_max_var = tk.DoubleVar(value=4233.5)
        ttk.Entry(scrollable_frame, textvariable=self.wave_max_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        # Absorption range
        ttk.Label(scrollable_frame, text="ABSORPTION RANGE", style='Header.TLabel').grid(row=row, column=0, columnspan=2, pady=(20, 10))
        row += 1
        
        ttk.Label(scrollable_frame, text="Min Absorption:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.abs_min_var = tk.DoubleVar(value=4232.05)
        ttk.Entry(scrollable_frame, textvariable=self.abs_min_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        ttk.Label(scrollable_frame, text="Max Absorption:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.abs_max_var = tk.DoubleVar(value=4232.8)
        ttk.Entry(scrollable_frame, textvariable=self.abs_max_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        # Molecule parameters
        ttk.Label(scrollable_frame, text="MOLECULE PARAMETERS", style='Header.TLabel').grid(row=row, column=0, columnspan=2, pady=(20, 10))
        row += 1
        
        # Molecule 1
        ttk.Label(scrollable_frame, text="Molecule 1 - b values:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.b1_var = tk.StringVar(value="2,2")
        ttk.Entry(scrollable_frame, textvariable=self.b1_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        ttk.Label(scrollable_frame, text="Molecule 1 - N values:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.N1_var = tk.StringVar(value="1e13/70,1e13/70")
        ttk.Entry(scrollable_frame, textvariable=self.N1_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        ttk.Label(scrollable_frame, text="Molecule 1 - v_rad values:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.vrad1_var = tk.StringVar(value="-11,4")
        ttk.Entry(scrollable_frame, textvariable=self.vrad1_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        # Molecule 2
        ttk.Label(scrollable_frame, text="Molecule 2 - b values:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.b2_var = tk.StringVar(value="2,2")
        ttk.Entry(scrollable_frame, textvariable=self.b2_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        ttk.Label(scrollable_frame, text="Molecule 2 - N values:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.N2_var = tk.StringVar(value="1e13,1e13")
        ttk.Entry(scrollable_frame, textvariable=self.N2_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        ttk.Label(scrollable_frame, text="Molecule 2 - v_rad values:").grid(row=row, column=0, sticky=tk.W, pady=5)
        self.vrad2_var = tk.StringVar(value="-11,4")
        ttk.Entry(scrollable_frame, textvariable=self.vrad2_var, width=25).grid(row=row, column=1, pady=5, sticky=tk.W)
        row += 1
        
        # Buttons
        button_frame = ttk.Frame(scrollable_frame)
        button_frame.grid(row=row, column=0, columnspan=2, pady=20)
        
        check_btn = tk.Button(button_frame, text="Check Continuum", command=self.check_continuum,
                              bg='#4a9eff', fg='white', font=('Arial', 11, 'bold'),
                              padx=20, pady=10, cursor='hand2')
        check_btn.pack(side=tk.LEFT, padx=5)
        
        fit_btn = tk.Button(button_frame, text="Run Fitting", command=self.run_fitting,
                           bg='#28a745', fg='white', font=('Arial', 11, 'bold'),
                           padx=20, pady=10, cursor='hand2')
        fit_btn.pack(side=tk.LEFT, padx=5)
        
    def create_plot_panel(self, parent):
        plot_frame = ttk.Frame(parent, padding="10")
        plot_frame.grid(row=0, column=1, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # Create matplotlib figure
        self.fig = Figure(figsize=(10, 8), facecolor='#2b2b2b')
        self.fig.subplots_adjust(hspace=0.3)
        
        self.ax1 = self.fig.add_subplot(211)
        self.ax2 = self.fig.add_subplot(212)
        
        for ax in [self.ax1, self.ax2]:
            ax.set_facecolor('#1e1e1e')
            ax.tick_params(colors='white')
            ax.spines['bottom'].set_color('white')
            ax.spines['top'].set_color('white')
            ax.spines['left'].set_color('white')
            ax.spines['right'].set_color('white')
            ax.xaxis.label.set_color('white')
            ax.yaxis.label.set_color('white')
            ax.title.set_color('white')
        
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas.draw()
        self.canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        
        # Toolbar
        toolbar = NavigationToolbar2Tk(self.canvas, plot_frame)
        toolbar.update()
        
    def create_console(self, parent):
        console_frame = ttk.Frame(parent, padding="10")
        console_frame.grid(row=1, column=0, columnspan=2, sticky=(tk.W, tk.E), pady=(10, 0))
        
        ttk.Label(console_frame, text="Console Output:", style='Header.TLabel').pack(anchor=tk.W)
        
        self.console = scrolledtext.ScrolledText(console_frame, height=8, bg='#1e1e1e',
                                                 fg='#00ff00', font=('Courier', 9),
                                                 insertbackground='white')
        self.console.pack(fill=tk.BOTH, expand=True)
        
    def parse_list(self, string_val):
        """Parse comma-separated string with eval for expressions like 1e13/70"""
        return [eval(x.strip()) for x in string_val.split(',')]
    
    def log_message(self, message):
        self.console.insert(tk.END, message + "\n")
        self.console.see(tk.END)
        self.root.update()
        
    def check_continuum(self):
        try:
            self.log_message("=" * 50)
            self.log_message("Checking continuum fit...")
            
            from main_run import continuum_fit
            
            star = self.star_var.get()
            molecules = [m.strip() for m in self.molecules_var.get().split(',')]
            file_no = self.file_no_var.get()
            wave_range = [self.wave_min_var.get(), self.wave_max_var.get()]
            absorption_range = (self.abs_min_var.get(), self.abs_max_var.get())
            
            wave, continuum_normalized_flux, flux, continuum = continuum_fit(
                star, molecules, file_no, wave_range, absorption_range
            )
            
            # Clear previous plots
            self.ax1.clear()
            self.ax2.clear()
            
            # Plot continuum fit
            self.ax1.plot(wave, flux, 'b-', label='Original Flux', linewidth=1.5)
            self.ax1.plot(wave, continuum, 'r-', label='Continuum Fit', linewidth=2)
            self.ax1.axvspan(absorption_range[0], absorption_range[1], 
                           color='yellow', alpha=0.2, label='Absorption Region')
            self.ax1.set_xlabel('Wavelength', color='white')
            self.ax1.set_ylabel('Flux', color='white')
            self.ax1.legend(facecolor='#2b2b2b', edgecolor='white', labelcolor='white')
            self.ax1.set_title('Continuum Fitting', color='white')
            self.ax1.grid(True, alpha=0.3, color='gray')
            
            # Plot normalized spectrum
            self.ax2.plot(wave, continuum_normalized_flux, 'g-', label='Normalized Flux', linewidth=1.5)
            self.ax2.axhline(1.0, color='cyan', linestyle='--', label='Continuum = 1.0', linewidth=2)
            self.ax2.axvspan(absorption_range[0], absorption_range[1], 
                           color='yellow', alpha=0.2)
            self.ax2.set_xlabel('Wavelength', color='white')
            self.ax2.set_ylabel('Normalized Flux', color='white')
            self.ax2.legend(facecolor='#2b2b2b', edgecolor='white', labelcolor='white')
            self.ax2.set_title('Normalized Spectrum', color='white')
            self.ax2.grid(True, alpha=0.3, color='gray')
            
            self.canvas.draw()
            
            self.log_message("✓ Continuum fit completed successfully!")
            self.log_message("Review the plots. If satisfied, click 'Run Fitting'.")
            
        except Exception as e:
            self.log_message(f"✗ Error: {str(e)}")
            messagebox.showerror("Error", f"An error occurred:\n{str(e)}")
    
    def run_fitting(self):
        try:
            self.log_message("=" * 50)
            self.log_message("Running spectral fitting model...")
            
            from main_run import astrovoigtfit_run
            
            star = self.star_var.get()
            molecules = [m.strip() for m in self.molecules_var.get().split(',')]
            file_no = self.file_no_var.get()
            wave_range = [self.wave_min_var.get(), self.wave_max_var.get()]
            absorption_range = (self.abs_min_var.get(), self.abs_max_var.get())
            
            # Parse parameters
            species_params = {
                0: {
                    'b': self.parse_list(self.b1_var.get()),
                    'N': self.parse_list(self.N1_var.get()),
                    'v_rad': self.parse_list(self.vrad1_var.get())
                },
                1: {
                    'b': self.parse_list(self.b2_var.get()),
                    'N': self.parse_list(self.N2_var.get()),
                    'v_rad': self.parse_list(self.vrad2_var.get())
                }
            }
            
            # Redirect stdout to capture print statements
            old_stdout = sys.stdout
            sys.stdout = StringIO()
            
            astrovoigtfit_run(star, molecules, wave_range, species_params, 
                            absorption_range, file_no)
            
            # Get output
            output = sys.stdout.getvalue()
            sys.stdout = old_stdout
            
            self.log_message(output)
            self.log_message("✓ Fitting completed successfully!")
            
            messagebox.showinfo("Success", "Fitting completed! Check console for results.")
            
        except Exception as e:
            sys.stdout = old_stdout
            self.log_message(f"✗ Error: {str(e)}")
            messagebox.showerror("Error", f"An error occurred:\n{str(e)}")

if __name__ == "__main__":
    root = tk.Tk()
    app = SpectralFittingGUI(root)
    root.mainloop()