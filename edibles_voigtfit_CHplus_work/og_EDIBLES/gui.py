import tkinter as tk
from tkinter import ttk, scrolledtext, messagebox
import matplotlib
matplotlib.use('TkAgg')  # Set backend before importing pyplot
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import io
from contextlib import redirect_stdout, redirect_stderr
import ast
import threading
import queue

try:
    from main_run import astrovoigtfit_run, continuum_fit
except ImportError:
    print("Warning: Could not import main_run module. Make sure it's in your Python path.")

class AstronomyGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Astronomy Voigt Fitting GUI")
        self.root.geometry("1400x800")
        
        # Set up proper cleanup on window close
        self.root.protocol("WM_DELETE_WINDOW", self.on_closing)
        
        # Queue for thread communication
        self.result_queue = queue.Queue()

        # Configure main grid
        root.rowconfigure(0, weight=1)
        root.columnconfigure(0, weight=1)
        main_frame = ttk.Frame(root, padding=10)
        main_frame.grid(row=0, column=0, sticky="nsew")
        main_frame.columnconfigure(0, weight=0)
        main_frame.columnconfigure(1, weight=1)
        main_frame.rowconfigure(0, weight=1)

        # Left panel for input
        self.create_input_panel(main_frame)

        # Right panel for plot + console
        right_frame = ttk.Frame(main_frame)
        right_frame.grid(row=0, column=1, sticky="nsew", padx=(10, 0))  # Added padding
        right_frame.columnconfigure(0, weight=1)
        right_frame.rowconfigure(0, weight=1)  # plot
        right_frame.rowconfigure(1, weight=1)  # console
        
        # Configure minimum sizes to ensure both are visible
        right_frame.grid_rowconfigure(0, minsize=300)  # Minimum height for plot
        right_frame.grid_rowconfigure(1, minsize=200)  # Minimum height for console

        self.create_plot_panel(right_frame)
        self.create_console_panel(right_frame)

        self.wave = None
        self.continuum_normalized_flux = None
        self.flux = None
        self.continuum = None
        
        # Flag to track if astrovoigtfit is running
        self.astrovoigt_running = False

    def create_input_panel(self, parent):
        input_frame = ttk.LabelFrame(parent, text="Parameters", padding=10)
        input_frame.grid(row=0, column=0, sticky="ns")
        input_frame.columnconfigure(1, weight=1)
        row = 0

        ttk.Label(input_frame, text="Star:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.star_var = tk.StringVar(value="HD 183143")
        ttk.Entry(input_frame, textvariable=self.star_var, width=30).grid(row=row, column=1, sticky="ew", pady=2)
        row += 1

        ttk.Label(input_frame, text="Molecules:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.molecules_var = tk.StringVar(value="['13CH+_4032','12CH+_4032']")
        ttk.Entry(input_frame, textvariable=self.molecules_var, width=30).grid(row=row, column=1, sticky="ew", pady=2)
        row += 1

        ttk.Label(input_frame, text="File Number:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.file_no_var = tk.StringVar(value="0")
        ttk.Entry(input_frame, textvariable=self.file_no_var, width=30).grid(row=row, column=1, sticky="ew", pady=2)
        row += 1

        ttk.Label(input_frame, text="Wave Range:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.wave_range_var = tk.StringVar(value="[4231.5, 4233.5]")
        ttk.Entry(input_frame, textvariable=self.wave_range_var, width=30).grid(row=row, column=1, sticky="ew", pady=2)
        row += 1

        ttk.Label(input_frame, text="Absorption Range:").grid(row=row, column=0, sticky=tk.W, pady=2)
        self.absorption_range_var = tk.StringVar(value="(4232.05, 4232.8)")
        ttk.Entry(input_frame, textvariable=self.absorption_range_var, width=30).grid(row=row, column=1, sticky="ew", pady=2)
        row += 1

        ttk.Label(input_frame, text="Species Parameters:").grid(row=row, column=0, sticky="nw", pady=2)
        self.species_params_text = tk.Text(input_frame, height=15, width=40)
        self.species_params_text.grid(row=row, column=1, sticky="ew", pady=2)
        default_params = """{
    0: {
        'b': [2, 2],
        'N': [1.4285714285714287e11, 1.4285714285714287e11],
        'v_rad': [-11, 4]
    },
    1: {
        'b': [2, 2],
        'N': [1e13, 1e13],
        'v_rad': [-11, 4]
    }
}"""
        self.species_params_text.insert('1.0', default_params)
        row += 1

        button_frame = ttk.Frame(input_frame)
        button_frame.grid(row=row, column=0, columnspan=2, pady=10)
        ttk.Button(button_frame, text="Continuum Fitting", command=self.run_continuum_fit).pack(side=tk.LEFT, padx=5)
        self.astrovoigt_button = ttk.Button(button_frame, text="Astrovoigtfit Fitting", command=self.run_astrovoigtfit_fit)
        self.astrovoigt_button.pack(side=tk.LEFT, padx=5)
        ttk.Button(button_frame, text="Clear Console", command=self.clear_console).pack(side=tk.LEFT, padx=5)

    def create_plot_panel(self, parent):
        plot_frame = ttk.LabelFrame(parent, text="Plots", padding=5)
        plot_frame.grid(row=0, column=0, sticky="nsew", pady=(0, 2))  # Reduced padding
        plot_frame.columnconfigure(0, weight=1)
        plot_frame.rowconfigure(0, weight=1)
        
        # Smaller figure size to fit better
        self.fig, (self.ax1, self.ax2) = plt.subplots(2, 1, figsize=(8, 4))
        self.fig.tight_layout(pad=2.0)
        self.canvas = FigureCanvasTkAgg(self.fig, plot_frame)
        self.canvas.get_tk_widget().grid(row=0, column=0, sticky="nsew")

    def create_console_panel(self, parent):
        console_frame = ttk.LabelFrame(parent, text="Console Output", padding=5)
        console_frame.grid(row=1, column=0, sticky="nsew", pady=(2, 0))  # Reduced padding
        console_frame.columnconfigure(0, weight=1)
        console_frame.rowconfigure(0, weight=1)
        
        # Smaller initial height but still scrollable
        self.console_text = scrolledtext.ScrolledText(console_frame, height=8, width=50)
        self.console_text.grid(row=0, column=0, sticky="nsew")

    def on_closing(self):
        """Handle GUI closing properly"""
        try:
            if self.astrovoigt_running:
                if messagebox.askokcancel("Quit", "AstroVoigtFit is still running. Force quit?"):
                    plt.close('all')  # Close all matplotlib figures
                    self.root.quit()
                    self.root.destroy()
                else:
                    return
            else:
                plt.close('all')  # Close all matplotlib figures
                self.root.quit()
                self.root.destroy()
        except Exception as e:
            print(f"Error during cleanup: {e}")
            self.root.destroy()

    def log_to_console(self, message):
        self.console_text.insert(tk.END, message + "\n")
        self.console_text.see(tk.END)
        self.root.update_idletasks()

    def clear_console(self):
        self.console_text.delete(1.0, tk.END)

    def get_parameters(self):
        try:
            star = self.star_var.get().strip()
            molecule = ast.literal_eval(self.molecules_var.get().strip())
            file_no = int(self.file_no_var.get().strip())
            wave_range = ast.literal_eval(self.wave_range_var.get().strip())
            absorption_range = ast.literal_eval(self.absorption_range_var.get().strip())
            species_params_raw = ast.literal_eval(self.species_params_text.get('1.0', tk.END).strip())

            def convert_value(v):
                if isinstance(v, str):
                    try:
                        if '/' in v:
                            return eval(v, {'__builtins__': None})
                        return float(v) if '.' in v or 'e' in v.lower() else int(v)
                    except:
                        return v
                return v

            species_params = {}
            for key, value in species_params_raw.items():
                species_key = int(key)
                species_value = {}
                for param_name, param_list in value.items():
                    converted_list = [convert_value(x) for x in param_list]
                    species_value[param_name] = converted_list
                species_params[species_key] = species_value

            return star, molecule, file_no, wave_range, absorption_range, species_params
        except Exception as e:
            messagebox.showerror("Parameter Error", f"Error parsing parameters: {str(e)}")
            return None

    def run_continuum_fit(self):
        params = self.get_parameters()
        if params is None:
            return
        star, molecule, file_no, wave_range, absorption_range, species_params = params
        try:
            self.log_to_console("Running continuum fit...")
            f = io.StringIO()
            orig_show = plt.show
            plt.show = lambda *args, **kwargs: None
            with redirect_stdout(f), redirect_stderr(f):
                self.wave, self.continuum_normalized_flux, self.flux, self.continuum = continuum_fit(
                    star, molecule, file_no, wave_range, absorption_range
                )
            plt.show = orig_show
            output = f.getvalue()
            if output:
                self.log_to_console(output)
            self.ax1.clear()
            self.ax2.clear()
            self.ax1.plot(self.wave, self.flux, 'b-', label='Original Flux')
            self.ax1.plot(self.wave, self.continuum, 'r-', label='Continuum Fit')
            self.ax1.axvspan(absorption_range[0], absorption_range[1],
                           color='gray', alpha=0.3, label='Absorption Region')
            self.ax1.set_xlabel('Wave')
            self.ax1.set_ylabel('Flux')
            self.ax1.legend()
            self.ax1.set_title('Continuum Fitting')
            self.ax2.plot(self.wave, self.continuum_normalized_flux, 'g-', label='Normalized Flux')
            self.ax2.axhline(1.0, color='k', linestyle='--', label='Continuum = 1.0')
            self.ax2.axvspan(absorption_range[0], absorption_range[1], color='gray', alpha=0.3)
            self.ax2.set_xlabel('Wavelength')
            self.ax2.set_ylabel('Normalized Flux')
            self.ax2.legend()
            self.ax2.set_title('Normalized Spectrum')
            self.fig.tight_layout()
            self.canvas.draw()
            self.log_to_console("Continuum fit completed successfully!")
        except Exception as e:
            self.log_to_console(f"Error in continuum fit: {str(e)}")
            messagebox.showerror("Error", f"Error in continuum fit: {str(e)}")

    def run_astrovoigtfit_fit(self):
        if self.astrovoigt_running:
            messagebox.showwarning("Warning", "AstroVoigtFit is already running!")
            return
            
        params = self.get_parameters()
        if params is None:
            return
            
        # Disable button during execution
        self.astrovoigt_button.config(state='disabled', text='Running...')
        self.astrovoigt_running = True
        
        # Run in separate thread to prevent GUI freezing
        thread = threading.Thread(target=self._run_astrovoigtfit_thread, args=params)
        thread.daemon = True
        thread.start()
        
        # Start checking for completion
        self.check_astrovoigt_completion()

    def _run_astrovoigtfit_thread(self, star, molecule, file_no, wave_range, absorption_range, species_params):
        """Run astrovoigtfit in a separate thread using Agg backend"""
        try:
            # Switch to Agg backend for thread-safe operation
            import matplotlib
            original_backend = matplotlib.get_backend()
            matplotlib.use('Agg')
            
            # Create temporary file for saving plots
            import tempfile
            temp_plot_file = tempfile.NamedTemporaryFile(suffix='.png', delete=False)
            temp_plot_path = temp_plot_file.name
            temp_plot_file.close()
            
            # Capture console output
            f = io.StringIO()
            with redirect_stdout(f), redirect_stderr(f):
                # Call the original function
                from main_run import astrovoigtfit_run
                astrovoigtfit_run(star, molecule, wave_range, species_params, absorption_range, file_no)
                
                # Save the current figure (if any)
                import matplotlib.pyplot as plt
                if plt.get_fignums():
                    plt.savefig(temp_plot_path)
                    plt.close('all')
                
            # Restore original backend
            matplotlib.use(original_backend)
            
            output = f.getvalue()
            self.result_queue.put(('success', output, temp_plot_path if os.path.exists(temp_plot_path) else None))
            
        except Exception as e:
            self.result_queue.put(('error', str(e)))

    def check_astrovoigt_completion(self):
        """Check if astrovoigtfit thread has completed"""
        try:
            result = self.result_queue.get_nowait()
            result_type = result[0]
            
            if result_type == 'success':
                output, plot_path = result[1], result[2]
                self.log_to_console("AstroVoigtFit completed successfully!")
                self.log_to_console(output)
                
                # Display results in GUI
                self.ax1.clear()
                self.ax2.clear()
                
                if plot_path:
                    try:
                        # Display saved plot
                        from matplotlib import image as mpimg
                        img = mpimg.imread(plot_path)
                        self.ax1.imshow(img)
                        self.ax1.axis('off')
                        self.ax1.set_title('AstroVoigtFit Results')
                        os.unlink(plot_path)  # Clean up temp file
                    except Exception as e:
                        self.log_to_console(f"Couldn't display plot: {str(e)}")
                
                self.ax2.text(0.5, 0.5, "Check console for detailed results", 
                            ha="center", va="center", fontsize=12,
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen"))
                self.ax2.axis('off')
                
            elif result_type == 'error':
                error_msg = result[1]
                self.log_to_console(f"Error in astrovoigtfit fitting: {error_msg}")
                messagebox.showerror("Error", f"Error in astrovoigtfit fitting: {error_msg}")
                
                # Show error in plots
                self.ax1.clear()
                self.ax2.clear()
                self.ax1.text(0.5, 0.5, f"Error occurred:\n{error_msg}\n\nCheck console for details", 
                            ha="center", va="center", fontsize=12, 
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="lightcoral"))
                self.ax1.axis('off')
            
            # Re-enable button and reset flag
            self.astrovoigt_button.config(state='normal', text='Astrovoigtfit Fitting')
            self.astrovoigt_running = False
            self.fig.tight_layout()
            self.canvas.draw()
            
        except queue.Empty:
            # Still running, check again in 100ms
            self.root.after(100, self.check_astrovoigt_completion)

    # def display_astrovoigt_plot(self, buf):
    #     """Display the saved plot in the GUI"""
    #     from matplotlib import image as mpimg
        
    #     self.ax1.clear()
    #     self.ax2.clear()
        
    #     # Display the saved plot
    #     img = mpimg.imread(buf)
    #     self.ax1.imshow(img)
    #     self.ax1.axis('off')
    #     self.ax1.set_title('AstroVoigtFit Results')
        
    #     # Add some informational text
    #     self.ax2.text(0.5, 0.5, "Fit results shown above\nCheck console for parameters", 
    #                 ha="center", va="center", fontsize=12,
    #                 bbox=dict(boxstyle="round,pad=0.3", facecolor="lightgreen"))
    #     self.ax2.axis('off')
        
    #     self.fig.tight_layout()
    #     self.canvas.draw()



def main():
    root = tk.Tk()
    app = AstronomyGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()