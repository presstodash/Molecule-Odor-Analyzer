from Molecule_Odor_Analyzer.Baza.Funkcije import SearchDescriptor, SearchSmiles, SearchName, getTanimoto, getScent
import customtkinter as ctk
import tkinter as tk
import CTkTable as ctktable
from tkinter import messagebox
import psycopg2
from PIL import Image, ImageTk


class ChemicalSearchApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Odor molecule analyzer")
        self.root.geometry("1200x600")

        self.category_var = tk.StringVar()
        self.search_var = tk.StringVar()

        self.create_widgets()

    def create_widgets(self):

        self.frame = ctk.CTkFrame(self.root)
        self.frame.pack(expand=True, fill="both")

        # Dynamic frame for search display
        self.frames = {"Descriptor": ctk.CTkFrame(self.frame),
                       "Smiles": ctk.CTkFrame(self.frame),
                       "Name": ctk.CTkFrame(self.frame)}

        for frame in self.frames.values():
            frame.grid(row = 1, column = 0, sticky = 'nsew')
            
        for frame in self.frames.values():
            frame.grid_forget()

        #dynamic sizing
        self.frame.grid_rowconfigure(1, weight=1)

        self.columncount = 5
        
        for i in range(self.columncount):
            self.frame.grid_columnconfigure(i, weight=1)

        
        # Dropdown menu for categories
        self.category_label = ctk.CTkLabel(self.frame, text="Select Category:")
        self.category_label.grid(row=0, column=0, pady=10, padx=10, sticky=tk.W)

        categories = ["Descriptor", "Smiles", "Name"]

        # Create the CTkComboBox
        self.category_dropdown = ctk.CTkComboBox(self.frame, values=categories)
        self.category_dropdown.grid(row=0, column=1, pady=10, padx=10)
        self.category_dropdown.configure(state="readonly")

        # Entry
        self.search_label = ctk.CTkLabel(self.frame, text="Enter Search Value:")
        self.search_label.grid(row=0, column=2, pady=10, padx=10, sticky=tk.W)

        self.search_entry = ctk.CTkEntry(self.frame, textvariable=self.search_var)
        self.search_entry.grid(row=0, column=3, pady=10, padx=10)

        # Search button
        self.search_button = ctk.CTkButton(self.frame, text="Search", command=self.perform_search)
        self.search_button.grid(row=0, column=4, columnspan=2, pady=10)

        # Trace changes in the CTkComboBox
        #self.category_var.trace('w', lambda *args: self.category_var.set(category_dropdown.get()))

        self.category_dropdown.set(categories[0])


    def perform_search(self):
        category = self.category_dropdown.get().strip()
        search_value = self.search_var.get().strip()

        if not search_value:
            messagebox.showwarning("ser")
            return

        if not category:
            messagebox.showwarning("cat")
            return

        if not category or not search_value:
            messagebox.showwarning("Warning", "Please select a category and enter a search value.")
            return

        try:
            if category == "Descriptor":
                result = SearchDescriptor(search_value)
            elif category == "Smiles":
                result, result2, img = SearchSmiles(search_value)
            elif category == "Name":
                result, result2, img = SearchName(search_value)

            if not result:
                messagebox.showinfo("Result", "No matching records found.")
            else:
                # Hide all frames
                for frame in self.frames.values():
                    frame.grid_remove()

                # Show the selected frame
                if category in ["Descriptor","Smiles","Name"]:

                    match category:
                        case "Descriptor":

                            self.frames[category].grid(columnspan=self.columncount, sticky='nsew', pady=10, padx=10,)

                            self.scrollable_frame = ctk.CTkScrollableFrame(master=self.frames[category])
                            self.scrollable_frame.pack(expand=True, fill="both")
                            
                            self.table = ctktable.CTkTable(master=self.scrollable_frame, row=len(result), column=5, values=result)
                            self.table.pack(expand=True, fill="both")
                        
                        
                        case "Smiles":
                            #messagebox.showinfo("Result", f"Matching Record(s):\n{result}")
                            
                            self.frames[category].grid(rowspan=self.columncount, columnspan=self.columncount, sticky='nsew', pady=10, padx=10,)
                                                       
                            self.molecule_display = ctk.CTkFrame(master=self.frames[category], border_width=5)
                            self.molecule_display.pack(expand=True, fill="both")
                            
                            for i in range(10):
                                self.molecule_display.grid_columnconfigure(i, weight=1)
                            
                            for i in range(6):
                                self.molecule_display.grid_rowconfigure(i, weight=1)
                                
                            self.ime_label = ctk.CTkLabel(self.molecule_display, text="Ime Molekule:")
                            self.ime_label.grid(row=0, column=0, pady=10, padx=10, sticky=tk.W)
                            
                            self.ime_value = ctk.CTkLabel(self.molecule_display, text=result[1])
                            self.ime_value.grid(row=0, column=1, pady=10, padx=10)
                            
                            
                            self.smiles_label = ctk.CTkLabel(self.molecule_display, text="SMILES Kod:")
                            self.smiles_label.grid(row=1, column=0, pady=10, padx=10, sticky=tk.W)
                            
                            self.ime_value = ctk.CTkLabel(self.molecule_display, text=result[2])
                            self.ime_value.grid(row=1, column=1, pady=10, padx=10)
                            
                            
                            self.formula_label = ctk.CTkLabel(self.molecule_display, text="Molekulska formula:")
                            self.formula_label.grid(row=2, column=0, pady=10, padx=10, sticky=tk.W)
                            
                            self.ime_value = ctk.CTkLabel(self.molecule_display, text=result[3])
                            self.ime_value.grid(row=2, column=1, pady=10, padx=10)
                            
                            
                            self.masa_label = ctk.CTkLabel(self.molecule_display, text="Relativna molekulska masa:")
                            self.masa_label.grid(row=3, column=0, pady=10, padx=10, sticky=tk.W)
                            
                            self.ime_value = ctk.CTkLabel(self.molecule_display, text=result[4])
                            self.ime_value.grid(row=3, column=1, pady=10, padx=10)
                            
                            if img is not None: 
                                self.slika_label = ctk.CTkLabel(self.molecule_display, text="Grafički prikaz:")
                                self.slika_label.grid(row=4, column=0, pady=10, padx=10, sticky=tk.W)
                                
                                img = ctk.CTkImage(light_image=img, dark_image=img, size=(300, 300))
                                label_image = ctk.CTkLabel(self.molecule_display, image=img, text="")
                                label_image.grid(row=4, column=1, pady=10, padx=10)
                 
                            #consider making dynamic number in the future
                            numberofdisplayed = 10
                            
                            tanimotoresult = getTanimoto(result[0],numberofdisplayed)

                            self.scrollable_frame = ctk.CTkScrollableFrame(master=self.molecule_display, border_width=5)
                            self.scrollable_frame.grid(row=0, column = 3, columnspan = 8, rowspan=3, sticky="nsew", pady=10, padx=10)

                            self.tanimototable = ctktable.CTkTable(master=self.scrollable_frame, row=len(tanimotoresult), column=3, values=tanimotoresult)
                            self.tanimototable.pack(expand=True, fill="both")
                            
                            
                            scentresult = getScent(result[0],numberofdisplayed)

                            self.scrollable_frame2 = ctk.CTkScrollableFrame(master=self.molecule_display, border_width=5)
                            self.scrollable_frame2.grid(row=3, column = 3, columnspan = 8, rowspan=3, sticky="nsew", pady=10, padx=10)

                            self.scenttable = ctktable.CTkTable(master=self.scrollable_frame2, row=len(scentresult), column=3, values=scentresult)
                            self.scenttable.pack(expand=True, fill="both")
                        
                        case "Name":
                            #messagebox.showinfo("Result", f"Matching Record(s):\n{result}")
                        
                            
                            
                            self.frames[category].grid(rowspan=self.columncount, columnspan=self.columncount, sticky='nsew', pady=10, padx=10,)
                            
                            self.molecule_display = ctk.CTkFrame(master=self.frames[category])
                            self.molecule_display.pack(expand=True, fill="both")
                            
                            for i in range(10):
                                self.molecule_display.grid_columnconfigure(i, weight=1)
                            
                            for i in range(6):
                                self.molecule_display.grid_rowconfigure(i, weight=1)
                                
                            self.ime_label = ctk.CTkLabel(self.molecule_display, text="Ime Molekule:")
                            self.ime_label.grid(row=0, column=0, pady=10, padx=10, sticky=tk.W)
                            
                            self.ime_value = ctk.CTkLabel(self.molecule_display, text=result[1])
                            self.ime_value.grid(row=0, column=1, pady=10, padx=10)
                            
                            
                            self.smiles_label = ctk.CTkLabel(self.molecule_display, text="SMILES Kod:")
                            self.smiles_label.grid(row=1, column=0, pady=10, padx=10, sticky=tk.W)
                            
                            self.ime_value = ctk.CTkLabel(self.molecule_display, text=result[2])
                            self.ime_value.grid(row=1, column=1, pady=10, padx=10)
                            
                            
                            self.formula_label = ctk.CTkLabel(self.molecule_display, text="Molekulska formula:")
                            self.formula_label.grid(row=2, column=0, pady=10, padx=10, sticky=tk.W)
                            
                            self.ime_value = ctk.CTkLabel(self.molecule_display, text=result[3])
                            self.ime_value.grid(row=2, column=1, pady=10, padx=10)
                            
                            
                            self.masa_label = ctk.CTkLabel(self.molecule_display, text="Relativna molekulska masa:")
                            self.masa_label.grid(row=3, column=0, pady=10, padx=10, sticky=tk.W)
                            
                            self.ime_value = ctk.CTkLabel(self.molecule_display, text=result[4])
                            self.ime_value.grid(row=3, column=1, pady=10, padx=10)
                            
                            if img is not None: 
                                self.slika_label = ctk.CTkLabel(self.molecule_display, text="Grafički prikaz:")
                                self.slika_label.grid(row=4, column=0, pady=10, padx=10, sticky=tk.W)

                                img = ctk.CTkImage(light_image=img, dark_image=img, size=(300, 300))
                                label_image = ctk.CTkLabel(self.molecule_display, image=img, text="")
                                label_image.grid(row=4, column=1, pady=10, padx=10)
                            
                            #consider making dynamic number in the future
                            numberofdisplayed = 10
                            
                            tanimotoresult = getTanimoto(result[0],numberofdisplayed)

                            self.scrollable_frame = ctk.CTkScrollableFrame(master=self.molecule_display, border_width=5)
                            self.scrollable_frame.grid(row=0, column = 3, columnspan = 8, rowspan=3, sticky="nsew", pady=10, padx=10)

                            self.tanimototable = ctktable.CTkTable(master=self.scrollable_frame, row=len(tanimotoresult), column=3, values=tanimotoresult)
                            self.tanimototable.pack(expand=True, fill="both")
                            
                            
                            scentresult = getScent(result[0],numberofdisplayed)

                            self.scrollable_frame2 = ctk.CTkScrollableFrame(master=self.molecule_display, border_width=5)
                            self.scrollable_frame2.grid(row=3, column = 3, columnspan = 8, rowspan=3, sticky="nsew", pady=10, padx=10)

                            self.scenttable = ctktable.CTkTable(master=self.scrollable_frame2, row=len(scentresult), column=3, values=scentresult)
                            self.scenttable.pack(expand=True, fill="both")
                        
                        case _:
                            print("No category found.")
                    


            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")


if __name__ == "__main__":
    root = tk.Tk()
    app = ChemicalSearchApp(root)
    root.mainloop()

