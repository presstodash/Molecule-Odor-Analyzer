from Molecule_Odor_Analyzer.Baza.Funkcije import SearchDescriptor, SearchSmiles, SearchName
import customtkinter as ctk
import tkinter as tk
from tkinter import messagebox
import psycopg2


class ChemicalSearchApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Odor molecule analyzer")
        self.root.geometry("1000x600")

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
                result = SearchSmiles(search_value)
            elif category == "Name":
                result = SearchName(search_value)

            if not result:
                messagebox.showinfo("Result", "No matching records found.")
            else:
                messagebox.showinfo("Result", f"Matching Record(s):\n{result}")

            # Hide all frames
            for frame in self.frames.values():
                frame.grid_remove()

            # Show the selected frame
            if category in ["Descriptor","Smiles","Name"]:
                self.frames[category].grid(columnspan=self.columncount, sticky='nsew', pady=10, padx=10,)
            
        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")


if __name__ == "__main__":
    root = tk.Tk()
    app = ChemicalSearchApp(root)
    root.mainloop()

