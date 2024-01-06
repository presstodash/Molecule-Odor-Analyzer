from Baza.Funkcije import SearchDescriptor, SearchSmiles, SearchName
import customtkinter as ctk
import tkinter as tk
from tkinter import messagebox
import psycopg2


class ChemicalSearchApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Odor molecule analyzer")
        self.root.geometry("400x200")

        self.category_var = tk.StringVar()
        self.search_var = tk.StringVar()

        self.create_widgets()

    def create_widgets(self):

        frame = ctk.CTkFrame(self.root)
        frame.pack(expand=True, fill="both")

        # Dropdown menu for categories
        category_label = ctk.CTkLabel(frame, text="Select Category:")
        category_label.grid(row=0, column=0, pady=10, padx=10, sticky=tk.W)

        categories = ["Descriptor", "Smiles", "Name"]

        # Create the CTkComboBox
        category_dropdown = ctk.CTkComboBox(frame, values=categories)
        category_dropdown.grid(row=0, column=1, pady=10, padx=10)
        category_dropdown.configure(state="readonly")

        # Entry
        search_label = ctk.CTkLabel(frame, text="Enter Search Value:")
        search_label.grid(row=1, column=0, pady=10, padx=10, sticky=tk.W)

        search_entry = ctk.CTkEntry(frame, textvariable=self.search_var)
        search_entry.grid(row=1, column=1, pady=10, padx=10)

        # Search button
        search_button = ctk.CTkButton(frame, text="Search", command=self.perform_search)
        search_button.grid(row=2, column=0, columnspan=2, pady=10)

        # Trace changes in the CTkComboBox
        #def category_change():
         #   self.category_var.set(category_dropdown.get())

        #category_dropdown.set(categories[0])
        #category_dropdown.bind("<<ComboboxSelected>>", category_change)

        category_dropdown.bind("<Configure>", lambda event, var=self.category_var: var.set(category_dropdown.get()))

    def perform_search(self):
        category = self.category_var.get().strip()
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

        except Exception as e:
            messagebox.showerror("Error", f"An error occurred: {str(e)}")


if __name__ == "__main__":
    root = tk.Tk()
    app = ChemicalSearchApp(root)
    root.mainloop()
