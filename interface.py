import tkinter as tk
from PIL import Image, ImageTk

def place_point(event):
    x, y = event.x, event.y
    # Supprimer tous les points existants
    canvas.delete("point")
    canvas.create_oval(x-3, y-3, x+3, y+3, fill="blue", outline="", tags="point")
    label.config(text=f"Clicked at ({x}, {y})")

root = tk.Tk()
root.title("Interface avec image")
screen_width = root.winfo_screenwidth()
screen_height = root.winfo_screenheight()

# Ajout d'une image
image = Image.open("imageTest.jpg")  # Insérez le chemin de votre image ici
photo = ImageTk.PhotoImage(image)

# Obtenir les dimensions de l'image
image_width = image.width
image_height = image.height

# Calculer les coordonnées du centre de l'image
center_x = screen_width // 2
center_y = screen_height // 2

# Canvas pour dessiner
canvas = tk.Canvas(root, width=screen_width-20, height=screen_height-20, bg="white")
canvas.pack()
canvas.create_image(center_x - image_width // 2, center_y - image_height // 2, anchor=tk.NW, image=photo)

# Afficher un point au centre de l'image par défaut
canvas.create_oval(center_x-3, center_y-3, center_x+3, center_y+3, fill="blue", outline="", tags="point")

# Ajout d'un label
label = tk.Label(root, text=f"Clicked at ({center_x}, {center_y})")
label.pack()

canvas.bind("<Button-1>", place_point)

# plein écran
root.attributes("-fullscreen", True)

root.mainloop()
