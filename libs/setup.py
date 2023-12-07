class Setup:
    def __init__(self, height, width, weight):
        self.H = height
        self.W  = width
        self.weight = weight

    def display_status(self):
        print("*------------ Setup ------------*")
        print("weight = ", self.weight)
        print("Height = ", self.H)
        print("Width  = ", self.W)
        