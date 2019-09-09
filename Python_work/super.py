class Rectangle:
	def __init__(self, length, width):
		self.length = length
		self.width = width

	def area(self):
		return self.length * self.width

	def perimeter(self):
		return 2 * self.length + 2 * self.width

class Square:
	def __init__(self, length):
		self.length = length

	def area(self):
		return self.length * self.length

	def perimeter(self):
		return 4 * self.length

	
"""

In this example, you have two shapes that are related to each other. A square
is a special kind of rectangle. The code, however, doesn't reflect that relationship and thus has code that is essentially repated.

By using inheritance, you can reduce the amount of code you write 
while simultaneously reflecting the real-world relationship between
rectangles and squares


"""


# Here we declare that the Square class inherits from the Rectangle class

class Square(Rectangle):
	"""
	Here, you've used super() to call the __init__() of the rectangle 
	class, allowing you to use it in the Square class without repeating 
	code.

	Below, the core functionality remains after making changes


	"""
	def __init__(self, length):
		super().__init__(length,length)


"""
Here, you've used super() to call the __init__() of the Rectangle
class, allowing you to use it in Square class without repeating code.

What can super() do for you?

Like in other object-orientated languages, it allows you 
to call methods of the superclass in your subclass. The 
primary use case of this is to extend the functionality 
of the inherited method 


A super() deep dive
------------------

Before heading into multiple inheritance, let's take a quick detour 
into the mechanics of super().

while the examples above call super() without any parameters, super() 
can also take two parameterrs - the first is the subclass, and the 
second parameter is an object that is an instance of that subclass

"""

class Triangle:
	def __init__(self, base, height):
		self.base = base
		self.height = height

	def(area, self):
		return 0.5 * self.base * self.height

class RightPyramid(Triangle, Square):
	def __init__(self, base, slant_height):
		self.base = base
		self.slant_height = slant_height

	def area(self):
		base_area = super().area()
		perimeter = super().perimeter()
		return 0.5 * perimeter * self.slant_height + base_area
	

		

