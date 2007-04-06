#
# inh_example.py
#
# Example code used by the epytext manual.
#
# These classes are used to illustrate the different inheritance
# styles. (grouped, listed, and included).
#
"""
Examples for the epytext manual.
"""
__docformat__='epytext'

class Animal:
    def eat(self, food): "Consume the given food object."
    def sleep(self, time): "Sleep for the given period of time."

class Bug(Animal):
    def infest(self, code): "Add new bugs to the given program."
    def hide(self):
        """
        Temporarily stop breaking a program, in order to evade
        capture.
        """
class Bird(Animal):
    def fly(self, dest): "Fly to the given destination."

class Fish(Animal):
    def swim(self, dest): "Swim to the given destination."
    
class Mammal(Animal):
    def run(self, dest): "Run to the given destination."

class Primate(Mammal):
    def climb(self, tree): "Climb up the given tree."
    def grab(self, object): "Grab hold of the given object."

class Human(Primate):
    def talk(self, animal):
        """
        Talk to the given animal.  Depending on what kind of creature
        C{animal} is, it may or may not be responsive.
        """

class Programmer(Human):
    def hack(self, code): "Improve the given program."
    def squish(self, bug, code):
        """
        Remove the given bug from the given program.
        @type bug: L{Bug}
        @param bug: The bug that should be removed from C{code}.
        """

