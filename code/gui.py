import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk


class GeneWindow(Gtk.Window):

    def __init__(self):
        Gtk.Window.__init__(self, title="Genes")

        self.button = Gtk.Button(label="Test")
        self.button.connect("clicked", lambda w: print("Hello!"))
        self.add(self.button)


window = GeneWindow()
window.connect("destroy", Gtk.main_quit)
window.show_all()
Gtk.main()
