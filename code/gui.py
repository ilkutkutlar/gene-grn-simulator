import gi

gi.require_version("Gtk", "3.0")
from gi.repository import Gtk


class GeneWindow(Gtk.Window):

    def get_add_species_button(self):
        btn = Gtk.Button()
        btn.set_label("Add Species")
        btn.set_hexpand(True)
        btn.set_valign(Gtk.Align.START)
        return btn

    def get_add_reaction_button(self):
        btn = Gtk.Button()
        btn.set_label("Add Reaction")
        btn.set_hexpand(True)
        btn.set_valign(Gtk.Align.START)
        return btn

    def get_species_list(self):
        m_list = Gtk.ListBox()
        return m_list

    def get_reactions_list(self):
        m_list = Gtk.ListBox()
        return m_list

    def get_species_container(self):
        container = Gtk.VBox()
        container.set_vexpand(False)
        container.pack_start(self.get_add_species_button(), True, True, 0)
        container.add(self.get_species_list())
        return container

    def get_reactions_container(self):
        container = Gtk.VBox()
        container.set_vexpand(False)
        container.pack_start(self.get_add_reaction_button(), True, True, 0)
        container.add(self.get_reactions_list())
        return container

    def __init__(self):
        Gtk.Window.__init__(self, title="Genes")
        self.main_container = Gtk.Box()
        self.add(self.main_container)
        self.main_container.add(self.get_species_container())
        self.main_container.add(self.get_reactions_container())


window = GeneWindow()
window.connect("destroy", Gtk.main_quit)
window.show_all()
Gtk.main()
