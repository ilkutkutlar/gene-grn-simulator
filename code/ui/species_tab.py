from PyQt5.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QListWidget

from ui.gene_controller import GeneController


class SpeciesTab(QWidget):
    def __init__(self):
        super().__init__()
        self.main_layout = QVBoxLayout()

        self.species_list = QListWidget()
        self._populate_list()
        self.buttons_layout = QHBoxLayout()
        self._init_buttons()

        self.main_layout.addWidget(self.species_list)
        self.main_layout.addLayout(self.buttons_layout)
        self.setLayout(self.main_layout)

    def _init_buttons(self):
        self.add_button = QPushButton()
        self.add_button.setText("Add")
        self.remove_button = QPushButton()
        self.remove_button.setText("Remove")
        self.buttons_layout.addWidget(self.add_button)
        self.buttons_layout.addWidget(self.remove_button)

    def _populate_list(self):
        for s in GeneController.get_instance().get_species():
            self.species_list.addItem(s)

    def _add_species_click_handler(self):
        pass

    def _remove_species_click_handler(self):
        pass