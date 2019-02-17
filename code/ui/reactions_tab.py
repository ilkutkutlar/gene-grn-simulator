from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QListWidget, QLabel, QGridLayout, QScrollArea, QHBoxLayout, QPushButton, \
    QVBoxLayout


class ReactionsTab(QWidget):
    def __init__(self):
        super().__init__()
        self.main_layout = QVBoxLayout()
        self.grid = QGridLayout()

        self._init_reactions_list()
        self._init_reaction_details()
        self._init_buttons()

        self.grid.addWidget(self.reactions_list, 0, 0)
        self.grid.addWidget(self.aux, 0, 1)

        self.main_layout.addLayout(self.grid)
        self.main_layout.addLayout(self.buttons_layout)
        self.setLayout(self.main_layout)

    def _init_reaction_details(self):
        self.aux = QScrollArea()
        self.reaction_details = QLabel()
        self.reaction_details.setAlignment(Qt.AlignLeft)
        self.reaction_details.setMinimumWidth(200)
        self.reaction_details.setWordWrap(True)
        self.aux.setWidget(self.reaction_details)

    def _init_buttons(self):
        self.buttons_layout = QHBoxLayout()
        self.add_button = QPushButton("Add")
        # self.add_button.clicked.connect(self._add_species_click_handler)
        self.remove_button = QPushButton("Remove")
        # self.remove_button.clicked.connect(self._remove_species_click_handler)
        self.buttons_layout.addWidget(self.add_button)
        self.buttons_layout.addWidget(self.remove_button)

    def _init_reactions_list(self):
        self.reactions_list = QListWidget()
        self.reactions_list.setMinimumWidth(200)
        self.reactions_list.setMaximumWidth(200)
