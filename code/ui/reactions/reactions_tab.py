from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QWidget, QListWidget, QLabel, QGridLayout, QScrollArea, QHBoxLayout, QPushButton, \
    QVBoxLayout

from ui.reactions.add_reaction_dialog import AddReactionDialog
from ui.gene_controller import GeneController


class ReactionsTab(QWidget):
    def __init__(self):
        super().__init__()
        self.main_layout = QVBoxLayout()
        self.grid = QGridLayout()

        self._init_reactions_list()
        self._init_reaction_details()
        self._init_buttons()
        self.update_list()

        self.grid.addWidget(self.reactions_list, 0, 0)
        self.grid.addWidget(self.aux, 0, 1)

        self.main_layout.addLayout(self.grid)
        self.main_layout.addLayout(self.buttons_layout)
        self.setLayout(self.main_layout)

    def _init_reaction_details(self):

        self.aux = QScrollArea()
        layout = QVBoxLayout()
        self.reaction_details = QLabel()
        self.reaction_details.setAlignment(Qt.AlignLeft)
        self.reaction_details.setMinimumWidth(200)
        self.reaction_details.setWordWrap(True)
        layout.addWidget(self.reaction_details)
        self.aux.setLayout(layout)

    def _add_reaction_click_handler(self):
        dialog = AddReactionDialog()
        dialog.finished.connect(lambda: self.update_list())
        dialog.exec_()

    def _remove_reaction_click_handler(self):
        pass

    def _reaction_list_clicked(self):
        index = self.reactions_list.currentRow()
        chosen_reaction = GeneController.get_instance().get_reactions()[index]
        self.reaction_details.setText(str(chosen_reaction))

    def _init_buttons(self):
        self.buttons_layout = QHBoxLayout()

        self.add_button = QPushButton("Add")
        self.add_button.clicked.connect(self._add_reaction_click_handler)

        self.remove_button = QPushButton("Remove")
        self.remove_button.clicked.connect(self._remove_reaction_click_handler)

        self.buttons_layout.addWidget(self.add_button)
        self.buttons_layout.addWidget(self.remove_button)

    def _init_reactions_list(self):
        self.reactions_list = QListWidget()
        self.reactions_list.setMinimumWidth(200)
        self.reactions_list.setMaximumWidth(200)
        self.reactions_list.itemClicked.connect(self._reaction_list_clicked)

    def update_list(self):
        self.reactions_list.clear()
        for s in GeneController.get_instance().get_reactions():
            self.reactions_list.addItem(s.name)
