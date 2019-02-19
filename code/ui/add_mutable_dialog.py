from PyQt5.QtWidgets import QDialog, QVBoxLayout, QFormLayout, QLineEdit, QLabel, QPushButton

from reverse_engineering.mutable import Mutable
from ui.gene_controller import GeneController


class AddMutableDialog(QDialog):

    def __init__(self):
        super().__init__()
        self.layout = QVBoxLayout()
        self._init_fields()
        self._init_add_button()
        self.setLayout(self.layout)

    def _init_add_button(self):
        button = QPushButton("Add")
        button.clicked.connect(self._add_clicked_handler)
        self.layout.addWidget(button)

    def _add_clicked_handler(self):
        s = self.species.text().strip()

        lo = self.lb.text().strip()
        lo = float(lo) if lo != "" else 0.0

        u = self.ub.text().strip()
        u = float(u) if u != "" else 0.0

        i = self.inc.text().strip()
        i = float(i) if i != "" else 0.0

        r = self.reaction.text().strip()

        GeneController.get_instance().mutables[s] = Mutable(lo, u, i, r)
        self.close()

    def _init_fields(self):
        fields = QFormLayout()

        self.species = QLineEdit()
        fields.addRow(QLabel("Mutable name"), self.species)

        self.reaction = QLineEdit()
        fields.addRow(QLabel("Reaction name (if applicable)"), self.reaction)

        self.lb = QLineEdit()
        fields.addRow(QLabel("Lower bound"), self.lb)

        self.ub = QLineEdit()
        fields.addRow(QLabel("Upper bound"), self.ub)

        self.inc = QLineEdit()
        fields.addRow(QLabel("Increments"), self.inc)

        self.layout.addLayout(fields)

