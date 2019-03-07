import re

from PyQt5.QtWidgets import QDialog, QLabel, QFormLayout, QLineEdit, QPushButton, QVBoxLayout

from constraint_satisfaction.constraint import Constraint
from ui.gene_controller import GeneController


class AddConstraintDialog(QDialog):

    def __init__(self):
        super().__init__()

        self.layout = QVBoxLayout()

        self._init_fields()
        self._init_button()

        self.setLayout(self.layout)

    def _button_clicked(self):
        constraint = self.constraint.text().strip()
        m = re.match("(.*)(<= | >=)(.*)", constraint)

        if m:
            species = m.group(1).strip()
            sign = m.group(2).strip()
            value = float(m.group(3).strip())
        else:
            return      # TODO: Error!

        time = self.time.text().strip()
        time = time.split("-")
        t0 = float(time[0].strip())
        t1 = float(time[1].strip())

        if sign == "<=":
            cons = lambda v: v - value
        elif sign == ">=":
            cons = lambda v: value - v
        else:
            # TODO: error
            cons = lambda v: v

        c = Constraint(species, cons, (t0, t1))
        c.pretty_print = species + sign + str(value) + " for time: " + str(t0) + "s - " + str(t1) + "s"
        GeneController.get_instance().add_constraint(c)
        self.close()

    def _init_button(self):
        self.button = QPushButton("Add Constraint")
        self.button.clicked.connect(self._button_clicked)
        self.layout.addWidget(self.button)

    def _init_fields(self):
        fields = QFormLayout()

        self.constraint = QLineEdit()
        fields.addRow(QLabel("Constraint (e.g. x >= 10)"), self.constraint)

        self.time = QLineEdit()
        fields.addRow(QLabel("Time period (e.g. 0 - 20)"), self.time)

        self.layout.addLayout(fields)
