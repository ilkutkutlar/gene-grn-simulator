# 100 - y => min is 100
# y - 100 => max is 100
import re

from PyQt5.QtWidgets import QDialog, QVBoxLayout, QListWidget, QHBoxLayout, QPushButton, QFormLayout, QLabel, QLineEdit, \
    QWidget, QInputDialog

# c1 = Constraint("y", lambda y: 200 - y, (40, 60))
# c2 = Constraint("z", lambda x: x - 150, (0, 20))

# m = Mutable(0.5, 50, 0.5, "one")
# t = ReverseEngineering.find_network(net, s, {"rate": m}, [c1, c2], {z: (100 - z) for z in range(0, 101)})
# ode.visualise(ode.simulate())
from reverse_engineering.constraint import Constraint
from reverse_engineering.reverse_engineering import Mutable
from ui.add_mutable_dialog import AddMutableDialog
from ui.gene_controller import GeneController


class ReverseEngineeringModifyTab(QWidget):
    def __init__(self):
        super().__init__()
        self.main_layout = QVBoxLayout()

        self._init_mutables()
        self._init_constraints()
        self._init_run_button()
        self.setLayout(self.main_layout)

    def _update_mutables_list(self):
        self.mutables_list.clear()

        for m in GeneController.get_instance().mutables:
            self.mutables_list.addItem(m)

    def _add_mutable_clicked(self):
        dia = AddMutableDialog()
        dia.finished.connect(self._update_mutables_list)
        dia.exec_()

    def _add_constraint_clicked(self):
        dia = QDialog()

        layout = QLabel()
        fields = QFormLayout()


        (text, ok) = QInputDialog.getText(self, "Add Constraint", "Format: x > 10")

        if ok:
            m = re.match("(.*)(< | > | <= | >=)(.*)", text)
            species = m.group(1).strip()
            sign = m.group(2).strip()
            value = m.group(3).strip()

            Constraint()

            print(species)
            print(sign)
            print(value)

        self._update_list()

    def _init_mutables(self):
        self.mutables_list = QListWidget()
        mutables_buttons_layout = QHBoxLayout()
        self.add_mutable_button = QPushButton("Add Mutable")
        mutables_buttons_layout.addWidget(self.add_mutable_button)
        self.add_mutable_button.clicked.connect(self._add_mutable_clicked)

        self.remove_mutable_button = QPushButton("Remove Mutable")
        mutables_buttons_layout.addWidget(self.remove_mutable_button)

        self.main_layout.addWidget(self.mutables_list)
        self.main_layout.addLayout(mutables_buttons_layout)

    def _init_constraints(self):
        self.constraints_list = QListWidget()
        constraints_buttons_layout = QHBoxLayout()
        self.add_constraint_button = QPushButton("Add Constraint")
        constraints_buttons_layout.addWidget(self.add_constraint_button)
        self.add_constraint_button.clicked.connect(self._add_constraint_clicked)

        self.remove_constraint_button = QPushButton("Remove Constraint")
        constraints_buttons_layout.addWidget(self.remove_constraint_button)

        self.main_layout.addWidget(self.constraints_list)
        self.main_layout.addLayout(constraints_buttons_layout)

    def _run_button_click_handler(self):
        pass

    def _init_run_button(self):
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._run_button_click_handler)
        self.main_layout.addWidget(self.run_button)
