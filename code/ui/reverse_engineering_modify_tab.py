# 100 - y => min is 100
# y - 100 => max is 100

from PyQt5.QtWidgets import QVBoxLayout, QListWidget, QHBoxLayout, QPushButton, QWidget, QMessageBox

# m = Mutable(0.5, 50, 0.5, "one")
# t = ReverseEngineering.find_network(net, s, {"rate": m}, [c1, c2], {z: (100 - z) for z in range(0, 101)})
# ode.visualise(ode.simulate())
from models.simulation_settings import SimulationSettings
from reverse_engineering.reverse_engineering import ReverseEngineering
from simulation.ode_simulator import OdeSimulator
from ui.add_constraint_dialog import AddConstraintDialog
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

        for m in GeneController.get_instance().get_mutables():
            self.mutables_list.addItem(m + ": " +
                                       str(GeneController.get_instance().get_mutables()[m]))

    def _update_constraints_list(self):
        self.constraints_list.clear()

        for c in GeneController.get_instance().get_constraints():
            self.constraints_list.addItem(str(c))

    def _add_mutable_clicked(self):
        dia = AddMutableDialog()
        dia.finished.connect(self._update_mutables_list)
        dia.exec_()

    def _add_constraint_clicked(self):
        dia = AddConstraintDialog()
        dia.finished.connect(self._update_constraints_list)
        dia.exec_()

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
        self._update_mutables_list()

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
        self._update_constraints_list()


    def _run_button_click_handler(self):
        g = GeneController.get_instance()
        s = SimulationSettings(0, 100, 100, ["x", "y", "z"])
        t = ReverseEngineering.find_network(g.network,
                                            s, g.get_mutables(), g.get_constraints(),
                                            {z: (100 - z) for z in range(0, 101)})
        # if t:
        ode = OdeSimulator(g.network, s)
        ode.visualise(ode.simulate())
        # else:
        #     error_message = QMessageBox()
        #     error_message.setIcon(QMessageBox.Warning)
        #     error_message.setWindowTitle("Error")
        #     error_message.setStandardButtons(QMessageBox.Ok)
        #     error_message.setText("No matching network found within the given parameters.")
        #     error_message.show()
        #     print("Error")

    def _init_run_button(self):
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._run_button_click_handler)
        self.main_layout.addWidget(self.run_button)
