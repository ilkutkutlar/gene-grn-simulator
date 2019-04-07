import matplotlib.image as image
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QVBoxLayout, QPushButton, QWidget, QMessageBox, QTabWidget

from constraint_satisfaction.constraint_satisfaction import ConstraintSatisfaction
from network_visualiser import NetworkVisualiser
from simulation.ode_simulator import OdeSimulator
from ui.constraint_satisfaction.constraints_tab import ConstraintsTab
from ui.constraint_satisfaction.mutables_tab import MutablesTab
from ui.gene_presenter import GenePresenter
from ui.simulation.deterministic_simulation_dialog import DeterministicSimulationDialog


class ReverseEngineeringModifyTab(QWidget):
    def __init__(self):
        super().__init__()
        self.main_layout = QVBoxLayout()

        self.tabs = QTabWidget()
        self.tabs.addTab(MutablesTab(), "Mutables")
        self.tabs.addTab(ConstraintsTab(), "Constraints")

        self.main_layout.addWidget(self.tabs)
        self._init_run_button()
        self.setLayout(self.main_layout)

    def _run_button_click_handler(self):

        def handler(s):
            g = GenePresenter.get_instance()
            schedule = ConstraintSatisfaction.generate_schedule(100)
            t = ConstraintSatisfaction.find_network(g.network, s,
                                                    g.get_mutables(), g.get_constraints(),
                                                    schedule)

            if t:
                im = NetworkVisualiser.visualise_as_image(t, "gene")
                plt.figure()
                plt.imshow(im)
                plt.show()

                OdeSimulator.visualise(t, s, OdeSimulator.simulate(t, s))
            else:
                error_message = QMessageBox()
                error_message.setIcon(QMessageBox.Warning)
                error_message.setWindowTitle("Error")
                error_message.setStandardButtons(QMessageBox.Ok)
                error_message.setText("No matching network found within the given parameters.")
                error_message.exec_()
                error_message.show()
                print("Error")

        DeterministicSimulationDialog(handler)

    def _init_run_button(self):
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._run_button_click_handler)
        self.main_layout.addWidget(self.run_button)
