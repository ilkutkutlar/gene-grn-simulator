import matplotlib.image as image
import matplotlib.pyplot as plt
from PyQt5.QtWidgets import QVBoxLayout, QPushButton, QWidget, QMessageBox, QTabWidget, QComboBox, QGroupBox, \
    QHBoxLayout, QLabel, QLineEdit, QFormLayout

from constraint_satisfaction.constraint_satisfaction import ConstraintSatisfaction
from network_visualiser import NetworkVisualiser
from simulation.ode_simulator import OdeSimulator
from ui.constraint_satisfaction.constraints_tab import ConstraintsTab
from ui.constraint_satisfaction.mutables_tab import MutablesTab
from ui.gene_presenter import GenePresenter
from ui.simulation.deterministic_simulation_dialog import DeterministicSimulationDialog


class ConstraintSatisfactionModifyTab(QWidget):
    def __init__(self):
        super().__init__()
        self.main_layout = QVBoxLayout()

        self.tabs = QTabWidget()
        self.tabs.addTab(MutablesTab(), "Mutables")
        self.tabs.addTab(ConstraintsTab(), "Constraints")

        self.main_layout.addWidget(self.tabs)
        self.main_layout.addWidget(self._get_run_panel())
        self._method_combo_index_changed()
        self.setLayout(self.main_layout)

    def _run_button_click_handler(self):

        def handler(s):
            g = GenePresenter.get_instance()

            if self.method_combo.currentIndex() == 0:
                t = ConstraintSatisfaction.find_network(g.network, s,
                                                        g.get_mutables(), g.get_constraints())
            else:
                temperature = float(self.temperature_edit.text())
                s = ConstraintSatisfaction.generate_schedule(1000)
                t = ConstraintSatisfaction.find_closest_network(g.network, s,
                                                                g.get_mutables(), g.get_constraints(), s)

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

    def _get_annealing_panel(self):
        annealing_panel = QWidget()
        self.temperature_edit = QLineEdit()
        self.temperature_edit
        annealing_layout = QFormLayout()
        annealing_layout.addRow(QLabel("Temperature"), self.temperature_edit)
        annealing_panel.setLayout(annealing_layout)
        annealing_panel.setVisible(False)
        return annealing_panel

    def _get_best_first_panel(self):
        best_first_panel = QWidget()
        self.give_up_edit = QLineEdit()
        best_first_layout = QFormLayout()
        best_first_layout.addRow(QLabel("Give up after"), self.give_up_edit)
        best_first_panel.setLayout(best_first_layout)
        best_first_panel.setVisible(False)
        return best_first_panel

    def _method_combo_index_changed(self):
        if self.method_combo.currentIndex() == 0:
            self.best_first_panel.setVisible(True)
            self.annealing_panel.setVisible(False)
        else:
            self.best_first_panel.setVisible(False)
            self.annealing_panel.setVisible(True)

    def _get_method_combo(self):
        self.method_combo = QComboBox()
        self.method_combo.addItems(["Best-first search (exact match)", "Simulated annealing (closest match)"])
        self.method_combo.currentIndexChanged.connect(self._method_combo_index_changed)
        return self.method_combo

    def _get_run_button(self):
        self.run_button = QPushButton("Run")
        self.run_button.clicked.connect(self._run_button_click_handler)
        return self.run_button

    def _get_run_options_panel(self):
        run_options_box = QVBoxLayout()
        mode_layout = QFormLayout()
        mode_layout.addRow(QLabel("Search method"), self._get_method_combo())
        run_options_box.addLayout(mode_layout)
        self.annealing_panel = self._get_annealing_panel()
        self.best_first_panel = self._get_best_first_panel()
        run_options_box.addWidget(self.annealing_panel)
        run_options_box.addWidget(self.best_first_panel)
        return run_options_box

    def _get_run_panel(self):
        run_group_box = QGroupBox()
        run_panel = QVBoxLayout()
        run_panel.addLayout(self._get_run_options_panel())
        run_panel.addWidget(self._get_run_button())
        run_group_box.setLayout(run_panel)
        return run_group_box
