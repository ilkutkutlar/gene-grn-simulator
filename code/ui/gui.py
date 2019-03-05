from PyQt5.QtWidgets import QTabWidget, QFileDialog
from PyQt5.QtWidgets import QVBoxLayout
from PyQt5.QtWidgets import QWidget, QApplication, QMainWindow, QAction, QMessageBox

from input_output.sbml_saver import SbmlSaver
from simulation.ode_simulator import OdeSimulator
from ui.gene_controller import GeneController
from ui.open_sbml_dialog import OpenSbmlDialog
from ui.reactions.reactions_tab import ReactionsTab
from ui.reverse_engineering.reverse_engineering_modify_tab import ReverseEngineeringModifyTab
from ui.simulation.deterministic_simulation_dialog import DeterministicSimulationDialog
from ui.simulation.stochastic_simulation_dialog import StochasticSimulationDialog
from ui.species.species_tab import SpeciesTab


def show_error_message(message):
    error_message = QMessageBox()
    error_message.setIcon(QMessageBox.Warning)
    error_message.setWindowTitle("Error")
    error_message.setStandardButtons(QMessageBox.Ok)
    error_message.setText(message)

    button = error_message.exec_()
    if button == QMessageBox.Ok:
        return True
    else:
        return False


class GeneWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self._init_ui()

    def _init_ui(self):
        layout = QVBoxLayout()

        tabs = QTabWidget()

        self.species_tab = SpeciesTab()
        self.reactions_tab = ReactionsTab()
        self.rev_eng_modify_tab = ReverseEngineeringModifyTab()
        tabs.addTab(self.species_tab, "Species")
        tabs.addTab(self.reactions_tab, "Reactions")
        tabs.addTab(self.rev_eng_modify_tab, "Reverse Engineering")

        self._init_menubar()

        layout.addWidget(tabs)
        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)
        self.setWindowTitle("Gene")
        self.show()

    def _deterministic_simulation_clicked(self):
        def handler(s):
            net = GeneController.get_instance().network
            OdeSimulator.visualise(net, s, OdeSimulator.simulate(net, s))

        DeterministicSimulationDialog(handler)

    def _stochastic_simulation_clicked(self):
        StochasticSimulationDialog()

    def _save_file_as_sbml_clicked(self):
        net = GeneController.get_instance().network
        d = QFileDialog()
        filename = d.getSaveFileName(self, "Save file", ".", "XML Files (*.xml)")

        if filename:
            SbmlSaver.save_network_to_file(net, filename[0] + ".xml")

    def _init_menubar(self):
        self.menubar = self.menuBar()

        # File menu
        file = self.menubar.addMenu("File")

        open_file = QAction("Open SBML file", self)
        open_file.triggered.connect(lambda _: OpenSbmlDialog(self))
        file.addAction(open_file)

        save_file = QAction("Save as SBML file", self)
        save_file.triggered.connect(self._save_file_as_sbml_clicked)
        file.addAction(save_file)

        # Simulate menu
        simulate = self.menubar.addMenu("Simulate")

        deterministic = QAction("Deterministic (ODE)", self)
        deterministic.triggered.connect(self._deterministic_simulation_clicked)
        simulate.addAction(deterministic)

        stochastic = QAction("Stochastic (Gillespie Algorithm)", self)
        stochastic.triggered.connect(self._stochastic_simulation_clicked)
        simulate.addAction(stochastic)


app = QApplication([])
g = GeneWindow()
app.exec_()
